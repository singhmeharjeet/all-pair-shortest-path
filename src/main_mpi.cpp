#include <assert.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <tuple>
#include <vector>

#include "../core/csv.h"
#include "../core/get_time.h"
#include "../core/utils.h"

#ifdef PRINT
const bool _print = true;
#else
const bool _print = false;
#endif

#define ROOT	0
#define MPI_TAG 1
#define DONE	true

// struct to hold the necessary information of a process
typedef struct {
	int rank;											   // global rank of a process
	int row, col;										   // row and col positions of a process in the Cartesian topology
	int processes, processes_per_row, sub_mtrx_dimension;  // Number of processes in total, number of processes in a dimension of the Cartesian topology, dimension size of a submatrix assigned to a process
	MPI_Comm comm;										   // Cartesian topology communicator
	MPI_Comm row_comm;									   // Row communicator of the topology
	MPI_Comm col_comm;									   // Column communicaator of the topology
} GRID_INFO;

const auto INF = std::numeric_limits<double>::max();  // define infinite
using Edge = std::tuple<int, double, int>;

// Weighted directed graph with no negative-weighted edge
class Graph {
	std::vector<double> mtrx;
	std::string e_file;
	int mtrx_dimension;

	// To access a 1D array like 2D
	inline const int at(const int i, const int j) const {
		// for speed purposes
		return i * mtrx_dimension + j;
	}

	bool validateInputs() const {
		if (mtrx_dimension <= 0) {
			std::cerr << "Number of nodes should be greater than 0" << std::endl;
			return false;
		}
		if (e_file.empty()) {
			std::cerr << "Edges file path should not be empty" << std::endl;
			return false;
		}
		return true;
	}

	public:
	Graph(const std::string &e_file, const int num_nodes) {
		this->e_file = e_file;
		this->mtrx_dimension = num_nodes;
	}

	const int get_mtrx_dimension() const {
		return this->mtrx_dimension;
	}
	const std::string &get_file() const {
		return this->e_file;
	}
	std::vector<double> &get_matrix() {
		return (this->mtrx);
	}

	void read() {
		if (!validateInputs()) {
			return;
		}

		io::CSVReader<3> edges_file(e_file);

		mtrx.resize(mtrx_dimension * mtrx_dimension, INF);	// Initialize matrix with INF

		// Set diagonal values to 0
		for (int i = 0; i < mtrx_dimension; i++) {
			mtrx[at(i, i)] = 0;
		}

		// Read in the rest of the values from the csv file
		Edge edge;
		while (edges_file.read_row(std::get<0>(edge), std::get<1>(edge), std::get<2>(edge))) {
			mtrx[at(std::get<0>(edge), std::get<2>(edge))] = std::get<1>(edge);
		}
	}

	void print() const {
		if (mtrx.empty()) {
			std::cerr << "Inital matrix is empty, Please Read first" << std::endl;
			return;
		}

		if (!validateInputs()) {
			return;
		}

		std::cout << "\t";
		for (int x = 0; x < mtrx_dimension; x++) {
			std::cout << x << "\t";
		}
		std::cout << std::endl;

		for (int i = 0; i < mtrx_dimension; i++) {
			std::cout << i << "\t";
			for (int j = 0; j < mtrx_dimension; j++) {
				const auto elem = mtrx[at(i, j)];
				if (elem > 1000000) {
					std::cout << "INF\t";
					continue;
				}
				std::cout << elem << "\t";
			}
			std::cout << std::endl;
		}
	}
};

const int at(int i, int j, int n);
void initialize(int argc, char **argv, GRID_INFO *grid, int num_nodes);
void deinitialize(GRID_INFO *grid);
bool is_valid_fox(int p, int n);
double *process_mtrx(const GRID_INFO &grid, double *time, double *mtrx_A, int n);
void floyd_warshall(double *A, double *B, double *C, int n);
inline void floyd_warshall(int mtrx_dimension, std::vector<double> &mtrx);
void print_mtrx(double *mtrx, int mtrx_dimension);
void read_args(int argc, char **argv, std::string &e_file, int &num_nodes);
void read_csv(const GRID_INFO &grid, Graph &graph);
bool send_processes(GRID_INFO &grid, Graph &graph, double *mtrx_A = nullptr);
void fix_final_mtrx(double *mtrx_to_fix, double *mtrx_F, int mtrx_dimension, int processes_per_row, int sub_mtrx_dimensio);

int main(int argc, char *argv[]) {
	std::cout << std::fixed << std::setprecision(2);

	std::string e_file;
	int num_nodes;
	read_args(argc, argv, e_file, num_nodes);

	GRID_INFO grid;
	Graph graph(e_file, num_nodes);

	initialize(argc, argv, &grid, num_nodes);
	read_csv(grid, graph);

	auto sub_mtrx_dimension = grid.sub_mtrx_dimension;
	double mtrx_A[sub_mtrx_dimension * sub_mtrx_dimension];

	// Root sends submatrics to other processes
	bool wasSerial;
	if (grid.rank == ROOT) {
		wasSerial = send_processes(grid, graph, mtrx_A);
	}
	if (wasSerial) {  // if there's only 1 process
		return 0;
	}

	// Receive submatrix sent by root
	if (grid.rank != ROOT) {
		MPI_Recv(
			mtrx_A, sub_mtrx_dimension * sub_mtrx_dimension,
			MPI_DOUBLE,
			0,
			MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	double time;
	double *mtrx_C = process_mtrx(grid, &time, mtrx_A, graph.get_mtrx_dimension());	 // process assigned submatrix

	double *mtrx_F = new double[graph.get_mtrx_dimension() * graph.get_mtrx_dimension()];
	double *mtrx_gather = new double[graph.get_mtrx_dimension() * graph.get_mtrx_dimension()];

	// Gather all processed submatrices from all processes
	MPI_Gather(mtrx_C, sub_mtrx_dimension * sub_mtrx_dimension, MPI_DOUBLE, mtrx_gather, sub_mtrx_dimension * sub_mtrx_dimension, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

	// Gather version using send/receive
	// if (grid.rank == ROOT) {
	// 	// fill process 0 submatrix
	// 	for (int i = 0; i < sub_mtrx_dimension; i++){
	// 		for (int j = 0; j < sub_mtrx_dimension; j++){
	// 			mtrx_F[at(i, j, graph.get_mtrx_dimension())] = mtrx_C[at(i, j, sub_mtrx_dimension)];
	// 		}
	// 	}

	// 	for (int i = 1; i < grid.processes; i++){

	// 		MPI_Recv(
	// 			mtrx_C, sub_mtrx_dimension * sub_mtrx_dimension,
	// 			MPI_DOUBLE,
	// 			i,
	// 			MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	// 		// Move down in rows only if the process surpasses processes_per row.
	// 		int row_start = (i / grid.processes_per_row) * sub_mtrx_dimension;
	// 		// Move to the right in cols only if the process hasn't supasses processes per_row
	// 		int col_start = (i % grid.processes_per_row) * sub_mtrx_dimension;

	// 		for (int i = 0; i < sub_mtrx_dimension; i++){
	// 			for (int j = 0; j < sub_mtrx_dimension; j++){
	// 				mtrx_F[at(i + row_start, j + col_start, graph.get_mtrx_dimension())] = mtrx_C[at(i, j, sub_mtrx_dimension)];
	// 			}
	// 		}
	// 	}
	// } else {
	// 	MPI_Send(
	// 					mtrx_C, sub_mtrx_dimension * sub_mtrx_dimension,
	// 					MPI_DOUBLE,
	// 					0,
	// 					MPI_TAG, MPI_COMM_WORLD);
	// }

	delete[] mtrx_C;

	// Print the final results
	if (grid.rank == ROOT) {
		std::cout << "\nSolution Matrix after Floyd Warshall" << std::endl;
		fix_final_mtrx(mtrx_gather, mtrx_F, graph.get_mtrx_dimension(), grid.processes_per_row, grid.sub_mtrx_dimension);
		// print_mtrx(mtrx_F, graph.get_mtrx_dimension());
		std::cout << "\nTime taken: " << std::setprecision(6) << time << " seconds\n";
	}

	delete[] mtrx_F;

	deinitialize(&grid);
	return 0;
}

// Treat 1D array like 2D
const int at(int i, int j, int n) {
	// for speed purposes
	return i * n + j;
}

// Initialize MPI environment
void initialize(int argc, char **argv, GRID_INFO *grid, int n) {
	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &(grid->processes));
	MPI_Comm_rank(MPI_COMM_WORLD, &(grid->rank));

	if (!is_valid_fox(grid->processes, n) && grid->rank == ROOT) {
		fprintf(stderr, "Fox algorithm can't be applied with a matrix of size %d and %d processes.\nAborting...\n", n, grid->processes);
		MPI_Abort(MPI_COMM_WORLD, 1);
		exit(1);
	}

	grid->processes_per_row = sqrt(grid->processes);
	grid->sub_mtrx_dimension = n / grid->processes_per_row;

	// Create a topology for the processes
	int dim[2] = {grid->processes_per_row, grid->processes_per_row};
	int periods[2] = {1, 1};
	MPI_Cart_create(MPI_COMM_WORLD, 2, dim, periods, 1, &(grid->comm));

	// Get the row and column position of the process in the topology
	int coords[2];
	MPI_Comm_rank(grid->comm, &(grid->rank));
	MPI_Cart_coords(grid->comm, grid->rank, 2, coords);
	grid->row = coords[0];
	grid->col = coords[1];

	// Create row communciator for processes in the same row to communicate
	coords[0] = 0;
	coords[1] = 1;
	MPI_Cart_sub(grid->comm, coords, &(grid->row_comm));

	// Create column communicator for proceseses in the same column to communicate
	coords[0] = 1;
	coords[1] = 0;
	MPI_Cart_sub(grid->comm, coords, &(grid->col_comm));
}

// Close the MPI environment
void deinitialize(GRID_INFO *grid) {
	MPI_Comm_free(&(grid->comm));
	MPI_Comm_free(&(grid->row_comm));
	MPI_Comm_free(&(grid->col_comm));
	MPI_Finalize();
}

// Check if fox can be applied on the matrix with given parameters
bool is_valid_fox(int p, int n) {
	// p = number of processes
	// n = matrix dimension
	// q = number of processors per row
	int q = sqrt(p);
	if (q * q == p && n % q == 0) {
		return true;
	}
	return false;
}

// Process the submatrix assigned to each process
double *process_mtrx(const GRID_INFO &grid, double *time, double *mtrx_A, int n) {
	const int src = (grid.row + 1) % grid.processes_per_row;
	const int dst = (grid.row - 1 + grid.processes_per_row) % grid.processes_per_row;

	const auto size_sqr = grid.sub_mtrx_dimension * grid.sub_mtrx_dimension;
	const auto max_iters = n < 100 ? n : 100;

	double temp_A[size_sqr];
	double mtrx_B[size_sqr];
	double *mtrx_C = new double[size_sqr];

	memcpy(mtrx_C, mtrx_A, size_sqr * sizeof(double));

	*time = MPI_Wtime();
	for (auto iter = 1; iter < max_iters; iter += 1) {	// Iterations needed to converge to the right answer
		memcpy(mtrx_B, mtrx_C, size_sqr * sizeof(double));

		for (int stage = 0; stage < grid.processes_per_row; stage++) {
			int bcast_root = (grid.row + stage) % grid.processes_per_row;
			if (bcast_root == grid.col) {  // If you're the root in this current stage, broadcast your submatrix to other processses in the same row and call Floyd-Warshall
				MPI_Bcast(mtrx_A, size_sqr, MPI_DOUBLE, bcast_root, grid.row_comm);
				floyd_warshall(mtrx_A, mtrx_B, mtrx_C, grid.sub_mtrx_dimension);
			} else {  // Otherwise, receive the broadcasted submatrix and call Floyd-Warshall
				MPI_Bcast(temp_A, size_sqr, MPI_DOUBLE, bcast_root, grid.row_comm);
				floyd_warshall(temp_A, mtrx_B, mtrx_C, grid.sub_mtrx_dimension);
			}
			MPI_Sendrecv_replace(mtrx_B, size_sqr, MPI_DOUBLE, dst, MPI_TAG, src, MPI_TAG, grid.col_comm, MPI_STATUS_IGNORE);  // Send your submatrix to the process above while receiving the submatrix from the process below
		}
	}
	*time = MPI_Wtime() - *time;
	return mtrx_C;
}

// Floyd Warshall version to be used for cases with only one process
inline void floyd_warshall(int mtrx_dimension, std::vector<double> &mtrx) {
	int i, j, k;

	for (k = 0; k < mtrx_dimension; k++) {
		for (i = 0; i < mtrx_dimension; i++) {
			// above picked source
			for (j = 0; j < mtrx_dimension; j++) {
				const auto sum = mtrx[at(i, k, mtrx_dimension)] + mtrx[at(k, j, mtrx_dimension)];
				auto &val = mtrx[at(i, j, mtrx_dimension)];

				if (sum < val) {
					val = sum;
				}
			}
		}
	}
}

// Floyd Warshall to be used for cases with multiple processes
inline void floyd_warshall(double *A, double *B, double *C, int n) {
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
				C[at(i, j, n)] = std::min(C[at(i, j, n)], A[at(i, k, n)] + B[at(k, j, n)]);
}

void print_mtrx(double *mtrx, int mtrx_dimension) {
	std::cout << "\t";
	for (int x = 0; x < mtrx_dimension; x++) {
		std::cout << x << "\t";
	}
	std::cout << std::endl;

	for (int i = 0; i < mtrx_dimension; i++) {
		std::cout << i << "\t";
		for (int j = 0; j < mtrx_dimension; j++) {
			const auto elem = mtrx[at(i, j, mtrx_dimension)];
			if (elem > 1000000) {
				std::cout << "INF\t";
				continue;
			}
			std::cout << elem << "\t";
		}
		std::cout << std::endl;
	}
}

// Read user inputs
void read_args(int argc, char **argv, std::string &e_file, int &num_nodes) {
	cxxopts::Options options(
		"main_mpi", "Calculate All Pair Shortest Path using MPI execution");
	options.add_options(
		"", {
				{"file", "Input graph file path", cxxopts::value<std::string>()->default_value("./input_graphs/10Edges.csv")},
				{"nodes", "Number of Nodes in the file", cxxopts::value<int>()->default_value("10")},
			});

	auto cl_options = options.parse(argc, argv);
	e_file = cl_options["file"].as<std::string>();
	num_nodes = cl_options["nodes"].as<int>();
}

// Read in the graph and display information
void read_csv(const GRID_INFO &grid, Graph &graph) {
	if (grid.rank == ROOT) {
		std::cout << "Input Graph File: " << graph.get_file() << std::endl
				  << "Number of Nodes: " << graph.get_mtrx_dimension() << std::endl
				  << "Processes: " << grid.processes << std::endl
				  << "Processes/row: " << grid.processes_per_row << std::endl
				  << std::endl;

		std::cout << "Inital matrix before Floyd Warshall" << std::endl;

		graph.read();
		// graph.print();
	}
}

// Split the matrix and send submatrices to other proccesses
bool send_processes(GRID_INFO &grid, Graph &graph, double *mtrx_A) {
	// Early exit if there is only one process (ROOT)
	if (grid.processes == 1) {
		timer t;
		t.start();
		floyd_warshall(graph.get_mtrx_dimension(), graph.get_matrix());	 // serial floyd warshall
		auto time = t.stop();

		std::cout << "\nSolution Matrix after Floyd Warshall for 1 process:" << std::endl;
		graph.print();
		std::cout << "\nTime taken: " << std::setprecision(6) << time << " seconds\n";

		deinitialize(&grid);
		return DONE;  // was serial
	}

	auto processes_per_row = grid.processes_per_row;
	auto sub_mtrx_dimension = grid.sub_mtrx_dimension;

	// Break the Adj Matrix into submatrices and send them to all processes
	double sub_mtrx[sub_mtrx_dimension * sub_mtrx_dimension];

	for (int sub_mtrx_row = 0; sub_mtrx_row < processes_per_row; sub_mtrx_row++) {	// Break and send to each process in the topology
		for (int sub_mtrx_col = 0; sub_mtrx_col < processes_per_row; sub_mtrx_col++) {
			for (int i = 0; i < sub_mtrx_dimension; i++) {	// Submatrix
				for (int j = 0; j < sub_mtrx_dimension; j++) {
					sub_mtrx[at(i, j, sub_mtrx_dimension)] =
						graph.get_matrix()[at((sub_mtrx_row * sub_mtrx_dimension) + i, (sub_mtrx_col * sub_mtrx_dimension) + j, processes_per_row * sub_mtrx_dimension)];
				}
			}
			// Sending & Receiving the processes is checked
			if (sub_mtrx_row * processes_per_row + sub_mtrx_col == ROOT) {
				memcpy(mtrx_A, sub_mtrx, sub_mtrx_dimension * sub_mtrx_dimension * sizeof(double));
			} else {
				MPI_Send(sub_mtrx, sub_mtrx_dimension * sub_mtrx_dimension, MPI_DOUBLE, (sub_mtrx_row * processes_per_row) + sub_mtrx_col, MPI_TAG, MPI_COMM_WORLD);
			}
		}
	}

	return !DONE;  // have many processes
}

// Reorder the submatrices recived from other processes in the final matrix
void fix_final_mtrx(double *mtrx_to_fix, double *mtrx_F, int mtrx_dimension, int processes_per_row, int sub_mtrx_dimension) {
	int i, j, k, l;
	int a = 0, b = 0;
	int q = processes_per_row;
	int m = sub_mtrx_dimension;
	int n = mtrx_dimension;

	for (k = 0; k < q; k++) {
		for (l = 0; l < q; l++) {
			for (i = k * m; i < (k + 1) * m; i++) {
				for (j = l * m; j < (l + 1) * m; j++) {
					mtrx_F[at(i, j, n)] = mtrx_to_fix[at(a, b, n)] == INF ? 0 : mtrx_to_fix[at(a, b, n)];
					if (++b == n) {
						++a;
						b = 0;
					}
				}
			}
		}
	}
}