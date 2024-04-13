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

typedef struct {
	int rank;
	int row, col;
	int processes, chunks_per_row, chunk_size;
	MPI_Comm comm;
	MPI_Comm row_comm;
	MPI_Comm col_comm;
} GRID_INFO;

const auto INF = std::numeric_limits<double>::max();
using Edge = std::tuple<int, double, int>;

class Graph {
	std::vector<double> sol;
	std::string e_file;
	int row_size;

	inline const int at(const int i, const int j) const {
		// for speed purposes
		return i * row_size + j;
	}

	bool validateInputs() const {
		if (row_size <= 0) {
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
		this->row_size = num_nodes;
	}

	const int get_row_size() const {
		return this->row_size;
	}
	const std::string &get_file() const {
		return this->e_file;
	}
	std::vector<double> &get_matrix() {
		return (this->sol);
	}

	void read() {
		if (!validateInputs()) {
			return;
		}

		io::CSVReader<3> edges_file(e_file);

		sol.resize(row_size * row_size, INF);

		for (int i = 0; i < row_size; i++) {
			sol[at(i, i)] = 0;
		}

		Edge edge;
		while (edges_file.read_row(std::get<0>(edge), std::get<1>(edge), std::get<2>(edge))) {
			sol[at(std::get<0>(edge), std::get<2>(edge))] = std::get<1>(edge);
		}
	}

	void print() const {
		if (sol.empty()) {
			std::cerr << "Solution matrix is empty, Please Read first" << std::endl;
			return;
		}

		if (!validateInputs()) {
			return;
		}

		std::cout << "\t";
		for (int x = 0; x < row_size; x++) {
			std::cout << x << "\t";
		}
		std::cout << std::endl;

		for (int i = 0; i < row_size; i++) {
			std::cout << i << "\t";
			for (int j = 0; j < row_size; j++) {
				const auto elem = sol[at(i, j)];
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
inline void floyd_warshall(int row_size, std::vector<double> &sol);
void print_mtrx(double *mtrx, int row_size);
void read_args(int argc, char **argv, std::string &e_file, int &num_nodes);
void read_csv(const GRID_INFO &grid, Graph &graph);
bool send_chunks(GRID_INFO &grid, Graph &graph, double *mtrx_A = nullptr);

int main(int argc, char *argv[]) {
	std::cout << std::fixed << std::setprecision(2);

	std::string e_file;
	int num_nodes;
	read_args(argc, argv, e_file, num_nodes);

	GRID_INFO grid;
	Graph graph(e_file, num_nodes);

	initialize(argc, argv, &grid, num_nodes);
	read_csv(grid, graph);

	auto chunk_size = grid.chunk_size;
	double mtrx_A[chunk_size * chunk_size];

	auto wasSerial = send_chunks(grid, graph, mtrx_A);
	if (wasSerial) {
		return 0;
	}

	if (grid.rank != ROOT) {
		MPI_Recv(
			mtrx_A, chunk_size * chunk_size,
			MPI_DOUBLE,
			0,
			MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	if (grid.rank == ROOT) {
		std::cout << "ROOT received the chunks\n";
		print_mtrx(mtrx_A, chunk_size);
	} else if (grid.rank == 1) {
		std::cout << "1 received the chunks\n";
		print_mtrx(mtrx_A, chunk_size);
	} else if (grid.rank == 2) {
		std::cout << "2 received the chunks\n";
		print_mtrx(mtrx_A, chunk_size);
	} else if (grid.rank == 3) {
		std::cout << "3 received the chunks\n";
		print_mtrx(mtrx_A, chunk_size);
	}

	double time;
	double *mtrx_C = process_mtrx(grid, &time, mtrx_A, graph.get_row_size());

	double *mtrx_F = new double[graph.get_row_size() * graph.get_row_size()];
	MPI_Gather(mtrx_C, chunk_size * chunk_size, MPI_INT, mtrx_F, chunk_size * chunk_size, MPI_INT, ROOT, MPI_COMM_WORLD);
	delete[] mtrx_C;

	// timer t;
	// t.start();
	// graph.floydWarshall();
	// auto end = t.stop();

	if (grid.rank == ROOT) {
		std::cout << "\nSolution Matrix after Floyd Warshall" << std::endl;
		print_mtrx(mtrx_F, graph.get_row_size());
		fprintf(stderr, "\nExecution Time: %10.3lf milliseconds.\n\n", time * 1000);
	}

	delete[] mtrx_F;

	deinitialize(&grid);
	return 0;
}

const int at(int i, int j, int n) {
	// for speed purposes
	return i * n + j;
}
void initialize(int argc, char **argv, GRID_INFO *grid, int n) {
	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &(grid->processes));
	MPI_Comm_rank(MPI_COMM_WORLD, &(grid->rank));

	if (!is_valid_fox(grid->processes, n) && grid->rank == ROOT) {
		fprintf(stderr, "Fox algorithm can't be applied with a matrix of size %d and %d processes.\nAborting...\n", n, grid->processes);
		MPI_Abort(MPI_COMM_WORLD, 1);
		exit(1);
	}

	grid->chunks_per_row = sqrt(grid->processes);
	grid->chunk_size = n / grid->chunks_per_row;

	int dim[2] = {grid->chunks_per_row, grid->chunks_per_row};
	int periods[2] = {1, 1};
	MPI_Cart_create(MPI_COMM_WORLD, 2, dim, periods, 1, &(grid->comm));

	int coords[2];
	MPI_Comm_rank(grid->comm, &(grid->rank));
	MPI_Cart_coords(grid->comm, grid->rank, 2, coords);
	grid->row = coords[0];
	grid->col = coords[1];

	coords[0] = 0;
	coords[1] = 1;
	MPI_Cart_sub(grid->comm, coords, &(grid->row_comm));

	coords[0] = 1;
	coords[1] = 0;
	MPI_Cart_sub(grid->comm, coords, &(grid->col_comm));
}
void deinitialize(GRID_INFO *grid) {
	MPI_Comm_free(&(grid->comm));
	MPI_Comm_free(&(grid->row_comm));
	MPI_Comm_free(&(grid->col_comm));
	MPI_Finalize();
}
bool is_valid_fox(int p, int n) {
	// p = number of processes
	// n = row size of the matrix
	// q = number of processors per row
	int q = sqrt(p);
	if (q * q == p && n % q == 0) {
		return true;
	}
	return false;
}

double *process_mtrx(const GRID_INFO &grid, double *time, double *mtrx_A, int n) {
	const int src = (grid.row + 1) % grid.chunks_per_row;
	const int dst = (grid.row - 1 + grid.chunks_per_row) % grid.chunks_per_row;

	const auto size_sqr = grid.chunk_size * grid.chunk_size;

	double temp_A[size_sqr];
	double mtrx_B[size_sqr];
	double *mtrx_C = new double[size_sqr];

	// mtrx_A has the chunk of the matrix
	memcpy(mtrx_C, mtrx_A, size_sqr * sizeof(double));

	*time = MPI_Wtime();
	for (auto iter = 1; iter < n; iter <<= 1) {
		memcpy(mtrx_B, mtrx_C, size_sqr * sizeof(double));

		for (int stage = 0; stage < grid.chunks_per_row; stage++) {
			int bcast_root = (grid.row + stage) % grid.chunks_per_row;
			if (bcast_root == grid.col) {
				MPI_Bcast(mtrx_A, size_sqr, MPI_DOUBLE, bcast_root, grid.row_comm);
				floyd_warshall(mtrx_A, mtrx_B, mtrx_C, grid.chunk_size);
			} else {
				MPI_Bcast(temp_A, size_sqr, MPI_DOUBLE, bcast_root, grid.row_comm);
				floyd_warshall(temp_A, mtrx_B, mtrx_C, grid.chunk_size);
			}
			MPI_Sendrecv_replace(mtrx_B, size_sqr, MPI_DOUBLE, dst, MPI_TAG, src, MPI_TAG, grid.col_comm, MPI_STATUS_IGNORE);
		}
	}
	*time = MPI_Wtime() - *time;
	return mtrx_C;
}

inline void floyd_warshall(int row_size, std::vector<double> &sol) {
	int i, j, k;

	for (k = 0; k < row_size; k++) {
		for (i = 0; i < row_size; i++) {
			// above picked source
			for (j = 0; j < row_size; j++) {
				const auto sum = sol[at(i, k, row_size)] + sol[at(k, j, row_size)];
				auto &val = sol[at(i, j, row_size)];

				if (sum < val) {
					val = sum;
				}
			}
		}
	}
}
inline void floyd_warshall(double *A, double *B, double *C, int n) {
	// double (*A)[n] = (double (*)[n]) _A;
	// double (*B)[n] = (double (*)[n]) _B;
	// double (*C)[n] = (double (*)[n]) _C;

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
				C[at(i, j, n)] = std::min(C[at(i, j, n)], A[at(i, k, n)] + B[at(k, j, n)]);
}
void print_mtrx(double *mtrx, int row_size) {
	std::cout << "\t";
	for (int x = 0; x < row_size; x++) {
		std::cout << x << "\t";
	}
	std::cout << std::endl;

	for (int i = 0; i < row_size; i++) {
		std::cout << i << "\t";
		for (int j = 0; j < row_size; j++) {
			const auto elem = mtrx[at(i, j, row_size)];
			if (elem > 1000000) {
				std::cout << "INF\t";
				continue;
			}
			std::cout << elem << "\t";
		}
		std::cout << std::endl;
	}
}
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
void read_csv(const GRID_INFO &grid, Graph &graph) {
	if (grid.rank == ROOT) {
		std::cout << "Input Graph File: " << graph.get_file() << std::endl
				  << "Number of Nodes: " << graph.get_row_size() << std::endl
				  << "Chunks: " << grid.processes << std::endl
				  << "Chunk/row: " << grid.chunks_per_row << std::endl
				  << std::endl;

		std::cout << "Solution Matrix before Floyd Warshall" << std::endl;

		graph.read();
		graph.print();
	}
}
bool send_chunks(GRID_INFO &grid, Graph &graph, double *mtrx_A) {
	if (grid.rank == ROOT) {
		// Early exit if there is only one process (ROOT)
		if (grid.processes == 1) {
			timer t;
			t.start();
			floyd_warshall(graph.get_row_size(), graph.get_matrix());  // serial floyd warshall
			auto time = t.stop();

			std::cout << "\nSolution Matrix after Floyd Warshall" << std::endl;
			graph.print();
			std::cout << "\nTime taken: " << std::setprecision(6) << time << " seconds\n";

			deinitialize(&grid);
			return DONE;  // was serial
		}

		auto chunks_per_row = grid.chunks_per_row;
		auto chunk_size = grid.chunk_size;

		// Break the Adj Matrix into chucks(blocks) and send them to all processes
		double sub_mtrx[chunk_size * chunk_size];

		for (int chunk_row = 0; chunk_row < chunks_per_row; chunk_row++) {
			for (int chunk_col = 0; chunk_col < chunks_per_row; chunk_col++) {
				for (int i = 0; i < chunk_size; i++) {
					for (int j = 0; j < chunk_size; j++) {
						sub_mtrx[at(i, j, chunk_size)] = graph.get_matrix()[at(
							(chunk_row * chunk_size) + i,
							(chunk_col * chunk_size) + j,
							chunks_per_row * chunk_size)];
					}
				}
				// Sending & Receiving the chunks is checked; Dont change now
				// It will send to root also to save the complexity of handling root separately
				if (chunk_row * chunks_per_row + chunk_col == ROOT) {
					memcpy(mtrx_A, sub_mtrx, chunk_size * chunk_size * sizeof(double));
				} else {
					MPI_Send(
						sub_mtrx, chunk_size * chunk_size,
						MPI_DOUBLE,
						(chunk_row * chunks_per_row) + chunk_col,
						MPI_TAG, MPI_COMM_WORLD);
				}
			}
		}
	}
	return !DONE;  // was parallel
}