#include <assert.h>
#include <stdlib.h>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <tuple>
#include <vector>
#include <mpi.h>
#include <math.h>

#include "../core/csv.h"
#include "../core/get_time.h"
#include "../core/utils.h"

#ifdef PRINT
const bool _print = true;
#else
const bool _print = false;
#endif

#define ROOT 0
#define MPI_TAG 1
#define TRUE 1
#define FALSE 0   
#define MIN(A, B) (A < B) ? A : B

typedef struct {
    int rank;
    int row, col;
    int p, q;
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
	Graph(const std::string& e_file, const int num_nodes) {
		this->e_file = e_file;
		this->row_size = num_nodes;
	}

	int get_row_size() {
		return this->row_size;
	}

	std::vector<double> get_matrix() {
		return this->sol;
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
void setup_grid(GRID_INFO *grid);
int check_fox(int p, int n);
void send_sub_mtrx(int *_mtrx, int n, int q);
void *process_mtrx(GRID_INFO *grid, double *time, int *mtrx_A, int n);
void floyd_warshall(int *A, int *B, int *C, int n);
void print_mtrx(int *mtrx, int row_size);

int main(int argc, char* argv[]) {
	std::cout << std::fixed << std::setprecision(2);

	cxxopts::Options options(
		"main_serial", "Calculate All Pair Shortest Path using serial execution");
	options.add_options(
		"", {
				{"file", "Input graph file path", cxxopts::value<std::string>()->default_value("./input_graphs/10Edges.csv")},
				{"nodes", "Number of Nodes in the file", cxxopts::value<int>()->default_value("10")},
			});

	auto cl_options = options.parse(argc, argv);
	const std::string e_file = cl_options["file"].as<std::string>();
	const int num_nodes = cl_options["nodes"].as<int>();

	std::cout << "Input Graph File: " << e_file << std::endl;
	std::cout << "Number of Nodes: " << num_nodes << std::endl;

	// set up MPI environment
    MPI_Init(&argc, &argv);
    GRID_INFO grid;
    setup_grid(&grid);	

	int n = num_nodes; 
	// double *mtrx[n] = new int [n * n];
	std::vector<double> mtrx;
	if (grid.rank == ROOT) {
		Graph graph(e_file, num_nodes);
		graph.read();
		std::cout << "Solution Matrix before Floyd Warshall" << std::endl;
		graph.print();

		if (!is_valid_fox(grid.p, n)) {
			fprintf(stderr, "Fox algorithm can't be applied with a matrix of size %d and %d processes.\nAborting...\n", n, grid.p);
            MPI_Abort(MPI_COMM_WORLD, 0);
            exit(1);
		}
		mtrx = graph.get_matrix();

	}

	if (grid.rank == ROOT && grid.p > 1) {
		send_sub_mtrx(mtrx, n, grid.q);
	}

	const int m = n / grid.q;
	double *mtrx_A = new double[m*m];
	if (grid.p > 1) {
		MPI_Recv(mtrx_A, m * m, MPI_INT, ROOT, MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	}
	else {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				mtrx_A[at(i,j,n)] = mtrx[at(i,j,n)];
			}
		}
	}

	double time;
	double *mtrx_C = process_mtrx(&grid, &time, mtrx_A, n);
	if (grid.p > 1) {
		delete[] mtrx_A;
	}

	double *mtrx_F = new double[n * n];
	MPI_Gather(mtrx_C, m*m, MPI_INT, mtrx_F, m*m, MPI_INT, ROOT, MPI_COMM_WORLD);
	delete[] mtrx_C; 

	// timer t;
	// t.start();
	// graph.floydWarshall();
	// auto end = t.stop();

	if (grid.rank == ROOT) {
		std::cout << "\nSolution Matrix after Floyd Warshall" << std::endl;
		print_final(mtrx_F, n);
		fprintf(stderr, "\nExecution Time: %10.3lf milliseconds.\n\n", time * 1000);
	}
    delete[] mtrx_F;
    MPI_Comm_free(&grid.comm);
    MPI_Comm_free(&grid.row_comm);
    MPI_Comm_free(&grid.col_comm);
    MPI_Finalize();
    return 0;
}

const int at(int i, int j, int n){
	// for speed purposes
	return i * n + j;
}

void setup_grid (GRID_INFO *grid) {
	MPI_Comm_size(MPI_COMM_WORLD, &(grid->p));
	MPI_Comm_rank(MPI_COMM_WORLD, &(grid->rank));

	grid->q = sqrt(grid->p);

	int dim[2] = {grid->q, grid->q};
	int periods[2] = {1, 1};
	MPI_Cart_create(MPI_COMM_WORLD, 2, dim, periods, 1, &(grid->comm));

	int coords[2];
	MPI_Comm_rank(grid->comm, &(grid->rank));
	MPI_Cart_coords(grid->comm, grid->rank, 2, coords);
	grid->row = coords[0];
	grid-> col = coords[1];

	coords[0] = 0;
	coords[1] = 1;
	MPI_Cart_sub(grid->comm, coords, &(grid->row_comm));

	coords[0] = 1;
	coords[1] = 0;
	MPI_Cart_sub(grid->comm, coords, &(grid->col_comm));
}

int is_valid_fox(int p, int n) {
	int q = sqrt(p);
	if (q*q == p && n % q == 0) {
		return TRUE;
	}
	return FALSE;
}

void send_sub_mtrx(std::vector<double> mtrx, int n, int q) {
	int m = n / q;
	double (*sub_mtrx) = new double[m*m];

	int dst = 0;
	for (int k = 0; k < q; k++) {
		for (int l = 0; l < q; l++) {
			for (int i = 0; i < m; i++) {
				for (int j = 0; j < m; j++) {
					sub_mtrx[at(i, j, m)] = mtrx[at(i + (k * m), j + (l * m), m)];
				}
			}
			MPI_Send(sub_mtrx, m * m, MPI_INT, dst++, MPI_TAG, MPI_COMM_WORLD);	// assign sub matrix for all processes from root to the last process
		}
	}
	delete[] sub_mtrx;
}

void *process_mtrx(GRID_INFO *grid, double *time, double *mtrx_A, int n) {
	int m = n / grid->q;
	int src = (grid->row + 1) % grid->q;
	int dst = (grid->row - 1 + grid->q) % grid->q;

	double *temp_A = new double[m * m];
	//assert(temp_A != nullptr);
	double *mtrx_B = new double[m * m];
	//assert(mtrx_B != nullptr);
	double *mtrx_C = new double[m * m];
	//assert(mtrx_C != nullptr);
	memcpy(mtrx_C, mtrx_A, m * m * sizeof(double));

	int iter;
	*time = MPI_Wtime();
	for (iter = 1; iter < n; iter <<= 1) {
		memcpy(mtrx_B, mtrx_C, m * m * sizeof(int));
		for (int stage = 0; stage < grid->q; stage++) {
			int bcast_root = (grid->row + stage) % grid->q;
			if (bcast_root == grid->col) {
				MPI_Bcast(mtrx_A, m * m, MPI_INT, bcast_root, grid->row_comm);
				floyd_warshall(mtrx_A, mtrx_B, mtrx_C, m);
			} else {
				MPI_Bcast(temp_A, m * m, MPI_INT, bcast_root, grid->row_comm);
				floyd_warshall(temp_A, mtrx_B, mtrx_C, m);
			}
			MPI_Sendrecv_replace(mtrx_B, m * m, MPI_INT, dst, MPI_TAG, src, MPI_TAG, grid->col_comm, MPI_STATUS_IGNORE);
		}
	}
	*time = MPI_Wtime() - *time;
	delete[] temp_A;
	delete[] mtrx_B;
	return mtrx_C;
}

inline void floyd_warshall(double *A, double *B, double *C, int n) {
	// double (*A)[n] = (double (*)[n]) _A;
	// double (*B)[n] = (double (*)[n]) _B;
	// double (*C)[n] = (double (*)[n]) _C;

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
				C[at(i,j,n)] = MIN(C[at(i,j,n)], A[at(i,k,n)] + B[at(k,j,n)]);
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
			const auto elem = mtrx_F[at(i,j,row_size)];
			if (elem > 1000000) {
				std::cout << "INF\t";
				continue;
			}
			std::cout << elem << "\t";
		}
		std::cout << std::endl;
	}

}
