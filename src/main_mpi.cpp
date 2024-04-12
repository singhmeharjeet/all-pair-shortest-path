#include <assert.h>
#include <stdlib.h>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <tuple>
#include <vector>
#include <mpi.h>

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
	std::vector<std::vector<double>> sol;	// 2D vector 
	std::string e_file;
	int row_size;

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

	void read() {
		if (!validateInputs()) {
			return;
		}

		io::CSVReader<3> edges_file(e_file);

		// Initialize the 2D vector
		for (int i = 0; i < row_size; i++) {
			std::vector<double> row(row_size, INF);
			sol.push_back(row);
		}
		// Set all path from any vertex i to itself to 0
		for (int i = 0; i < row_size; i++) {
			sol[i][i] = 0;
		}

		Edge edge;
		while (edges_file.read_row(std::get<0>(edge), std::get<1>(edge), std::get<2>(edge))) {
			sol[std::get<0>(edge)][std::get<2>(edge)] = std::get<1>(edge);
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
				const auto elem = sol[i][j];
				if (elem > 1000000) {
					std::cout << "INF\t";
					continue;
				}
				std::cout << elem << "\t";
			}
			std::cout << std::endl;
		}
	}

	void floydWarshall() {
		if (!validateInputs()) {
			return;
		}

		int i, j, k;

		for (k = 0; k < row_size; k++) {
			// Pick all vertices as source one by one
			for (i = 0; i < row_size; i++) {
				// Pick all vertices as destination for the
				// above picked source
				for (j = 0; j < row_size; j++) {
					// If vertex k is on the shortest path from
					// i to j, then update the value of
					// sol_matrix[i][j]
					const auto sum = sol[i][k] + sol[k][j];
					auto& val = sol[i][j];

					if (sum < val) {
						val = sum;
					}
				}
			}
		}
	}
};

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

	Graph graph(e_file, num_nodes);

	graph.read();

	if (_print) {
		std::cout << "Solution Matrix before Floyd Warshall" << std::endl;
		graph.print();
	}

	timer t;
	t.start();
	graph.floydWarshall();
	auto end = t.stop();

	if (_print) {
		std::cout << "\nSolution Matrix after Floyd Warshall" << std::endl;
		graph.print();
	}

	std::cout << "\nTime taken: " << std::setprecision(6) << end << " seconds\n";

	return 0;
}
