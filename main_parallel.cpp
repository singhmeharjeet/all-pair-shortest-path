#include <assert.h>
#include <stdlib.h>

#include <iomanip>
#include <iostream>
#include <thread>
#include <tuple>
#include <vector>

#include "core/csv.h"
#include "core/utils.h"

using Edge = std::tuple<int, double, int>;

class Graph {
	std::vector<std::vector<double>> sol;
	std::string e_file;
	int num_nodes;

	bool validateInputs() const {
		if (num_nodes <= 0) {
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
		this->num_nodes = num_nodes;
	}

	void read() {
		if (!validateInputs()) {
			return;
		}

		io::CSVReader<3> edges_file(e_file);

		sol.resize(num_nodes);
		for (int i = 0; i < num_nodes; i++) {
			sol[i].resize(num_nodes, std::numeric_limits<double>::max());
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
		for (int x = 0; x < sol.size(); x++) {
			std::cout << x << "\t";
		}
		std::cout << std::endl;

		auto y = 0;
		for (const auto& row : sol) {
			std::cout << y++ << "\t";
			for (const auto& elem : row) {
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

		int V = sol.size();

		int i, j, k;

		for (k = 0; k < V; k++) {
			// Pick all vertices as source one by one
			for (i = 0; i < V; i++) {
				// Pick all vertices as destination for the
				// above picked source
				for (j = 0; j < V; j++) {
					// If vertex k is on the shortest path from
					// i to j, then update the value of
					// sol_matrix[i][j]
					if (sol[i][k] + sol[k][j] < sol[i][j])
						sol[i][j] = sol[i][k] + sol[k][j];
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
				{"edgesFile", "Input graph file path", cxxopts::value<std::string>()->default_value("./input_graphs/10Edges.csv")},
				{"numNodes", "Input graph file path", cxxopts::value<int>()->default_value("11")},
			});

	auto cl_options = options.parse(argc, argv);
	const std::string e_file = cl_options["edgesFile"].as<std::string>();
	const int num_nodes = cl_options["numNodes"].as<int>();

	Graph graph(e_file, num_nodes);

	graph.read();
	std::cout << "Solution Matrix before Floyd Warshall" << std::endl;
	graph.print();
	graph.floydWarshall();
	std::cout << "\nSolution Matrix after Floyd Warshall" << std::endl;
	graph.print();

	return 0;
}
