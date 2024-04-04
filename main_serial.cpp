#include <stdlib.h>

#include <iomanip>
#include <iostream>
#include <thread>
#include <vector>

#include "core/csv.h"
#include "core/utils.h"

using Edge = std::tuple<int, double, int>;

struct Graph {
	std::vector<std::vector<double>> adj_matrix;

	Graph(const std::string& e_file, const int num_nodes) {
		io::CSVReader<3> edges_file(e_file);

		Edge edge;

		adj_matrix.resize(num_nodes);
		for (auto& row : adj_matrix) {
			row.resize(num_nodes, 0);
		}

		while (edges_file.read_row(std::get<0>(edge), std::get<1>(edge), std::get<2>(edge))) {
			// std::cout << std::get<0>(edge) << " " << std::get<1>(edge) << " " << std::get<2>(edge) << std::endl;
			adj_matrix[std::get<0>(edge)][std::get<2>(edge)] = std::get<1>(edge);
		}
	}

	void print() const {
		std::cout << "Adjacency Matrix: " << std::endl;
		for (const auto& row : adj_matrix) {
			for (const auto& elem : row) {
				std::cout << elem << " ";
			}
			std::cout << std::endl;
		}
	}
};

int main(int argc, char* argv[]) {
	cxxopts::Options options(
		"main_serial", "Calculate All Pair Shortest Path using serial execution");
	options.add_options(
		"", {
				{"edgesFile", "Input graph file path", cxxopts::value<std::string>()->default_value("./input_graphs/10Edges.csv")},
				{"numNodes", "Input graph file path", cxxopts::value<int>()->default_value("11")},
			});

	std::cout << std::fixed;
	auto cl_options = options.parse(argc, argv);
	const std::string e_file = cl_options["edgesFile"].as<std::string>();
	const int num_nodes = cl_options["numNodes"].as<int>();

	Graph graph(e_file, num_nodes);

	// graph.print();

	return 0;
}
