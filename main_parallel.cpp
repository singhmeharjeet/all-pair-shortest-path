#include <stdlib.h>

#include <iomanip>
#include <iostream>
#include <thread>
#include <vector>

#include "core/csv.h"
#include "core/utils.h"

int main(int argc, char *argv[]) {
	cxxopts::Options options(
		"main_serial", "Calculate All Pair Shortest Path using serial execution");
	options.add_options(
		"", {
				{"edgesFile", "Input graph file path", cxxopts::value<std::string>()->default_value("./input_graphs/10Edges.csv")},
				{"nodesFile", "Input graph file path", cxxopts::value<std::string>()->default_value("./input_graphs/10Nodes.csv")},
			});

	std::cout << std::fixed;
	auto cl_options = options.parse(argc, argv);
	const std::string e_file = cl_options["edgesFile"].as<std::string>();
	const std::string n_file = cl_options["nodesFile"].as<std::string>();

	io::CSVReader<3> edges_file(e_file);
	io::CSVReader<3> nodes_file(n_file);

	edges_file.read_header(io::ignore_no_column, "Edges/src", "Edges/w", "Edges/dest");
	nodes_file.read_header(io::ignore_extra_column, "Pos/x", "Pos/y", "Node/id");

	// edge = { src, w, dest }
	std::tuple<int, double, int> edge;
	// node = { x, y, id }
	std::tuple<double, double, int> node;

	auto cols = edges_file.get_column_names();
	for (int i = 0; i < 3; i++) {
		std::cout << cols[i] << " ";
	}
	std::cout << std::endl;

	while (edges_file.read_row(std::get<0>(edge), std::get<1>(edge), std::get<2>(edge))) {
		std::cout << std::get<0>(edge) << " " << std::get<1>(edge) << " " << std::get<2>(edge) << std::endl;
	}

	cols = nodes_file.get_column_names();
	for (int i = 0; i < 3; i++) {
		std::cout << cols[i] << " ";
	}
	std::cout << std::endl;
	while (nodes_file.read_row(std::get<0>(node), std::get<1>(node), std::get<2>(node))) {
		std::cout << std::get<0>(node) << " " << std::get<1>(node) << " " << std::get<2>(node) << std::endl;
	}

	return 0;
}
