#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "../core/cxxopts.h"

void make(const int number, const std::string& outputPath) {
	// make a csv file to read the {src: int, weight: double, dest: int} edges
	// It should be a double linked graph circular

	std::ofstream file;
	try {
		file.open(outputPath);
	} catch (const std::exception& e) {
		std::cerr << "Error: " << e.what() << std::endl;
		return;
	}

	// the weight is random between 1 and 2

	for (int i = 0; i < number; i++) {
		file << i << "," << (double)rand() / RAND_MAX + 1 << "," << (i + 1) % number << std::endl;

		file << (i + 1) % number << "," << (double)rand() / RAND_MAX + 1 << "," << i << std::endl;
	}

	file.close();
}

int main(int argc, char** argv) {
	srand(42);
	cxxopts::Options options(
		"datamaker", "Make Data");
	options.add_options(
		"", {
				{"nodes", "Size of Nodes", cxxopts::value<int>()->default_value("10")},
			});

	auto cl_options = options.parse(argc, argv);
	const int nodes = cl_options["nodes"].as<int>();

	const std::string path = "./input_graphs/" + std::to_string(nodes) + "Edges.csv";

	make(nodes, path);
	std::cout << "Data made" << std::endl;
}