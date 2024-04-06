#include <assert.h>
#include <stdlib.h>

#include <iomanip>
#include <iostream>
#include <thread>
#include <tuple>
#include <vector>

#include "core/csv.h"
#include "core/utils.h"

// 					  { src, weight, dest }
using Edge = std::tuple<int, double, int>;
using Arr = std::vector<double>;

const double INF = std::numeric_limits<double>::max();

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
				if (elem == INF) {
					std::cout << "INF\t";
					continue;
				}
				std::cout << elem << "\t";
			}
			std::cout << std::endl;
		}
	}

	void PL_APSP(int size, int portion, int start, int end) {
		for (int k = 0; k < size; ++k) {
			for (int i = start; i < end; ++i) {
				for (int j = 0; j < size; ++j) {
					if (sol[at(i, k)] != INF && sol[at(k, j)] != INF) {
						int sum = sol[at(i, k)] + sol[at(k, j)];
						if (sol[at(i, j)] > sum) {
							sol[at(i, j)] = sum;
						}
					}
				}
			}
		}
	}

	void computeShortestPaths() {
		std::vector<std::thread> threads;
		int num_threads = std::thread::hardware_concurrency();
		const auto N = sol.size();
		int portion = N / num_threads;

		for (int t = 0; t < num_threads; ++t) {
			int start = t * portion;
			int end = (t == num_threads - 1) ? N : (t + 1) * portion;
			threads.emplace_back(&Graph::PL_APSP, N, portion, start, end);
		}

		for (auto& t : threads) {
			t.join();
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
	graph.computeShortestPaths();
	std::cout << "\nSolution Matrix after Floyd Warshall" << std::endl;
	graph.print();

	return 0;
}
