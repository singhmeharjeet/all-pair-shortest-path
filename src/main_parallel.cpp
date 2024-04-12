#include <assert.h>
#include <stdlib.h>

#include <iomanip>
#include <iostream>
#include <mutex>
#include <thread>
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

const int NUM_THREADS = 4;
CustomBarrier barrier(NUM_THREADS);
const auto INF = std::numeric_limits<double>::max();
std::mutex mtx;
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

	void worker(int tid, int k, int i_start, int i_end, int j_start, int j_end) {
		for (int i = i_start; i < i_end; i++) {
			for (int j = j_start; j < j_end; j++) {
				sol[at(i, j)] = std::min(sol[at(i, j)], sol[at(i, k)] + sol[at(k, j)]);
			}
		}
	}

	void floydWarshall() {
		if (!validateInputs()) {
			return;
		}

		const int chunk_size = row_size / NUM_THREADS;

		std::vector<std::thread> threads(NUM_THREADS);
		for (int k = 0; k < row_size; k++) {
			for (int i = 0; i < NUM_THREADS; i++) {
				threads[i] = std::thread(
					&Graph::worker,
					this,
					i, k, i * chunk_size, (i + 1) * chunk_size, 0, row_size);
			}
			for (auto& t : threads) {
				t.join();
			}
		}
	}
};

int main(int argc, char* argv[]) {
	std::cout << std::fixed << std::setprecision(2);

	cxxopts::Options options(
		"main_parallel", "Calculate All Pair Shortest Path using parallel execution");
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
