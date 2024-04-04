#include <stdlib.h>

#include <iomanip>
#include <iostream>
#include <thread>
#include <vector>

#include "core/csv.h"
#include "core/utils.h"

using Edge = std::tuple<int, double, int>;

typedef enum {
	ADJ,
	SOL
} TYPE;
struct Graph {
	std::vector<std::vector<double>> adj_matrix;
	std::vector<std::vector<double>> sol_matrix;

	Graph(const std::string& e_file, const int num_nodes) {
		io::CSVReader<3> edges_file(e_file);

		Edge edge;

		adj_matrix.resize(num_nodes);
		for (auto& row : adj_matrix) {
			row.resize(num_nodes, std::numeric_limits<double>::max());
		}
		for (int i = 0; i < num_nodes; i++) {
			adj_matrix[i][i] = 0;
		}

		while (edges_file.read_row(std::get<0>(edge), std::get<1>(edge), std::get<2>(edge))) {
			adj_matrix[std::get<0>(edge)][std::get<2>(edge)] = std::get<1>(edge);
		}
	}

	void print(TYPE t = ADJ) const {
		auto y = 0;
		if (t == ADJ) {
			std::cout << "Adjacency Matrix: " << std::endl;
			std::cout << "\t";
			for (int x = 0; x < adj_matrix.size(); x++) {
				std::cout << x << "\t";
			}
			std::cout << std::endl;
			for (const auto& row : adj_matrix) {
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
		} else {
			std::cout << "Solution Matrix: " << std::endl;
			std::cout << "\t";
			for (int x = 0; x < adj_matrix.size(); x++) {
				std::cout << x << "\t";
			}
			std::cout << std::endl;

			for (const auto& row : sol_matrix) {
				std::cout << y++ << "\t";
				for (const auto& elem : row) {
					std::cout << elem << "\t";
				}
				std::cout << std::endl;
			}
		}
	}

	void floydWarshall() {
		int V = adj_matrix.size();
		double dist[V][V];

		int i, j, k;

		/* Initialize the solution matrix same as input graph
		   matrix. Or we can say the initial values of shortest
		   distances are based
		   on shortest paths considering no intermediate vertex.
		 */
		for (i = 0; i < V; i++)
			for (j = 0; j < V; j++)
				dist[i][j] = adj_matrix[i][j];

		/* Add all vertices one by one to the set of
		  intermediate vertices.
		  ---> Before start of a iteration, we have shortest
		  distances between all pairs of vertices such that the
		  shortest distances consider only the vertices in set
		  {0, 1, 2, .. k-1} as intermediate vertices.
		  ----> After the end of a iteration, vertex no. k is
		  added to the set of
		  intermediate vertices and the set becomes {0, 1, 2, ..
		  k} */
		for (k = 0; k < V; k++) {
			// Pick all vertices as source one by one
			for (i = 0; i < V; i++) {
				// Pick all vertices as destination for the
				// above picked source
				for (j = 0; j < V; j++) {
					// If vertex k is on the shortest path from
					// i to j, then update the value of
					// dist[i][j]
					if (dist[i][k] + dist[k][j] < dist[i][j])
						dist[i][j] = dist[i][k] + dist[k][j];
				}
			}
		}

		sol_matrix.resize(V);

		for (i = 0; i < V; i++) {
			sol_matrix[i].resize(V, 0);
		}

		for (i = 0; i < V; i++) {
			for (j = 0; j < V; j++) {
				sol_matrix[i][j] = dist[i][j];
				// sol_matrix[i][j] = 1;
				// std::cout << std::setw(2) << dist[i][j] << " ";
			}
			// std::cout << std::endl;
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

	graph.print(ADJ);
	graph.floydWarshall();
	graph.print(SOL);

	return 0;
}
