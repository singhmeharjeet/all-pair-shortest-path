<h1 align="center">All Pair Shortest Path</h1>

## Background

#### Problem:

The problem is to find the shortest path between all pairs of nodes in the graph.

#### Input:

The csv file represents an undirected and weighted graph for the computation of the problem. We are using 2 csv files to represent the graph. The first file contains the nodes and the second file contains the edges

-   The nodes file contains the node id and the node position (x, y, z).
-   The edges file contains the source node id, the destination node id and the weight of the edge.

#### Output:

It will be a 2D matrix of size n x n where n is the number of nodes in the graph. The element at the i-th row and j-th column will be the shortest path between the i-th node and the j-thnode.

Each path will have the following attributes:

-   **Steps**
    -   Format: [start id, .., .., end id]
    -   About: Next node in the array is the next step in the path
    -   Size: n
-   **Costs**
    -   Format: [cost(start, ..), .., cost(.., end)]
    -   About: Cost of each edge
    -   Size: n-1

## Implementations

-   #### Serial

-   #### Parallel

-   #### Distributed

## Analysis:

## Credits:

-   [CSV Parser](https://github.com/ben-strasser/fast-cpp-csv-parser)
-   [Command Line Parser](https://www.sfu.ca/computing/people/faculty/alaa-alameldeen.html)
-   [Data Set](https://github.com/neryabigon/Weighted-directed-graphs/tree/main)
-   Contributors:
    -   [Mehar](https://github.com/singhmeharjeet)
    -   [Vikaram](https://github.com/vikramsodhan)
    -   [Anh Vu](https://github.com/codingpog)
