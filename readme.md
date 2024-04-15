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

## Instructions

-  Create a build folder to store complied files:
``` mkdir build```

-  Create an out folder to store output files:
``` mkdir out ```

-  Use the sample graphs in the input_graphs folder

#### Serial
-  Compile serial code:
```make serial```

- To run the file "100Edges.csv" from input_graphs (the --nodes parameter is the matrix dimension, in this case, it's 100):
```./build/main_serial --nodes 100 --file "./input_graphs/100Edges.csv" > ./out/serial_100.txt```

- To run the file "10Edges.csv" from input_graphs (the --nodes parameter is the matrix dimension, in this case, it's 10):
```./build/main_serial --nodes 10 --file "./input_graphs/10Edges.csv" > ./out/serial_10.txt```

#### Parallel
-  Compile parallel code:
```make parallel```

- To run the file "100Edges.csv" from input_graphs (the --nodes parameter is the matrix dimension, in this case, it's 100):
```./build/main_parallel --nodes 100 --file "./input_graphs/100Edges.csv" > ./out/parallel_100.txt```

- To change the number of threads, change the "NUM_THREADS" value in main_parallel.cpp

#### Distributed
-  Compile distributed code:
```make mpi```

- To run the file "100Edges.csv" from input_graphs with 4 processes (the --nodes parameter is the matrix dimension, in this case, it's 100):
``` mpirun np - 4 ./build/main_mpi --nodes 100 --file "./input_graphs/100Edges.csv" > ./out/mpi_100.txt```


## Credits:

-   [CSV Parser](https://github.com/ben-strasser/fast-cpp-csv-parser)
-   [Command Line Parser](https://www.sfu.ca/computing/people/faculty/alaa-alameldeen.html)
-   [Data Set](https://github.com/neryabigon/Weighted-directed-graphs/tree/main)
-   [MPI Floyd-Warshall](https://github.com/nuno-azevedo/floyd-warshall-mpi/tree/master)
-   Contributors:
    -   [Mehar](https://github.com/singhmeharjeet)
    -   [Vikram](https://github.com/vikramsodhan)
    -   [Anh Vu](https://github.com/codingpog)
