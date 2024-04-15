<h1 align="center">All Pair Shortest Path</h1>

## Background

#### Problem:

The problem is to find the shortest path between all pairs of nodes in the graph.

#### Input:

A csv file contating the data of a directed and weighted graph is used for the computation of the problem. All the weights are assumed to be positive. The expected format of the file is `src: int, weight: double, dst: int`.

#### Output:

The output is a `N x N` matrix printed to the console which can also be directed to a file. At the end of the output file, the time taken by the algorithm is also mentioned.

## Instructions

-   Create a build folder to store complied files:
    ` mkdir build`
-   Create an out folder to store output files:
    `mkdir out`
-   Use the sample graphs in the input_graphs folder
-   `make` can compile the serial, parallel & distributed implementations

#### Serial

-   Compile serial code:
    `make serial`
-   To run the file "100Edges.csv" from input_graphs (the --nodes parameter is the matrix dimension, in this case, it's 100):
    `./build/main_serial --nodes 100 --file "./input_graphs/100Edges.csv" > ./out/serial_100.txt`
-   To run the file "10Edges.csv" from input_graphs (the --nodes parameter is the matrix dimension, in this case, it's 10):
    `./build/main_serial --nodes 10 --file "./input_graphs/10Edges.csv" > ./out/serial_10.txt`

#### Parallel

-   Compile parallel code:
    `make parallel`
-   To run the file "100Edges.csv" from input_graphs (the --nodes parameter is the matrix dimension, in this case, it's 100):
    `./build/main_parallel --nodes 100 --file "./input_graphs/100Edges.csv" > ./out/parallel_100.txt`
-   To change the number of threads, change the "NUM_THREADS" value in main_parallel.cpp

#### Distributed

-   Compile distributed code:
    `make mpi`
-   To run the file "100Edges.csv" from input_graphs with 4 processes (the --nodes parameter is the matrix dimension, in this case, it's 100):
    ` mpirun np - 4 ./build/main_mpi --nodes 100 --file "./input_graphs/100Edges.csv" > ./out/mpi_100.txt`

## Credits:

-   [CSV Parser](https://github.com/ben-strasser/fast-cpp-csv-parser)
-   [Command Line Parser](https://www.sfu.ca/computing/people/faculty/alaa-alameldeen.html)
-   [Data Set](https://github.com/neryabigon/Weighted-directed-graphs/tree/main)
-   [MPI Floyd-Warshall](https://github.com/nuno-azevedo/floyd-warshall-mpi/tree/master)
-   Contributors:
    -   [Mehar](https://github.com/singhmeharjeet)
    -   [Vikram](https://github.com/vikramsodhan)
    -   [Anh Vu](https://github.com/codingpog)
