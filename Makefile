all: serial parallel mpi

serial:
	g++ -g -std=c++11 ./src/main_serial.cpp -D PRINT -o ./build/main_serial -l pthread

parallel:
	g++ -g -std=c++11 ./src/main_parallel.cpp -D PRINT -o ./build/main_parallel -l pthread

mpi:
	mpic++ -g -std=c++11 ./src/main_mpi.cpp -D PRINT -o ./build/main_mpi -l pthread

datamaker:
	g++ -std=c++11 datamaker.cpp -o data; ./data --nodes 100

serial_10:
	./build/main_serial --nodes 10 --file "./input_graphs/10Edges.csv" > ./out/serial_10.txt

serial_500:
	./build/main_serial --nodes 500 --file "./input_graphs/500Edges.csv" > ./out/serial_500.txt

serial_1000:
	./build/main_serial --nodes 1000 --file "./input_graphs/1000Edges.csv" > ./out/serial_1000.txt

parallel_10:
	./build/main_parallel --nodes 10 --file "./input_graphs/10Edges.csv" > ./out/parallel_10.txt

parallel_500:
	./build/main_parallel --nodes 500 --file "./input_graphs/500Edges.csv" > ./out/parallel_500.txt

parallel_1000:
	./build/main_parallel --nodes 1000 --file "./input_graphs/1000Edges.csv" > ./out/parallel_1000.txt

mpi_10:
	./build/main_mpi --nodes 10 --file "./input_graphs/10Edges.csv" > ./out/mpi_10.txt

mpi_500:
	./build/main_mpi --nodes 500 --file "./input_graphs/500Edges.csv" > ./out/mpi_500.txt

mpi_1000:
	./build/main_mpi --nodes 1000 --file "./input_graphs/1000Edges.csv" > ./out/mpi_1000.txt

clean:
	rm -f -r ./build/*
	rm -f ./out/*.txt

