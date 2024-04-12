all: serial parallel mpi

serial:
	g++ -g -std=c++11 ./main_serial.cpp -D PRINT -o ./build/main_serial -l pthread

parallel:
	g++ -g -std=c++11 ./main_parallel.cpp -D PRINT -o ./build/main_parallel -l pthread

mpi:
	g++ -g -std=c++11 ./main_mpi.cpp -D PRINT -o ./build/main_mpi -l pthread

datamaker:
	g++ -std=c++11 datamaker.cpp -o data; ./data --nodes 100

serial_500:
	./build/main_serial --nodes 500 --file "./input_graphs/500Edges.csv" > ./out/serial_500.txt

serial_1000:
	./build/main_serial --nodes 1000 --file "./input_graphs/1000Edges.csv" > ./out/serial_1000.txt

parallel_500:
	./build/main_parallel --nodes 500 --file "./input_graphs/500Edges.csv" > ./out/parallel_500.txt

parallel_1000:
	./build/main_parallel --nodes 1000 --file "./input_graphs/1000Edges.csv" > ./out/parallel_1000.txt

clean:
	rm -f -r ./build/*
	rm -f ./out/*.txt

