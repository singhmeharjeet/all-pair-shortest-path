compile:
	g++ -g -std=c++11 ./main_serial.cpp -o main_serial
	g++ -g -std=c++11 ./main_parallel.cpp -o main_parallel
	g++ -g -std=c++11 ./main_mpi.cpp -o main_mpi

serial:
	@echo "Compiling & Running the serial file" 
	@echo ""
	g++ -g -std=c++11 ./main_serial.cpp -o main_serial && ./main_serial  && make clean

serial_print:
	@echo "Compiling & Running the serial file" 
	@echo ""
	g++ -g -std=c++11 ./main_serial.cpp -D PRINT -o main_serial && ./main_serial  && make clean

serial2:
	@echo "Compiling & Running the serial file" 
	@echo ""
	g++ -g -std=c++11 ./main_serial.cpp -o main_serial && ./main_serial --numNodes 12 --edgesFile "./input_graphs/12Edges.csv"  && make clean

parallel:
	@echo "Compiling & Running the parallel file" 
	@echo ""
	g++ -g -std=c++11 ./main_parallel.cpp -o main_parallel && ./main_parallel  && make clean

mpi:
	@echo "Compiling & Running the mpi file" 
	@echo ""
	g++ -g -std=c++11 ./main_mpi.cpp -o main_mpi && ./main_mpi  && make clean

clean:
	@echo ""
	@echo "Removing extra file"
	rm -f *.o main_serial main_parallel main_mpi data
