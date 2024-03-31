compile:
	g++ -g -std=c++11 ./main.cpp -o main
run:
	@echo "Compiling & Running the main file" 
	@echo ""
	g++ -g -std=c++11 ./main.cpp -o main && ./main  && make clean

clean:
	@echo ""
	@echo "Removing extra file"
	rm -f *.o main
