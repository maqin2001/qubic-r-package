all:
	g++ -std=c++0x -fopenmp src/*.cc src/*.c -O3 -o qubic

clean:
	rm -f qubic
