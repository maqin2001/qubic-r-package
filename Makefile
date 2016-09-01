all:
	g++ -std=c++0x -DNDEBUG -fopenmp src/*.cc src/*.c -O3 -o qubic

clean:
	rm -f qubic
