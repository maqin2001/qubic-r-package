all:
	g++ -std=c++0x -shared src/*.cc src/*.c -O3 -o qubic

clean:
	rm -f qubic
