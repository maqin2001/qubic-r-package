all:
	g++ -std=c++0x src/*.cc src/*.c -O3 -o qubic

clean:
	rm -f qubic
