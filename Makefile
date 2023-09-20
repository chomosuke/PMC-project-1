solution: solution.cc
	clang++ -O3 -fopenmp solution.cc -o solution -std=c++20

debug: solution.cc
	clang++ -O0 -fopenmp solution.cc -o debug -std=c++20 -g

PQ-Dijkstra: PQ-Dijkstra.c
	clang -O3 PQ-Dijkstra.c -o PQ-Dijkstra

clean:
	rm solution PQ-Dijkstra -f
