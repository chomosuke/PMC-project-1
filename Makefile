solution: solution.cpp
	clang++ -O3 -fopenmp solution.cpp -o solution -std=c++20

debug: solution.cpp
	clang++ -O0 -fopenmp solution.cpp -o debug -std=c++20 -g

PQ-Dijkstra: PQ-Dijkstra.c
	clang -O3 PQ-Dijkstra.c -o PQ-Dijkstra

clean:
	rm solution PQ-Dijkstra -f
