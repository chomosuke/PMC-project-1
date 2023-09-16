solution: solution.cpp
	clang++ -O3 solution.cpp -o solution -std=c++20

PQ-Dijkstra: PQ-Dijkstra.c
	clang -O3 PQ-Dijkstra.c -o PQ-Dijkstra

clean:
	rm solution PQ-Dijkstra -f
