solution: solution.cpp
	clang++ solution.cpp -o solution -std=c++20

PQ-Dijkstra: PQ-Dijkstra.c
	gcc -O3 PQ-Dijkstra.c -o PQ-Dijkstra
	@echo "**Note** This is only used to build the skelton."
	@echo "To compile your own code, use 'make solution'."

clean:
	rm solution -f
