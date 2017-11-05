
perc: perc.c stack.c
	mpicc -fopenmp -std=c99 -o perc perc.c stack.c graph.c
