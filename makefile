
perc: perc.c stack.c
	gcc-7 -fopenmp -std=c99 -o perc perc.c stack.c graph.c
