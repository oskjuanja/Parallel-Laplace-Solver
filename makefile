all: Jac MG

Jac: PDEsolvers.c PDEsolvers.h Jacobi.c
	gcc -lpthread PDEsolvers.c Jacobi.c -o $@

MG: PDEsolvers.c PDEsolvers.h multigrid.c
	gcc -lpthread PDEsolvers.c multigrid.c -o $@