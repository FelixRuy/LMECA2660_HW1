CC=gcc
CFLAGS=-O2 -Wall -g

CFL=0.5
RESOL=4
METHOD=I4
FILE_SOL=data/comp_sol.txt
FILE_DIAG=data/diag.tkt

exec : main.c vector.c PDE.c
	rm -f $(FILE_SOL)
	rm -f $(FILE_DIAG)
	rm -f data/analytical_solution.txt
	$(CC) $(CFLAGS) -o $@ $^
	./exec $(CFL) $(RESOL) $(METHOD) $(FILE_SOL) $(FILE_DIAG) 
	rm -f exec

execPlot :
	python3 plot_convect

