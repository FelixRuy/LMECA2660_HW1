CC=gcc
CFLAGS=-O2 -Wall -g

MODE=1	
INIITIAL_CONDITION=0
CFL=0.5
RESOL=8
METHOD=E4
FILE_SOL=data/comp_sol.txt
FILE_SOL_PHYS=data/comp_sol_phys.txt
FILE_DIAG=data/diag.txt

exec : main.c vector.c PDE.c
	rm -f $(FILE_SOL)
	rm -f $(FILE_DIAG)
	rm -f $(FILE_SOL_PHYS)
	rm -f data/analytical_solution.txt
	rm -f data/analytical_solution_phys.txt
	$(CC) $(CFLAGS) -o $@ $^
	./exec $(CFL) $(RESOL) $(METHOD) $(MODE) $(INIITIAL_CONDITION) $(FILE_SOL) $(FILE_DIAG) $(FILE_SOL_PHYS) 
	rm -f exec

execPlot :
	python3 plot_convect

