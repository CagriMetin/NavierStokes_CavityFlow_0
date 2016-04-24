
FC = gfortran
FFLAGS = -O3
LFLAGS = $(FFLAGS)

.PHONY: plots, clean, clobber

%.o : %.f90
	$(FC) $(FFLAGS) -c $<

run.exe: Cavity.o
	$(FC) $(LFLAGS) Cavity.o -o run.exe

solution.txt: run.exe
	@echo 
	@echo Running code...
	./run.exe

plots: solution.txt
	@echo 
	@echo Plotting results...
	python plot_solution.py

clean:
	rm -f *.o *.exe

clobber: clean
	rm -f *.txt *.png

