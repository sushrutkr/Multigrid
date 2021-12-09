FC=gfortran
FC1=ifort
FFLAGS=-O3
SRC=modules.F90 		\
	multigrid.F90 		\
	project.F90			\


OBJ=${SRC:.F90=.o}
FNAME=multigrid

%.o : %.F90
	$(FC) $(FFLAGS) -o $@ -c $<

mac : $(OBJ)
	$(FC) $(FFLAGS) -o $(FNAME) $(OBJ)

lab : $(OBJ)
	$(FC1) $(FFLAGS) -o $(FNAME) $(OBJ)

clean:
	@rm *.mod *.o 