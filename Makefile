FC=gfortran -O3 -c
LD=gfortran
SRC=fact.f90 chgm.f90 main.f90
OBJ=fact.o chgm.o main.o
psi:
	$(FC) $(SRC)
	$(LD) $(OBJ) -o psi.x
	rm -rf *.o
clean:
	rm -rf *.o *.x *.pdf

