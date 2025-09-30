FC = gfortran
#OPTF = -g -O2 -fbounds-check

OPTF = -g -O3 -pg -fopenmp 

all: main

const.o: const.f90
	$(FC) $(OPTF) -c const.f90
	
var.o: const.o  var.f90
	$(FC) $(OPTF) -c const.o var.f90

init.o : var.o const.o init.f90
	$(FC) $(OPTF) -c var.o const.o init.f90

grad.o: var.o const.o grad.f90
	$(FC) $(OPTF) -c var.o const.o grad.f90

solv.o: var.o const.o solv.f90
	$(FC) $(OPTF) -c var.o const.o solv.f90

pos.o: var.o const.o solv.o CIC.o interp.o grad.o pos.f90
	$(FC) $(OPTF) -c var.o CIC.o interp.o const.o solv.o grad.o pos.f90

CIC.o: const.o var.o CIC.f90
	$(FC) $(OPTF) -c const.o var.o CIC.f90

interp.o : const.o var.o interp.f90
	$(FC) $(OPTF) -c const.o interp.f90

main: var.o const.o init.o grad.o solv.o interp.o pos.o CIC.o main.f90
	$(FC) $(OPTF)  var.o const.o init.o grad.o interp.o solv.o pos.o main.f90 CIC.o -o main


profile_main:
	gprof main_test gmon.out > rapport_profiling_main.txt

clean:
	\rm main main_test *.o *.mod *.dSYM


