# this is simple make file

FC=gfortran
#FCFLAGS = -O3 -ffree-line-length-none
FCFLAGS = -Wall -Wtabs -Wconversion -fbounds-check -g -ffree-line-length-none
FCFLAGS += -I.

main.o: main.f90 CNT_class.o dataclass.o energyCorrection.o inputParameters.o kappaMatrix_module.o parallelForster_module.o physicalConstants.o prepareForster_module.o smallFns.o transitionTable_module.o arbitraryAngleForster_module.o
main: CNT_class.o dataclass.o energyCorrection.o inputParameters.o kappaMatrix_module.o parallelForster_module.o physicalConstants.o prepareForster_module.o smallFns.o transitionTable_module.o arbitraryAngleForster_module.o
CNT_class.o: CNT_class.f90 physicalConstants.o smallFns.o
dataclass.o: dataclass.f90 inputParameters.o prepareForster_module.o
energyCorrection.o: energyCorrection.f90
inputParameters.o: inputParameters.f90
kappaMatrix_module.o: kappaMatrix_module.f90 parallelForster_module.o arbitraryAngleForster_module.o
parallelForster_module.o: parallelForster_module.f90
physicalConstants.o: physicalConstants.f90
prepareForster_module.o: prepareForster_module.f90 CNT_class.o
smallFns.o: smallFns.f90
transitionTable_module.o: transitionTable_module.f90
arbitraryAngleForster_module.o: arbitraryAngleForster_module.f90



%: %.f90
	$(FC) -o $@ $^ $(FCFLAGS) -llapack -lblas

%.o: %.f90
	$(FC) -c $< $(FCFLAGS)

# Utility targets
.PHONY: clean

clean:
	rm -f *.o *.mod *.MOD *~ main