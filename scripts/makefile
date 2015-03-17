# this is simple make file

FC=gfortran
#FCFLAGS = -O3 -ffree-line-length-none
FCFLAGS = -Wall -Wtabs -Wconversion -fbounds-check -g -ffree-line-length-none
FCFLAGS += -I.

main.o: main.f90 comparams.o cnt_class.o data_class.o energyShift.o input_class.o kappaMatrix_module.o parallelForster_module.o physicalConstants.o prepareForster_module.o smallFns.o transitionTable_module.o arbitraryAngleForster_module.o
main: comparams.o cnt_class.o data_class.o energyShift.o input_class.o kappaMatrix_module.o parallelForster_module.o physicalConstants.o prepareForster_module.o smallFns.o transitionTable_module.o arbitraryAngleForster_module.o
comparams.o: comparams.f90
cnt_class.o: cnt_class.f90 physicalConstants.o smallFns.o
data_class.o: data_class.f90 comparams.o prepareForster_module.o
energyShift.o: energyShift.f90
input_class.o: input_class.f90
kappaMatrix_module.o: kappaMatrix_module.f90 parallelForster_module.o arbitraryAngleForster_module.o
parallelForster_module.o: parallelForster_module.f90
physicalConstants.o: physicalConstants.f90
prepareForster_module.o: prepareForster_module.f90 cnt_class.o input_class.o
smallFns.o: smallFns.f90
transitionTable_module.o: transitionTable_module.f90 comparams.o
arbitraryAngleForster_module.o: arbitraryAngleForster_module.f90



%: %.f90
	$(FC) -o $@ $^ $(FCFLAGS) -llapack -lblas

%.o: %.f90
	$(FC) -c $< $(FCFLAGS)

# Utility targets
.PHONY: clean

clean:
	rm -f *.o *.mod *.MOD *~ main