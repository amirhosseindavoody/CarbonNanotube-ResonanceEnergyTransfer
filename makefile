
FC=gfortran
FCFLAGS = -g -ffree-line-length-none
FCFLAGS += -Wall -Wconversion-extra -Wtabs -Wimplicit-interface -Wintrinsic-shadow -Wsurprising -Wfunction-elimination
FCFLAGS += -fbounds-check -ffree-line-length-none -fall-intrinsics

SRCDIR = ./src
OBJDIR = ./obj

main.o: comparams.o cnt_class.o input_class.o kappaMatrix_module.o transitionTable_module.o output_module.o
cnt_class.o: physicalConstants.o mathFunctionsMOD.o
kappaMatrix_module.o: physicalConstants.o parallelForster_module.o arbitraryAngleForster_module.o output_module.o cnt_class.o input_class.o prepareForster_module.o
parallelForster_module.o: comparams.o physicalConstants.o mathFunctionsMOD.o output_module.o cnt_class.o prepareForster_module.o 
prepareForster_module.o: cnt_class.o physicalConstants.o comparams.o
transitionTable_module.o: input_class.o parallelForster_module.o arbitraryAngleForster_module.o prepareForster_module.o physicalConstants.o cnt_class.o
output_module.o: mathFunctionsMOD.o comparams.o cnt_class.o
input_class.o: cnt_class.o physicalConstants.o
arbitraryAngleForster_module.o: input_class.o physicalConstants.o prepareForster_module.o cnt_class.o
energyShiftMOD.o: mathFunctionsMOD.o physicalConstants.o cnt_class.o
output_module.o: mathFunctionsMOD.o

main: main.o
	$(FC) -o $@.exe $(wildcard $(OBJDIR)/*.o) $(FCFLAGS) 

%.o: $(SRCDIR)/%.f90 | $(OBJDIR)
	$(FC) -c $< $(FCFLAGS) -J$(OBJDIR)
	@mv -f $@ $(OBJDIR)

$(OBJDIR):
	@mkdir -p $(OBJDIR)

# Utility targets
.PHONY: clean
clean:
	@rm -f *.o *.mod *.exe
	@rm -rf $(OBJDIR)