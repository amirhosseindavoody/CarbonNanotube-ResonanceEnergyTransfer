
FC=gfortran
FCFLAGS = -O3 -ffree-line-length-none
# FCFLAGS += -Wall -Wconversion-extra -Wtabs -Wimplicit-interface -Wintrinsic-shadow -Wsurprising -Wfunction-elimination
# FCFLAGS += -fbounds-check -ffree-line-length-none -fall-intrinsics

SRCDIR = ./src
OBJDIR = ./obj

main.o: cnt_class.o kappaMatrix_module.o transitionTable_module.o parse_input_file_mod.o input_cnt_mod.o
cnt_class.o: physicalConstants.o math_functions_mod.o
kappaMatrix_module.o: physicalConstants.o parallelForster_module.o arbitraryAngleForster_module.o output_module.o cnt_class.o prepareForster_module.o
parallelForster_module.o: comparams.o physicalConstants.o math_functions_mod.o output_module.o cnt_class.o prepareForster_module.o 
prepareForster_module.o: cnt_class.o physicalConstants.o comparams.o output_module.o
transitionTable_module.o: parallelForster_module.o arbitraryAngleForster_module.o prepareForster_module.o physicalConstants.o cnt_class.o
output_module.o: comparams.o cnt_class.o
input_cnt_mod.o: cnt_class.o physicalConstants.o
arbitraryAngleForster_module.o: physicalConstants.o prepareForster_module.o cnt_class.o
energyShiftMOD.o: math_functions_mod.o physicalConstants.o cnt_class.o
parse_input_file_mod.o: comparams.o physicalConstants.o output_module.o
comparams.o: cnt_class.o

main: main.o
	$(FC) $(FCFLAGS) -o $@.exe $(wildcard $(OBJDIR)/*.o) 

%.o: $(SRCDIR)/%.f90 | $(OBJDIR)
	$(FC) -c $(FCFLAGS) $< -J$(OBJDIR)
	@mv -f $@ $(OBJDIR)

$(OBJDIR):
	@mkdir -p $(OBJDIR)

# Utility targets
.PHONY: clean
clean:
	@rm -f *.o *.mod *.exe
	@rm -rf $(OBJDIR)