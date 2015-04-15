
FC=gfortran
FCFLAGS = -g -ffree-line-length-none
FCFLAGS += -Wall -Wconversion-extra -Wtabs -Wimplicit-interface -Wintrinsic-shadow -Wsurprising -Wfunction-elimination
FCFLAGS += -fbounds-check -ffree-line-length-none -fall-intrinsics

SRCDIR = ./src
OBJDIR = ./obj

arbitraryAngleForster_module.o: cnt_class.o comparams.o matrix_element_mod.o physicalConstants.o prepareForster_module.o transition_points_mod.o
cnt_class.o: physicalConstants.o math_functions_mod.o
comparams.o: cnt_class.o
energyShiftMOD.o: math_functions_mod.o physicalConstants.o cnt_class.o
input_cnt_mod.o: cnt_class.o physicalConstants.o write_log_mod.o
kappaMatrix_module.o: arbitraryAngleForster_module.o cnt_class.o comparams.o parallelForster_module.o physicalConstants.o transition_points_mod.o write_log_mod.o
main.o: cnt_class.o input_cnt_mod.o kappaMatrix_module.o write_log_mod.o parse_input_file_mod.o transitionTable_module.o
matrix_element_mod.o: comparams.o transition_points_mod.o write_log_mod.o
parallelForster_module.o: cnt_class.o comparams.o math_functions_mod.o physicalConstants.o prepareForster_module.o transition_points_mod.o write_log_mod.o  
parse_input_file_mod.o: comparams.o physicalConstants.o write_log_mod.o
prepareForster_module.o: cnt_class.o physicalConstants.o comparams.o write_log_mod.o
transitionTable_module.o: parallelForster_module.o arbitraryAngleForster_module.o prepareForster_module.o physicalConstants.o cnt_class.o write_log_mod.o
transition_points_mod.o: cnt_class.o write_log_mod.o

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