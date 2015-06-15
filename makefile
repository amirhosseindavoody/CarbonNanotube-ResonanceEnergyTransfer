
FC=gfortran
FCFLAGS = -O3 -ffree-line-length-none
# FCFLAGS += -fcheck=all -Wall -Wconversion-extra -Wtabs -Wimplicit-interface -Wintrinsic-shadow -Wsurprising -Wfunction-elimination
# FCFLAGS += -fbounds-check -ffree-line-length-none -fall-intrinsics

SRCDIR = ./src
OBJDIR = ./obj

a2a_matrix_element_mod.o: comparams.o physicalConstants.o transition_points_mod.o write_log_mod.o
a2ep_matrix_element_mod.o: comparams.o physicalConstants.o transition_points_mod.o write_log_mod.o
cnt_class.o: physicalConstants.o math_functions_mod.o
comparams.o: cnt_class.o
input_cnt_mod.o: cnt_class.o physicalConstants.o write_log_mod.o
kappaMatrix_module.o: cnt_class.o comparams.o parallel_geometry_mod.o physicalConstants.o transition_points_mod.o write_log_mod.o
main.o: cnt_class.o comparams.o input_cnt_mod.o kappaMatrix_module.o occupation_mod.o parse_input_file_mod.o prepareForster_module.o transitionTable_module.o write_log_mod.o
occupation_mod.o: comparams.o cnt_class.o physicalConstants.o
parallel_geometry_mod.o: cnt_class.o comparams.o math_functions_mod.o physicalConstants.o prepareForster_module.o transition_points_mod.o write_log_mod.o  
parse_input_file_mod.o: comparams.o physicalConstants.o transitionTable_module.o write_log_mod.o
prepareForster_module.o: cnt_class.o comparams.o physicalConstants.o
rotate_shift_mod.o: cnt_class.o physicalConstants.o
transition_points_mod.o: cnt_class.o comparams.o math_functions_mod.o physicalConstants.o write_log_mod.o
transitionTable_module.o: a2a_matrix_element_mod.o a2ep_matrix_element_mod.o cnt_class.o comparams.o physicalConstants.o prepareForster_module.o rotate_shift_mod.o transition_points_mod.o write_log_mod.o

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