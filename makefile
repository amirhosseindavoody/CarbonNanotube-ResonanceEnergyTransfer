
FC=gfortran
# FCFLAGS = -O3 -ffree-line-length-none
FCFLAGS += -fcheck=all -Wall -Wconversion-extra -Wtabs -Wimplicit-interface -Wintrinsic-shadow -Wsurprising -Wfunction-elimination -fbounds-check -ffree-line-length-none -fall-intrinsics

SRCDIR = ./src
OBJDIR = ./obj

a2a_kspace_matrix_element_mod.o: comparams.o physicalConstants.o transition_points_mod.o write_log_mod.o
a2ep_kspace_matrix_element_mod.o: comparams.o physicalConstants.o transition_points_mod.o write_log_mod.o
a2em_kspace_matrix_element_mod.o: comparams.o physicalConstants.o transition_points_mod.o write_log_mod.o
cnt_class.o: physicalConstants.o math_functions_mod.o
comparams.o: cnt_class.o
density_of_states_mod.o: cnt_class.o
em2a_kspace_matrix_element_mod.o: comparams.o physicalConstants.o transition_points_mod.o write_log_mod.o
em2ep_kspace_matrix_element_mod.o: comparams.o physicalConstants.o transition_points_mod.o write_log_mod.o
em2em_kspace_matrix_element_mod.o: comparams.o physicalConstants.o transition_points_mod.o write_log_mod.o	
ep2a_kspace_matrix_element_mod.o: comparams.o physicalConstants.o transition_points_mod.o write_log_mod.o
ep2ep_kspace_matrix_element_mod.o: comparams.o physicalConstants.o transition_points_mod.o write_log_mod.o
ep2em_kspace_matrix_element_mod.o: comparams.o physicalConstants.o transition_points_mod.o write_log_mod.o
geometric_matrix_element_mod.o: comparams.o physicalConstants.o
input_cnt_mod.o: cnt_class.o physicalConstants.o write_log_mod.o
main.o: cnt_class.o comparams.o density_of_states_mod.o input_cnt_mod.o occupation_mod.o parse_input_file_mod.o partition_function_mod.o target_exciton_mod.o transition_table_mod.o write_log_mod.o
occupation_mod.o: comparams.o cnt_class.o physicalConstants.o
parse_input_file_mod.o: comparams.o physicalConstants.o transition_table_mod.o write_log_mod.o
partition_function_mod.o: cnt_class.o comparams.o physicalConstants.o
rotate_shift_mod.o: cnt_class.o math_functions_mod.o physicalConstants.o
target_exciton_mod.o: cnt_class.o write_log_mod.o
transition_points_mod.o: cnt_class.o comparams.o math_functions_mod.o physicalConstants.o write_log_mod.o
transition_table_mod.o: a2a_kspace_matrix_element_mod.o a2ep_kspace_matrix_element_mod.o a2em_kspace_matrix_element_mod.o cnt_class.o comparams.o density_of_states_mod.o ep2a_kspace_matrix_element_mod.o ep2em_kspace_matrix_element_mod.o ep2ep_kspace_matrix_element_mod.o em2a_kspace_matrix_element_mod.o em2em_kspace_matrix_element_mod.o em2ep_kspace_matrix_element_mod.o geometric_matrix_element_mod.o physicalConstants.o partition_function_mod.o rotate_shift_mod.o transition_points_mod.o write_log_mod.o

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