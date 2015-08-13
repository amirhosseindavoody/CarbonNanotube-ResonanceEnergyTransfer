!*******************************************************************************
! This program calculates the excitonic band structure of single wall carbon nanotubes through simple tight-binding method.
! Amirhossein Davoody
! Last modified: 3/17/2014
!*******************************************************************************

program cnt_resonance_energy_transfer
	use cnt_class, only: cnt
	use comparams, only: starttime, endtime, cnt1, cnt2
	use input_cnt_mod, only: input_cnt
	use occupation_mod, only: calculate_occupation_table
	use parse_input_file_mod, only: parse_input_file
	use prepareForster_module, only: saveDOS
	use target_exciton_mod, only: set_target_exciton
	use transition_table_mod, only: calculate_transition_table, save_transition_rates
	use write_log_mod, only: writeLog
	
	implicit none

	character(len=100) :: logInput

	call CPU_time(starttime)

	call parse_input_file()
	
	call writeLog(new_line('A')//"************** Reading cnt1 ****************")
	call input_cnt(cnt1)

	call writeLog(new_line('A')//"************** Reading cnt2 ****************")
	call input_cnt(cnt2)

	! transfer rate from Ex0_A2 to other exciton types
	call writeLog(new_line('A')//"Ex0_A2 --> Ex0_A2")
	call set_target_exciton(cnt1, 'Ex0_A2')
	call set_target_exciton(cnt2, 'Ex0_A2')
	call calculate_transition_table(cnt1,cnt2)

	call writeLog(new_line('A')//"Ex0_A2 --> Ex0_Ep")
	call set_target_exciton(cnt1, 'Ex0_A2')
	call set_target_exciton(cnt2, 'Ex0_Ep')
	call calculate_transition_table(cnt1,cnt2)

	call writeLog(new_line('A')//"Ex0_A2 --> Ex0_Em")
	call set_target_exciton(cnt1, 'Ex0_A2')
	call set_target_exciton(cnt2, 'Ex0_Em')
	call calculate_transition_table(cnt1,cnt2)

	! transfer rate from Ex0_Ep to other exciton types
	call writeLog(new_line('A')//"Ex0_Ep --> Ex0_A2")
	call set_target_exciton(cnt1, 'Ex0_Ep')
	call set_target_exciton(cnt2, 'Ex0_A2')
	call calculate_transition_table(cnt1,cnt2)

	call writeLog(new_line('A')//"Ex0_Ep --> Ex0_Ep")
	call set_target_exciton(cnt1, 'Ex0_Ep')
	call set_target_exciton(cnt2, 'Ex0_Ep')
	call calculate_transition_table(cnt1,cnt2)

	call writeLog(new_line('A')//"Ex0_Ep --> Ex0_Em")
	call set_target_exciton(cnt1, 'Ex0_Ep')
	call set_target_exciton(cnt2, 'Ex0_Em')
	call calculate_transition_table(cnt1,cnt2)

	! transfer rate from Ex0_Em to other exciton types
	call writeLog(new_line('A')//"Ex0_Em --> Ex0_A2")
	call set_target_exciton(cnt1, 'Ex0_Em')
	call set_target_exciton(cnt2, 'Ex0_A2')
	call calculate_transition_table(cnt1,cnt2)

	call writeLog(new_line('A')//"Ex0_Em --> Ex0_Ep")
	call set_target_exciton(cnt1, 'Ex0_Em')
	call set_target_exciton(cnt2, 'Ex0_Ep')
	call calculate_transition_table(cnt1,cnt2)

	call writeLog(new_line('A')//"Ex0_Em --> Ex0_Em")
	call set_target_exciton(cnt1, 'Ex0_Em')
	call set_target_exciton(cnt2, 'Ex0_Em')
	call calculate_transition_table(cnt1,cnt2)

	call save_transition_rates()

	call CPU_time(endtime)
	write(logInput,'("Run time = ",f10.3," seconds.")'),endtime-starttime
	call writeLog(logInput)

end program cnt_resonance_energy_transfer

