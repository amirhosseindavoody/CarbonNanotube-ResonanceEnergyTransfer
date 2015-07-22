!*******************************************************************************
! This program calculates the excitonic band structure of single wall carbon nanotubes through simple tight-binding method.
! Amirhossein Davoody
! Last modified: 3/17/2014
!*******************************************************************************

program cnt_resonance_energy_transfer
	use cnt_class, only: cnt
	use comparams, only: starttime, endtime, cnt1, cnt2
	use free_electron_mod, only: calculate_free_electron_transition, save_free_electron_transition
	use input_cnt_mod, only: input_cnt
	use occupation_mod, only: calculate_occupation_table
	use parse_input_file_mod, only: parse_input_file
	use prepareForster_module, only: saveDOS
	use transition_table_mod, only: calculate_transition_table
	use write_log_mod, only: writeLog
	
	implicit none

	character(len=100) :: logInput

	call CPU_time(starttime)

	call parse_input_file()
	
	call writeLog(new_line('A')//"************** Reading cnt1 ****************")
	call input_cnt(cnt1)
	call calculate_free_electron_transition(cnt1)
! 	call save_free_electron_transition(cnt1)
! 	call exit()

	call writeLog(new_line('A')//"************** Reading cnt2 ****************")
	call input_cnt(cnt2)
	call calculate_free_electron_transition(cnt2)

 	call calculate_transition_table(cnt1,cnt2)

	call CPU_time(endtime)
	write(logInput,'("Run time = ",f10.3," seconds.")'),endtime-starttime
	call writeLog(logInput)

end program cnt_resonance_energy_transfer

