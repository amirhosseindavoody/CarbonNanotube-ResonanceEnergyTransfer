!*******************************************************************************
! This program calculates the excitonic band structure of single wall carbon nanotubes through simple tight-binding method.
! Amirhossein Davoody
! Last modified: 3/17/2014
!*******************************************************************************

program cnt_resonance_energy_transfer
	use cnt_class, only: cnt
	use kappaMatrix_module
	use transitionTable_module
	use output_module, only: writeLog
	use parse_input_file_mod, only: parse_input_file
	use comparams, only: starttime, endtime, cnt1, cnt2
	use input_cnt_mod, only: input_cnt
	implicit none

	character(len=100) :: logInput

	call CPU_time(starttime)

	call parse_input_file()
	
	call writeLog(new_line('A')//"************** Reading cnt1 ****************")
	call input_cnt(cnt1)

	call writeLog(new_line('A')//"************** Reading cnt2 ****************")
	call input_cnt(cnt2)

 	call calculateTransitionTable(cnt1,cnt2)
	
! 	call calculateKappaMatrix(cnt1,cnt2)

	call CPU_time(endtime)
	write(logInput,'("Run time = ",f10.3," seconds.")'),endtime-starttime
	call writeLog(logInput)

end program cnt_resonance_energy_transfer

