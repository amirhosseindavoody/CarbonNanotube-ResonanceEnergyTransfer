!*******************************************************************************
! This program calculates the excitonic band structure of single wall carbon nanotubes through simple tight-binding method.
! Amirhossein Davoody
! Last modified: 3/17/2014
!*******************************************************************************

program cnt_resonance_energy_transfer
	use comparams, only: starttime, endtime
	use cnt_class, only: cnt
	use input_class, only: inputCNT
	use kappaMatrix_module
	use transitionTable_module
	use output_module, only: open_output_directory, close_output_directory
	implicit none

	character(len=100) :: buffer
	type (cnt) :: cnt1, cnt2

	call CPU_time(starttime)

	if (command_argument_count() .ne. 2) then
		write(*,*) "ERROR! Input format: main.exe cnt1_folder cnt2_folder"
		call exit()
	end if 


	! read the information of the first CNT
	call get_command_argument(1,buffer)
	buffer = trim(buffer)
	call inputCNT(cnt1,buffer)
	! read the information of the second CNT
	call get_command_argument(2,buffer)
	buffer = trim(buffer)
	call inputCNT(cnt2,buffer)

	!if Ckappa is not equal for both cnt1 and cnt2 print a warning and exit
	if (cnt1%Ckappa .ne. cnt2%Ckappa) then
		write(*,*) "cnt1 and cnt2 have different Ckappa"
		call exit()
	end if

	call open_output_directory(cnt1, cnt2)
	
	call calculateTransitionTable(cnt1,cnt2)
	
! 	call calculateKappaMatrix(cnt1,cnt2)
	
	call CPU_time(endtime)
	call close_output_directory()


end program cnt_resonance_energy_transfer

