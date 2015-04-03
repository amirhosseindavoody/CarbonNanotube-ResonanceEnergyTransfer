!*******************************************************************************
! This program calculates the excitonic band structure of single wall carbon nanotubes through simple tight-binding method.
! Amirhossein Davoody
! Last modified: 3/17/2014
!*******************************************************************************

program cnt_resonance_energy_transfer
	use cnt_class, only: cnt
	use kappaMatrix_module
	use transitionTable_module
	use output_module, only: open_output_directory, close_output_directory
	use initial_final_mod, only: initialize, finalize
	implicit none

	type (cnt) :: cnt1, cnt2

	

	call initialize(cnt1,cnt2)
	
	call calculateTransitionTable(cnt1,cnt2)
	
! 	call calculateKappaMatrix(cnt1,cnt2)

	call finalize()

end program cnt_resonance_energy_transfer

