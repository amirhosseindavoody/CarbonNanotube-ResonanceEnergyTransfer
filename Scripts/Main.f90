!*******************************************************************************
! This program calculates the excitonic band structure of single wall carbon nanotubes through simple tight-binding method.
! Amirhossein Davoody
! Last modified: 3/17/2014
!*******************************************************************************

program cnt_resonance_energy_transfer
	use comparams
	use cnt_class
	use input_class
	use data_class
	use prepareForster_module
	use parallelForster_module
	use arbitraryAngleForster_module
	use kappaMatrix_module
	use transitionTable_module
	use energyShift_module
	use smallFunctions
	use output_module
	implicit none

	character(len=100) :: dirname
	type (cnt) :: cnt1, cnt2

	call CPU_time(starttime)

	! read the information of CNTs stored the the following folders
	write(dirname,'(A)') "CNT(07,05)-nkg(1001)-nr(0200)-E_th(1.5)-Kcm_max(1.5)-i_sub(1)-Ckappa(0.5)"
	call inputCNT(cnt1,dirname)
	write(dirname,'(A)') "CNT(07,05)-nkg(1001)-nr(0200)-E_th(1.5)-Kcm_max(1.5)-i_sub(1)-Ckappa(0.5)"
	call inputCNT(cnt2,dirname)


	! specify output directory
	write(outdir,"('Transfer-(',I2.2,',',I2.2,')to(',I2.2,',',I2.2,')-Ckappa(0.5)')") cnt1%n_ch, cnt1%m_ch, cnt2%n_ch, cnt2%m_ch
	call open_output_directory()
	print *, "I am here 3"

! 	call exit()
	
	call calculateTransitionTable(cnt1,cnt2)
	
! 	call calculateKappaMatrix(cnt1,cnt2)

	call close_output_directory()


end program cnt_resonance_energy_transfer

