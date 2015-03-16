!*******************************************************************************
! This program calculates the excitonic band structure of single wall carbon nanotubes through simple tight-binding method.
! Amirhossein Davoody
! Last modified: 3/17/2014
!*******************************************************************************

program cntForsterEnergyTransfer
  use cntClass
  use inputParameters
  use dataClass
	use prepareForster_module
	use parallelForster_module
	use arbitraryAngleForster_module
	use kappaMatrix_module
	use transitionTable_module
	use energyCorrection_module
	
  implicit none
	
  type (cnt) :: cnt1, cnt2
	
  cnt1 = cnt( n_ch1, m_ch1, nkg)
  cnt2 = cnt (n_ch2, m_ch2, nkg)
  
  call cnt1.printProperties()
  call cnt2.printProperties()
  
  call cnt1.calculateBands(i_sub1, E_th, Kcm_max)
  call cnt2.calculateBands(i_sub2, E_th, Kcm_max)
	
	cnt1.excitonDirectory = "CNT(07,05)-nkg(1001)-nr(0200)-E_th(1.5)-Kcm_max(1.5)-i_sub(1)-kappa(5.4495)"
	cnt2.excitonDirectory = "CNT(08,07)-nkg(1001)-nr(0200)-E_th(1.5)-Kcm_max(1.5)-i_sub(1)-kappa(6.336)"
  write(outputDirectory,"('ForsterRate (',I2.2,',',I2.2,') to (',I2.2,',',I2.2,')')") n_ch1, m_ch1, n_ch2, m_ch2
	
  !call saveCNTProperties(cnt1)
  !call saveCNTProperties(cnt2)
  
  call loadExcitonWavefunction(cnt1)
  call loadExcitonWavefunction(cnt2)
	
	!call correctEnergy(cnt1)
	!call correctEnergy(cnt2)
		
	call calculateTransitionTable(cnt1,cnt2)
	
	!call calculateKappaMatrix(cnt1,cnt2)
	
	print *, ''
  print *, 'Press Enter to continue ...'
  pause
  
end program cntForsterEnergyTransfer

