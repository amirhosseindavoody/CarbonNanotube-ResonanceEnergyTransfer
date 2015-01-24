!*******************************************************************************
! This program calculates the excitonic band structure of single wall carbon nanotubes through simple tight-binding method.
! Amirhossein Davoody
! Last modified: 3/17/2014
!*******************************************************************************

program cntForsterEnergyTransfer
  use cntClass
  use inputParameters
  use dataClass
  use forsterClass
  implicit none
  
  type (cnt) :: cnt1, cnt2
	
  cnt1 = cnt( n_ch1, m_ch1, nkg)
  cnt2 = cnt (n_ch2, m_ch2, nkg)
  
  call cnt1.printProperties()
  call cnt2.printProperties()
  
  call cnt1.calculateBands(i_sub1, E_th, Kcm_max)
  call cnt2.calculateBands(i_sub2, E_th, Kcm_max)
  
  !call saveCNTProperties(cnt1)
  !call saveCNTProperties(cnt2)
	
	!pause
	!stop
  
  call loadExcitonWavefunction(cnt1)
  call loadExcitonWavefunction(cnt2)
	
	!pause
	!stop
  
  call findCrossings(cnt1,cnt2)
	call findSameEnergy(cnt1,cnt2)
  call saveCrossingPoints(cnt1,cnt2)
	
	
	
	call saveDOS(cnt1,cnt2)
	
	!pause
	!stop
	
	call calculateTransferRateParallel(cnt1,cnt2)
	
	call calculateTransferRatePerpendicular(cnt1,cnt2)
	
  print *,'Finish!!!!'
  pause
  
end program cntForsterEnergyTransfer

