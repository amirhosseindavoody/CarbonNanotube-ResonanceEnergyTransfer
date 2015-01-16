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
  
  type (cnt) :: firstCNT, secondCNT
  
  
  
  firstCNT = cnt( n_ch1, m_ch1, nkg)
  secondCNT = cnt (n_ch2, m_ch2, nkg)
  
  call firstCNT.printProperties()
  call secondCNT.printProperties()
  
  call firstCNT.calculateBands(i_sub1, E_th, Kcm_max)
  call secondCNT.calculateBands(i_sub2, E_th, Kcm_max)
  
  call saveCNTProperties(firstCNT)
  call saveCNTProperties(secondCNT)
	
	!pause
	!stop
  
  call loadExcitonWavefunction(firstCNT)
  call loadExcitonWavefunction(secondCNT)
  
  call findCrossings(firstCNT,secondCNT)
  call saveCrossingPoints(firstCNT,secondCNT)
	
	!call calculateDOStest(firstCNT)
	!call saveDOS(firstCNT)
	
	call calculateTransferRate(firstCNT,secondCNT)
	
  print *,'Finish!!!!'
  pause
  
end program cntForsterEnergyTransfer

