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
  use perpendicularForster_module
	use parallelForster_module
	use arbitraryAngleForster_module
	
  implicit none
  
  type (cnt) :: cnt1, cnt2
	integer :: iTheta, nTheta
	real*8 :: dTheta
	
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
  call saveTransitionPoints(cnt1,cnt2)

	call saveDOS(cnt1,cnt2)
	
	!pause
	!stop
	
	call calculateParallelForsterRate(cnt1,cnt2)
	
	print *, ''
	print *, 'Press Enter to continue ...'
	pause
	
	!call calculatePerpendicularForsterRate(cnt1,cnt2)
	!
	!print *, ''
 ! print *, 'Press Enter to continue ...'
 ! pause
	
	nTheta = 10
	dTheta = pi/2.d0/dble(nTheta)
	
	do iTheta = 1, nTheta
	
		theta = dble(iTheta)*dTheta
		
		call calculateArbitraryForsterRate(cnt1,cnt2)
	end do
	
	print *, ''
  print *, 'Press Enter to continue ...'
  pause
  
end program cntForsterEnergyTransfer

