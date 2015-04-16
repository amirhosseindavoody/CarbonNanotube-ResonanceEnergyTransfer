!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate transition table
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
module transitionTable_module
	implicit none
	private
	public :: calculateTransitionTable

	real*8, dimension(:,:,:), allocatable :: transitionRate	

	real*8, public :: c2cMin
	real*8, public :: c2cMax
	integer, public :: nc2c
	real*8 :: c2cDistance
	integer :: ic2c
	real*8 :: dc2c
	
	real*8, public :: thetaMax
	real*8, public :: thetaMin
	integer, public :: nTheta
	real*8 :: theta
	integer :: iTheta
	real*8 :: dTheta
	
	
contains
	!**************************************************************************************************************************
	! calculate transition table
	!**************************************************************************************************************************
	
	subroutine calculateTransitionTable (cnt1,cnt2)
		use cnt_class, only: cnt
		use matrix_element_mod, only: calculate_kSpaceMatrixElement, calculate_geometricMatrixElement
		use parallel_geometry_mod, only: calculateParallelGeometryRate
		use physicalConstants, only: pi
		use transition_points_mod, only: findCrossings, findSameEnergy
		use write_log_mod, only: writeLog
		use unparallel_geometry_mod, only: calculateUnparallelGeometryRate

		type(cnt), intent(in) :: cnt1,cnt2
		character(len=100) :: logInput
		
		
		call writeLog(new_line('A')//"************** Start calculating transitionTable ****************")
		
		! set seperation properties
		if (nc2c .ne. 1) then
			dc2c = (c2cMax-c2cMin)/dble(nc2c-1)
		else
			dc2c = 0.d0
		end if

		write(logInput, '("c2cMin[nm] = ", F4.1)') c2cMin*1.d9
		call writeLog(trim(logInput))
		write(logInput, '("c2cMax[nm] = ", F4.1)') c2cMax*1.d9
		call writeLog(trim(logInput))
		write(logInput, '("nC2C = ", I3.3)') nC2C
		
		! set orientation properties
		if (nTheta .ne. 1) then
			dTheta = (thetaMax-thetaMin)/dble(nTheta-1)
		else
			dTheta = 0.d0
		end if

		call writeLog(trim(logInput))
		write(logInput, '("thetaMin = ", I3.3)') nint(thetaMin*180/pi)
		call writeLog(trim(logInput))
		write(logInput, '("thetaMax = ", I3.3)') nint(thetaMax*180/pi)
		call writeLog(trim(logInput))
		write(logInput, '("nTheta = ", I3.3)') nTheta
		call writeLog(trim(logInput))
	
		!calculate the crossing points and points with the same energy between cnt1 and cnt2
! 		call findSameEnergy(cnt1,cnt2)
! 		call calculate_kSpaceMatrixElement()
		call calculate_geometricMatrixElement(thetaMin, thetaMax, nTheta, c2cDistance)
		call exit()
			
		!allocate the transition rate table
		allocate(transitionRate(2,nTheta,nc2c))
		
		do ic2c = 1, nc2c
			do iTheta = 1, nTheta
		
				c2cDistance = c2cMin+dble(ic2c-1)*dc2c
				theta = thetaMin + dble(iTheta-1)*dTheta
			
				if (theta .eq. 0.d0) then
					!call calculateParallelGeometryRate(cnt1,cnt2, transitionRate(1,iTheta,ic2c), transitionRate(2,iTheta,ic2c), c2cDistance)
					transitionRate(1,iTheta,ic2c) = 0.d0
					transitionRate(2,iTheta,ic2c) = 0.d0
				else
					call calculateUnparallelGeometryRate(cnt1,cnt2, transitionRate(1,iTheta,ic2c), transitionRate(2,iTheta,ic2c), c2cDistance, theta)
				end if
				
				write(logInput,*) 'iTheta=', iTheta, ', nTheta=', nTheta, 'iC2C=', ic2c, ', nC2C=', nc2c
				call writeLog(logInput)		
			end do
		end do
		
		call saveTransitionRates()
		return				
	end subroutine calculateTransitionTable
	
	!**************************************************************************************************************************
	! save the calculated transition table
	!**************************************************************************************************************************
	subroutine saveTransitionRates()

		!write transition rates to the file
		open(unit=100,file='transitionRates12.dat',status="unknown")
		do ic2c = 1,nc2c
			do iTheta=1,nTheta
				write(100,10, advance='no') transitionRate(1,iTheta,ic2c)
			end do
			write(100,10)
		end do
		close(100)
		
		open(unit=100,file='transitionRates21.dat',status="unknown")
		do ic2c = 1,nc2c
			do iTheta=1,nTheta
				write(100,10, advance='no') transitionRate(2,iTheta,ic2c)
			end do
			write(100,10)
		end do
		close(100)
				
		open(unit=100,file='theta.dat',status="unknown")
		do iTheta=1,nTheta
			write(100,10, advance='no') thetaMin+dble(iTheta-1)*dTheta
		enddo
		write(100,10)
				
		open(unit=100,file='c2c.dat',status="unknown")
		do ic2c=1,nc2c
			write(100,10, advance='no') c2cMin+dble(ic2c-1)*dc2c
		enddo
		close(100)
				
10		FORMAT (E16.8)
		
		return
	end subroutine saveTransitionRates
end module transitionTable_module