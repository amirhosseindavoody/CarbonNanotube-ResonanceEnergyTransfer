!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate transition table
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
module transitionTable_module
	use cnt_class
	implicit none
	real*8, dimension(:,:,:), allocatable :: transitionRate
	
	private
	
	public :: calculateTransitionTable
    
contains
	!**************************************************************************************************************************
	! calculate transition table
    !**************************************************************************************************************************
	subroutine calculateTransitionTable (cnt1,cnt2)
		use input_class
		use parallelForster_module
		use arbitraryAngleForster_module
		use prepareForster_module
		use smallFunctions

		type(cnt), intent(in) :: cnt1,cnt2
		
		write(logInput,*) "Start calculating transitionTable"
		call writeLog()
		
		! set seperation properties
		c2cMin = 01.2d-9
		c2cMax = 01.2d-9
		nc2c = 1
		if (nc2c .ne. 1) then
			dc2c = (c2cMax-c2cMin)/dble(nc2c-1)
		else
			dc2c = 0.d0
		end if
		
		! set orientation properties
		thetaMax = pi/2.d0
		thetaMin = 0.d0
		nTheta = 100
		if (nTheta .ne. 1) then
			dTheta = (thetaMax-thetaMin)/dble(nTheta-1)
		else
			dTheta = 0.d0
		end if
		
		!calculate the crossing points and points with the same energy between cnt1 and cnt2
		call findCrossings(cnt1,cnt2)
		call findSameEnergy(cnt1,cnt2)
			
		!allocate the transition rate table
		allocate(transitionRate(2,nTheta,nc2c))
		
		do ic2c = 1, nc2c
			do iTheta = 1, nTheta
		
				c2cDistance = c2cMin+dble(ic2c-1)*dc2c
				theta = thetaMin + dble(iTheta-1)*dTheta
			
				if (theta .eq. 0.d0) then
					!call calculateParallelForsterRate(cnt1,cnt2, transitionRate(1,iTheta,ic2c), transitionRate(2,iTheta,ic2c))
					transitionRate(1,iTheta,ic2c) = 0.d0
					transitionRate(2,iTheta,ic2c) = 0.d0
				else
					call calculateArbitraryForsterRate(cnt1,cnt2, transitionRate(1,iTheta,ic2c), transitionRate(2,iTheta,ic2c))
				end if
				
				write(logInput,*) 'iTheta=', iTheta, ', nTheta=', nTheta, 'iC2C=', ic2c, ', nC2C=', nc2c
				call writeLog()		
			end do
		end do
			
		call saveTransitionRates()
		return				
	end subroutine calculateTransitionTable
	
	!**************************************************************************************************************************
	! save the calculated transition table
	!**************************************************************************************************************************
	subroutine saveTransitionRates()
		use comparams
		use input_class
		use smallFunctions
        
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
				
		10 FORMAT (E16.8)
       
    	return    
	end subroutine saveTransitionRates
end module transitionTable_module