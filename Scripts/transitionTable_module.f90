!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate transition table
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
module transitionTable_module
	use cntClass
  implicit none
	real*8, dimension(:,:,:), allocatable :: transitionRate
	
	private
	
	public :: calculateTransitionTable
    
	contains
		!**************************************************************************************************************************
		! calculate kappa matrix
    !**************************************************************************************************************************
		subroutine calculateTransitionTable (cnt1,cnt2)
			use inputParameters
			use parallelForster_module
			use arbitraryAngleForster_module
			
			type(cnt), intent(in) :: cnt1,cnt2
			
			! set seperation properties
			c2cMin = 30.0d-9
			c2cMax = 30.0d-9
			nc2c = 1
			if (nc2c .ne. 1) then
				dc2c = (c2cMax-c2cMin)/dble(nc2c-1)
			else
				dc2c = 0.d0
			end if
			
			! set orientation properties
			thetaMax = pi/2.d0
			thetaMin = 0.d0
			nTheta = 30
			if (nTheta .ne. 1) then
				dTheta = (thetaMax-thetaMin)/dble(nTheta-1)
			else
				dTheta = 0.d0
			end if
			
			
			
			allocate(transitionRate(2,nTheta,nc2c))
			
			do ic2c = 1, nc2c
				do iTheta = 1, nTheta
			
					c2cDistance = c2cMin+dble(ic2c-1)*dc2c
					theta = dble(iTheta-1)*dTheta
				
					if (iTheta .eq. 1) then
						call calculateParallelForsterRate(cnt1,cnt2, transitionRate(1,iTheta,ic2c), transitionRate(2,iTheta,ic2c))
					else
						call calculateArbitraryForsterRate(cnt1,cnt2, transitionRate(1,iTheta,ic2c), transitionRate(2,iTheta,ic2c))
					end if
					
					print *, 'iTheta=', iTheta, ', nTheta=', nTheta, 'iC2C=', ic2c, ', nC2C=', nc2c		
				end do
			end do
			
			call saveTransitionRates()
				
		end subroutine calculateTransitionTable
	
		!**************************************************************************************************************************
		! save total exciton density of states for a given cnt
		!**************************************************************************************************************************
		subroutine saveTransitionRates()
			use ifport
			use inputParameters
			character*100 :: dirname
			integer(4) :: istat
			logical(4) :: result
        
			!create and change the directory to that of the CNT
			write(dirname,"('ForsterRate (',I2.2,',',I2.2,') to (',I2.2,',',I2.2,')')") n_ch1, m_ch1, n_ch2, m_ch2
			result=makedirqq(dirname)
			if (result) print *,'Directory creation successful!!'
			istat=chdir(dirname)
			if (istat .ne. 0) then
				print *, 'Directory did not changed!!!'
				pause
				stop
			end if
        
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
        
			!change the directory back
			istat=chdir('..')
			if (istat .ne. 0) then
				print *, 'Directory did not changed!!!'
				pause
				stop
			end if
        
		end subroutine saveTransitionRates
	
end module transitionTable_module