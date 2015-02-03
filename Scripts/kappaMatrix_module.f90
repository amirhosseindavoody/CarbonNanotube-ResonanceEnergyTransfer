!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate kappa matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
module kappaMatrix_module
	use cntClass
  implicit none
	real*8, dimension(:,:), allocatable :: kappaMatrix
	
	private
	
	public :: calculateKappaMatrix
    
	contains
		!**************************************************************************************************************************
		! calculate kappa matrix
    !**************************************************************************************************************************
		subroutine calculateKappaMatrix (cnt1,cnt2)
			use inputParameters
			use parallelForster_module
			use arbitraryAngleForster_module
			
			
			type(cnt), intent(in) :: cnt1,cnt2
			integer :: iTheta1, iTheta2
			
			! set seperation properties
			c2cDistance = 10.d-9
  
			! set orientation properties
			thetaMax = 3.d0*pi/4.d0
			thetaMin = 0.d0
			nTheta = 4
			if (nTheta .ne. 1) then
				dTheta = (thetaMax-thetaMin)/dble(nTheta-1)
			else
				dTheta = 0.d0
			end if
	
			allocate(kappaMatrix(nTheta, nTheta))
	
			kappaMatrix = 0.d0 * kappaMatrix
	
			do iTheta1 = 1, nTheta
				do iTheta2 = 1, nTheta
	
				theta = dble(iTheta1-iTheta2)*dTheta
		
				if ((iTheta1-iTheta2) .eq. 0) then
					!call calculateParallelForsterRate(cnt1,cnt2, kappaMatrix(iTheta1,iTheta2), kappaMatrix(iTheta2,iTheta1))
				else
					call calculateArbitraryForsterRate(cnt1,cnt2, kappaMatrix(iTheta1,iTheta2), kappaMatrix(iTheta2,iTheta1))
				end if
			
				print *, 'iTheta1=', iTheta1, ', iTheta2=', iTheta2
				end do
			end do
			
			call saveKappaMatrix()
				
		end subroutine calculateKappaMatrix
	
		!**************************************************************************************************************************
		! save kappa matrix
    !**************************************************************************************************************************
		subroutine saveKappaMatrix()
			use ifport
			use inputParameters
			character*100 :: dirname
			integer(4) :: istat
			logical(4) :: result
			integer :: iTheta1, iTheta2
			
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
			open(unit=100,file='kappaMatrix.dat',status="unknown")
			do iTheta1 = 1,nTheta
				do iTheta2 = 1,nTheta
					write(100,10, advance='no') kappaMatrix(iTheta1,iTheta2)
				end do
				write(100,10)
			end do
			close(100)
  
10		FORMAT (E16.8)
        
			!change the directory back
			istat=chdir('..')
			if (istat .ne. 0) then
				print *, 'Directory did not changed!!!'
				pause
				stop
			end if
        
			end subroutine saveKappaMatrix
	
end module kappaMatrix_module