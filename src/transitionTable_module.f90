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

	integer, public :: partition_function_type
	
	
contains
	!**************************************************************************************************************************
	! calculate transition table
	!**************************************************************************************************************************
	
	subroutine calculateTransitionTable (cnt1,cnt2)
		use cnt_class, only: cnt
		use comparams, only: ppLen, Temperature
		use matrix_element_mod, only: calculate_kSpaceMatrixElement, calculate_finiteGeometricMatrixElement, calculate_infiniteGeometricMatrixElement_unparallel, kSpaceMatrixElement
		use physicalConstants, only: pi, kb, hb
		use prepareForster_module, only: calculatePartitionFunction, calculateDOS
		use transition_points_mod, only: findSameEnergy, sameEnergy
		use write_log_mod, only: writeLog

		type(cnt), intent(in) :: cnt1,cnt2
		character(len=200) :: logInput

		integer :: nSameEnergy, iC
		integer :: ix1, ix2, iKcm1, iKcm2
		real*8 :: partitionFunction1, partitionFunction2
		real*8 :: dos1, dos2
		complex*16 :: matrixElement, geometricMatrixElement
		
		
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
		call findSameEnergy(cnt1,cnt2)
		
		call calculate_kSpaceMatrixElement()
			
		!allocate the transition rate table
		allocate(transitionRate(2,nTheta,nc2c))
		transitionRate = 0.d0

		call calculatePartitionFunction(partition_function_type, cnt1, partitionFunction1)
		call calculatePartitionFunction(partition_function_type, cnt2, partitionFunction2)
		
		do ic2c = 1, nc2c
			c2cDistance = c2cMin+dble(ic2c-1)*dc2c

			do iTheta = 1, nTheta				
				theta = thetaMin + dble(iTheta-1)*dTheta
				
				if ((cnt1%length .lt. huge(1.d0)) .and. (cnt2%length .lt. huge(1.d0))) then
					write(logInput,*) 'Calculating finite transition rate: iTheta=', iTheta, ', nTheta=', nTheta, 'iC2C=', ic2c, ', nC2C=', nc2c
					call writeLog(logInput)

					nSameEnergy = size(sameEnergy,1)
									
					do iC = 1,nSameEnergy
													
						ix1 = sameEnergy(iC,1)
						ix2 = sameEnergy(iC,2)
						iKcm1 = sameEnergy(iC,3)
						iKcm2 = sameEnergy(iC,4)

						call calculate_finiteGeometricMatrixElement(iKcm1, iKcm2, theta, c2cDistance, geometricMatrixElement)
						
						matrixElement = geometricMatrixElement * kSpaceMatrixElement(iC)
						call calculateDOS(cnt1,iKcm1,ix1,dos1)
						call calculateDOS(cnt2,iKcm2,ix2,dos2)

						transitionRate(1,iTheta,ic2c) = transitionRate(1,iTheta,ic2c) + exp(-(cnt1%Ex_t(ix1,iKcm1))/kb/Temperature) * (abs(matrixElement)**2) * dos2 / hb / cnt1%length / partitionFunction1
						transitionRate(2,iTheta,ic2c) = transitionRate(2,iTheta,ic2c) + exp(-(cnt2%Ex_t(ix2,iKcm2))/kb/Temperature) * (abs(matrixElement)**2) * dos1 / hb / cnt2%length / partitionFunction2

					end do

				else
					write(logInput,*) 'Calculating infinite transition rate: iTheta=', iTheta, ', nTheta=', nTheta, 'iC2C=', ic2c, ', nC2C=', nc2c
					call writeLog(logInput)
					
					if (theta .eq. 0.d0) then
						theta = thetaMin + dble(iTheta)*dTheta/2.d0
					endif

					nSameEnergy = size(sameEnergy,1)
									
					do iC = 1,nSameEnergy
						
						ix1 = sameEnergy(iC,1)
						ix2 = sameEnergy(iC,2)
						iKcm1 = sameEnergy(iC,3)
						iKcm2 = sameEnergy(iC,4)

						call calculate_infiniteGeometricMatrixElement_unparallel(iKcm1, iKcm2, theta, c2cDistance, geometricMatrixElement)
						
						matrixElement = geometricMatrixElement * kSpaceMatrixElement(iC)
						call calculateDOS(cnt1,iKcm1,ix1,dos1)
						call calculateDOS(cnt2,iKcm2,ix2,dos2)

						transitionRate(1,iTheta,ic2c) = transitionRate(1,iTheta,ic2c) + exp(-(cnt1%Ex_t(ix1,iKcm1))/kb/Temperature) * (abs(matrixElement)**2) * dos2 * (sin(theta)) / hb / ppLen/ partitionFunction1
						transitionRate(2,iTheta,ic2c) = transitionRate(2,iTheta,ic2c) + exp(-(cnt2%Ex_t(ix2,iKcm2))/kb/Temperature) * (abs(matrixElement)**2) * dos1 * (sin(theta)) / hb / ppLen/ partitionFunction2

					end do

				end if

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