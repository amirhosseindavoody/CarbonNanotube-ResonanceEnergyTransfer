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
			use prepareForster_module
			use smallFunctions
			
			type(cnt), intent(in) :: cnt1,cnt2
			integer :: iTheta1, iTheta2
			
			write(logInput,*) "Calculating kappa matrix ..."
			call writeLog()
			
			! set seperation properties
			c2cDistance = 5.0d-9
			ppLen = 40.0d-9
  
			! set orientation properties
			nTheta = 10
			thetaMax = pi
			thetaMin = 0.d0
			if (nTheta .ne. 1) then
				dTheta = (thetaMax-thetaMin)/dble(nTheta)
			else
				dTheta = 0.d0
			end if
			
			ppLen  = ppLen*dble(nTheta-1)
	
			!! calculate kappa matrix for only one type of CNT in the system
			!allocate(kappaMatrix(nTheta, nTheta))
			!kappaMatrix = 0.d0 * kappaMatrix
			!
			!! calculate crossing points and same energy points between cnt1 and cnt2
			!call findCrossings(cnt1,cnt2)
			!call findSameEnergy(cnt1,cnt2)
			!
			!! calculate scattering between cnt1 and cnt2
			!do iTheta1 = 1, nTheta
			!	do iTheta2 = 1, nTheta
			!		theta = abs(dble(iTheta1-iTheta2))*dTheta
			!		if ((iTheta1-iTheta2) .ne. 0) then
			!			call calculateArbitraryForsterRate(cnt1,cnt2, kappaMatrix(iTheta1,iTheta2), kappaMatrix(iTheta2,iTheta1))
			!		end if
			!		print *, 'iTheta1=', iTheta1, ', iTheta2=', iTheta2
			!	end do
			!end do
			
			
			! calculate kappa matrix for two types of CNT in the system
			allocate(kappaMatrix(2*nTheta, 2*nTheta))
			kappaMatrix = 0.d0 * kappaMatrix
			
			! calculate crossing points and same energy points between cnt1 and cnt1
			call findCrossings(cnt1,cnt1)
			call findSameEnergy(cnt1,cnt1)
			
			! calculate scattering between cnt1 and cnt1
			do iTheta1 = 1, nTheta
				do iTheta2 = 1, nTheta
					theta = abs(dble(iTheta1-iTheta2))*dTheta
					if ((iTheta1-iTheta2) .ne. 0) then
						call calculateArbitraryForsterRate(cnt1,cnt1, kappaMatrix(iTheta1,iTheta2), kappaMatrix(iTheta2,iTheta1))
					end if
					write(logInput,*) 'iTheta1=', iTheta1, ', iTheta2=', iTheta2
					call writeLog()
				end do
			end do
			
			! calculate crossing points and same energy points between cnt2 and cnt2
			call findCrossings(cnt2,cnt2)
			call findSameEnergy(cnt2,cnt2)
			
			! calculate scattering between cnt2 and cnt2
			do iTheta1 = 1, nTheta
				do iTheta2 = 1, nTheta
					theta = abs(dble(iTheta1-iTheta2))*dTheta
					if ((iTheta1-iTheta2) .ne. 0) then
						call calculateArbitraryForsterRate(cnt2,cnt2, kappaMatrix(nTheta+iTheta1,nTheta+iTheta2), kappaMatrix(nTheta+iTheta2,nTheta+iTheta1))
					end if
					write(logInput,*) 'iTheta1=', iTheta1, ', iTheta2=', iTheta2
					call writeLog()
				end do
			end do
			
			! calculate crossing points and same energy points between cnt1 and cnt2
			call findCrossings(cnt1,cnt2)
			call findSameEnergy(cnt1,cnt2)			
			
			! calculate scattering between cnt1 and cnt2
			do iTheta1 = 1, nTheta
				do iTheta2 = 1, nTheta
					theta = abs(dble(iTheta1-iTheta2))*dTheta
					if ((iTheta1-iTheta2) .eq. 0) then
						call calculateParallelForsterRate(cnt1,cnt2, kappaMatrix(iTheta1,nTheta+iTheta2), kappaMatrix(nTheta+iTheta2,iTheta1))
					else
						call calculateArbitraryForsterRate(cnt1,cnt2, kappaMatrix(iTheta1,nTheta+iTheta2), kappaMatrix(nTheta+iTheta2,iTheta1))
					end if
					write(logInput,*) 'iTheta1=', iTheta1, ', iTheta2=', iTheta2
					call writeLog()
				end do
			end do
			
			call saveKappaMatrix()
				
		end subroutine calculateKappaMatrix
	
		!**************************************************************************************************************************
		! save kappa matrix
    	!**************************************************************************************************************************
		subroutine saveKappaMatrix()
			use inputParameters
			use smallFunctions

			character(len=200) :: command
			integer(4) :: istat
			integer :: iTheta1, iTheta2
			integer :: nKappaMatrix
			integer :: i
			
			!create and change the directory to that of the CNT
			write(command,'("mkdir ",A100)') outputDirectory
			call execute_command_line(command,wait=.true.,exitstat=i)

			istat=chdir(outputDirectory)
			if (istat .ne. 0) then
				write(logInput,*) "Directory did not changed!!!"
				call writeLog()
				call exit()
			end if
			
			nKappaMatrix = size(kappaMatrix,1)
			
			!write transition rates to the file
			open(unit=100,file='kappaMatrix.dat',status="unknown")
			do iTheta1 = 1,nKappaMatrix
				do iTheta2 = 1,nKappaMatrix
					write(100,10, advance='no') kappaMatrix(iTheta1,iTheta2)
				end do
				write(100,10)
			end do
			close(100)
  
10			FORMAT (E16.8)
        
			!change the directory back
			istat=chdir('..')
			if (istat .ne. 0) then
				write(logInput,* ) "Directory did not changed!!!"
				call writeLog()
				call exit()
			end if
        
		end subroutine saveKappaMatrix
	
end module kappaMatrix_module