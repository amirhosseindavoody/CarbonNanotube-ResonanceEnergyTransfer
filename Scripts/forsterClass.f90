module forsterClass
    use cntClass
    implicit none
    private
    
    integer, dimension(:,:), allocatable, public :: finalCrossingPoints, finalSameEnergy
		real*8, public :: totalTransitionRate12, totalTransitionRate21
    
    public  :: findCrossings, findSameEnergy, calculateDOS, calculateTransferRateParallel, calculateTransferRatePerpendicular
    
    contains
      !**************************************************************************************************************************
      ! find the points that the bands cross each other
      !**************************************************************************************************************************
      subroutine findCrossings(cnt1,cnt2)
        type(cnt), intent(in) :: cnt1,cnt2
        integer :: ix1,ix2
        integer :: iKcm
        integer :: tmp,i, nCrossing
        integer, dimension(:,:), allocatable :: crossingPoints
        real*8 :: rtmp1, rtmp2
        
        tmp = cnt1.nX * cnt2.nX * (2*cnt1.iKcm_max+1) * (2*cnt2.iKcm_max+1)
        allocate(crossingPoints(tmp,3))
        do i = 1,tmp
          crossingPoints(i,1) = 0    
          crossingPoints(i,2) = 0
          crossingPoints(i,3) = 0
        end do
        
        nCrossing = 0
        do ix1 = 1,cnt1.nX
          do ix2 = 1,cnt2.nX
            do iKcm = cnt1.iKcm_min+1 , cnt1.iKcm_max
              rtmp1 = (cnt1.Ex0_A2(ix1,iKcm)-cnt2.Ex0_A2(ix2,iKcm))
              rtmp2 = (cnt1.Ex0_A2(ix1,iKcm-1)-cnt2.Ex0_A2(ix2,iKcm-1))
              if ((rtmp1 * rtmp2) .le. 0.d0) then
                if ((abs(rtmp1) .le. abs(rtmp2)) .or. (abs(rtmp1) .eq. 0.d0)) then
                  nCrossing = nCrossing+1
                  crossingPoints(nCrossing,1) = ix1
                  crossingPoints(nCrossing,2) = ix2
                  crossingPoints(nCrossing,3) = iKcm
                else
                  nCrossing = nCrossing+1
                  crossingPoints(nCrossing,1) = ix1
                  crossingPoints(nCrossing,2) = ix2
                  crossingPoints(nCrossing,3) = iKcm-1
                endif
              end if
            end do
          end do
        end do
        
        allocate(finalCrossingPoints(nCrossing ,3))
        finalCrossingPoints(:,:) = crossingPoints(1:nCrossing,:)
        print *, "nCrossing = ", nCrossing
        
			end subroutine findCrossings
			
			!**************************************************************************************************************************
      ! find the points that the bands cross each other
      !**************************************************************************************************************************
      subroutine findSameEnergy(cnt1,cnt2)
        type(cnt), intent(in) :: cnt1,cnt2
        integer :: ix1,ix2
        integer :: iKcm1, iKcm2
        integer :: tmp,i, nSameEnergy
        integer, dimension(:,:), allocatable :: sameEnergy
        real*8 :: rtmp1, rtmp2
        
        tmp = cnt1.nX * cnt2.nX * (2*cnt1.iKcm_max+1) * (2*cnt2.iKcm_max+1)
        allocate(sameEnergy(tmp,4))
        do i = 1,tmp
          sameEnergy(i,1) = 0    
          sameEnergy(i,2) = 0
          sameEnergy(i,3) = 0
					sameEnergy(i,4) = 0
        end do
        
        nSameEnergy = 0
        do ix1 = 1,cnt1.nX
					do iKcm1 = cnt1.iKcm_min , cnt1.iKcm_max
						do ix2 = 1,cnt2.nX
							do iKcm2 = cnt2.iKcm_min+1 , cnt2.iKcm_max
								rtmp1 = (cnt1.Ex0_A2(ix1,iKcm1)-cnt2.Ex0_A2(ix2,iKcm2))
								rtmp2 = (cnt1.Ex0_A2(ix1,iKcm1)-cnt2.Ex0_A2(ix2,iKcm2-1))
								if ((rtmp1 * rtmp2) .le. 0.d0) then
									if ((abs(rtmp1) .le. abs(rtmp2)) .or. (abs(rtmp1) .eq. 0.d0)) then
										nSameEnergy = nSameEnergy+1
										sameEnergy(nSameEnergy,1) = ix1
										sameEnergy(nSameEnergy,2) = ix2
										sameEnergy(nSameEnergy,3) = iKcm1
										sameEnergy(nSameEnergy,4) = iKcm2
									else
										nSameEnergy = nSameEnergy+1
										sameEnergy(nSameEnergy,1) = ix1
										sameEnergy(nSameEnergy,2) = ix2
										sameEnergy(nSameEnergy,3) = iKcm1
										sameEnergy(nSameEnergy,4) = iKcm2-1
									endif
								end if
							end do
					  end do
					end do
				end do
				
        
        allocate(finalSameEnergy(nSameEnergy ,4))
        finalSameEnergy(:,:) = sameEnergy(1:nSameEnergy,:)
        print *, "nSameEnergy = ", nSameEnergy
        
      end subroutine findSameEnergy
      
      
      !**************************************************************************************************************************
      ! calculate the partition function for a given carbon nanotube
      !**************************************************************************************************************************
      subroutine calculatePartitionFunction(currcnt, partitionFunction)
        use physicalConstants, only : kb
        use inputParameters, only : Temperature
        type(cnt), intent(in) :: currcnt
				real*8, intent(out) :: partitionFunction
        integer :: ix, iKcm
        
        partitionFunction = 0.d0
        
        do ix = 1,currcnt.nX
          do iKcm = currcnt.iKcm_min,currcnt.iKcm_max
            partitionFunction = partitionFunction + currcnt.dk * exp(-currcnt.Ex0_A2(ix,iKcm)/kb/Temperature)    
          end do
        end do
        
      end subroutine calculatePartitionFunction
      
      !**************************************************************************************************************************
      ! calculate the density of states at a given point
      !**************************************************************************************************************************
      subroutine calculateDOS(currcnt,iKcm,iX,DOS)
        type(cnt), intent(in) :: currcnt
        integer, intent(in) :: iKcm,iX
        real*8, intent(out) :: DOS
        
        if (iKcm == currcnt.iKcm_min) then
          DOS = currcnt.dk/abs(currcnt.Ex0_A2(iX,iKcm)-currcnt.Ex0_A2(iX,iKcm+1))
				else if(iKcm == currcnt.iKcm_max) then
					DOS = currcnt.dk/abs(currcnt.Ex0_A2(iX,iKcm-1)-currcnt.Ex0_A2(iX,iKcm))
				else if(iKcm == 0) then
					DOS = currcnt.dk/abs(currcnt.Ex0_A2(iX,iKcm-1)-currcnt.Ex0_A2(iX,iKcm))
				else
					DOS = 2.d0*currcnt.dk/abs(currcnt.Ex0_A2(iX,iKcm-1)-currcnt.Ex0_A2(iX,iKcm+1))
				end if
			end subroutine calculateDOS
			
			!**************************************************************************************************************************
      ! calculate the matrix element for the crossing point number iC in two parallel tube
      !**************************************************************************************************************************
      subroutine calculateMatrixElementParallel(cnt1,cnt2,iC,matrixElementFinal)
        use inputParameters
        use physicalConstants, only : i1, pi, eps0, q0
				use smallFunctions
        type(cnt), intent(in) :: cnt1,cnt2
        integer, intent(in) :: iC
				complex*16, intent(out) :: matrixElementFinal
        integer :: ix1,ix2
        integer :: iKcm
        integer :: ikr1, ikr2
        integer :: imu1,imu2
        integer :: is,isp
				real*8 :: dk
				complex*16, dimension(2,2) :: matrixElementTemp
				real*8 :: Ck
				real*8 :: phi1, phi2, dPhi
				integer :: iPhi1, iPhi2, nPhi
				real*8 :: radius1, radius2
				real*8 :: tmpr
				
        
        ix1 = finalCrossingPoints(iC,1)
        ix2 = finalCrossingPoints(iC,2)
        iKcm = finalCrossingPoints(iC,3)
				
				radius1 = cnt1.radius
				radius2	= cnt2.radius
				dk = cnt1.dk
				
				matrixElementFinal = (0.d0,0.d0)
				matrixElementTemp(1,1) = (0.d0,0.d0)
				matrixElementTemp(1,2) = (0.d0,0.d0)
				matrixElementTemp(2,1) = (0.d0,0.d0)
				matrixElementTemp(2,2) = (0.d0,0.d0)
        
				nPhi = 1000
				Ck = 0.d0
				dPhi = 2.d0*pi/dble(nPhi)
				do iPhi1 = 0,nPhi
					do iPhi2 = 0,nPhi
						phi1=dble(iPhi1)*dPhi
						phi2=dble(iPhi2)*dPhi
						tmpr = 2.d0*dble(iKcm)*dk*sqrt((radius1*sin(phi1)-radius2*sin(phi2))**2+(c2cDistance-radius1*cos(phi1)+radius2*cos(phi2))**2)
						Ck = Ck + sqrt(2/pi) * dPhi * dPhi * bessk0(abs(tmpr))
					end do
				end do
				
				do ikr1 = cnt1.ikr_low, cnt1.ikr_high
					do ikr2 = cnt2.ikr_low, cnt2.ikr_high					
						do is = 1,2
							do isp = 1,2
								matrixElementTemp(1,1) = matrixElementTemp(1,1) +	conjg(cnt1.Cc(1,ikr1+iKcm,is))*cnt1.Cv(1,ikr1-iKcm,is)*cnt2.Cc(1,ikr2+iKcm,isp)*conjg(cnt2.Cv(1,ikr2-iKcm,isp))
								matrixElementTemp(1,2) = matrixElementTemp(1,2) + conjg(cnt1.Cc(1,ikr1+iKcm,is))*cnt1.Cv(1,ikr1-iKcm,is)*cnt2.Cc(2,-ikr2+iKcm,isp)*conjg(cnt2.Cv(2,-ikr2-iKcm,isp))
								matrixElementTemp(2,1) = matrixElementTemp(2,1) + conjg(cnt1.Cc(2,-ikr1+iKcm,is))*cnt1.Cv(2,-ikr1-iKcm,is)*cnt2.Cc(1,ikr2+iKcm,isp)*conjg(cnt2.Cv(1,ikr2-iKcm,isp))
								matrixElementTemp(2,2) = matrixElementTemp(2,2) + conjg(cnt1.Cc(2,-ikr1+iKcm,is))*cnt1.Cv(2,-ikr1-iKcm,is)*cnt2.Cc(2,-ikr2+iKcm,isp)*conjg(cnt2.Cv(2,-ikr2-iKcm,isp))
              end do  
						end do
						matrixElementFinal = matrixElementFinal + (matrixElementTemp(1,1) + matrixElementTemp(1,2) + matrixElementTemp(2,1) + matrixElementTemp(2,2)) &
																											* conjg(cnt1.Psi0_A2(ikr1,ix1,iKcm))*cnt2.Psi0_A2(ikr2,ix2,iKcm) / (2.d0,0.d0)
						matrixElementTemp(1,1) = (0.d0,0.d0)
						matrixElementTemp(1,2) = (0.d0,0.d0)
						matrixElementTemp(2,1) = (0.d0,0.d0)
						matrixElementTemp(2,2) = (0.d0,0.d0)
					end do
				end do				
				
        matrixElementFinal = matrixElementFinal * dcmplx(q0*q0 * Ck / (4.d0*pi*eps0*(4.d0*pi**2)*2.d0*pi/dk))
        
      end subroutine calculateMatrixElementParallel
			
			!**************************************************************************************************************************
      ! calculate scattering rate from cnt1 to cnt2 when they are parallel
      !**************************************************************************************************************************
      subroutine calculateTransferRateParallel(cnt1, cnt2)
				use physicalConstants, only : kb, pi, hb
				use inputParameters, only : Temperature
        type(cnt), intent(in) :: cnt1, cnt2
				integer :: iC
				integer :: ix1, ix2, iKcm
				integer :: nCrossing
				real*8 :: partitionFunction1, partitionFunction2
				real*8 :: dos1, dos2
				complex*16 :: matrixElement
				
        call calculatePartitionFunction(cnt1, partitionFunction1)
				call calculatePartitionFunction(cnt2, partitionFunction2)
				
				print *, "partitionFunction1 = ", partitionFunction1
				print *, "partitionFunction2 = ", partitionFunction2
				
				!pause
				!stop
				
				
				totalTransitionRate12 = 0.d0
				totalTransitionRate21 = 0.d0
				nCrossing = size(finalCrossingPoints,1)
				
				do iC = 1,nCrossing
					
					ix1 = finalCrossingPoints(iC,1)
					ix2 = finalCrossingPoints(iC,2)
					iKcm = finalCrossingPoints(iC,3)
					
					if (iKcm .ne. 0) then
					
						call calculateMatrixElementParallel(cnt1,cnt2,iC,matrixElement)
						call calculateDOS(cnt1,iKcm,ix1,dos1)
						call calculateDOS(cnt2,iKcm,ix2,dos2)

						totalTransitionRate12 = totalTransitionRate12 + exp(-cnt1.Ex0_A2(ix1,iKcm)/kb/Temperature) * dble(conjg(matrixElement) * matrixElement) * dos1 * 2.d0 * pi / hb / partitionFunction1
						totalTransitionRate21 = totalTransitionRate21 + exp(-cnt2.Ex0_A2(ix2,iKcm)/kb/Temperature) * dble(conjg(matrixElement) * matrixElement) * dos2 * 2.d0 * pi / hb / partitionFunction2
					end if
										
				end do
				
				write(*,10) "totalTransitionRate12 = ", totalTransitionRate12
				write(*,10) "totalTransitionRate21 = ", totalTransitionRate21
				
10			FORMAT (A,E16.8)
				
			end subroutine calculateTransferRateParallel
			
			!**************************************************************************************************************************
      ! calculate the matrix element for the crossing point number iC in two perpendicular tube
      !**************************************************************************************************************************
      subroutine calculateMatrixElementPerpendicular(cnt1,cnt2,ix1, ix2, iKcm1, iKcm2 ,matrixElementFinal)
        use inputParameters
        use physicalConstants, only : i1, pi, eps0, q0
				use smallFunctions
        type(cnt), intent(in) :: cnt1,cnt2
				complex*16, intent(out) :: matrixElementFinal
        integer, intent(in) :: ix1,ix2
        integer, intent(in) :: iKcm1, iKcm2
				
        integer :: ikr1, ikr2
        integer :: imu1,imu2
        integer :: is,isp
				real*8 :: dk, K1, K2
				complex*16 :: matrixElementTemp
				complex*16 :: Jk
				real*8 :: phi1, phi2, dPhi
				integer :: iPhi1, iPhi2, nPhi
				integer :: ix, iy, nx, ny
				real*8 :: x, y, dx, dy, x_max, y_max
				real*8 :: radius1, radius2
				real*8 :: deltaR, tmpr
				
				radius1 = cnt1.radius
				radius2	= cnt2.radius
				dk = cnt1.dk
				K1 = 2.d0*dble(iKcm1)*dk
				K2 = 2.d0*dble(iKcm2)*dk
        
				nPhi = 20
				dPhi = 2.d0*pi/dble(nPhi)
				
				x_max = 20.d-9
				y_max = x_max
				
				dx = 1.0d-10
				dy = dx
				
				nx = int(x_max/dx)
				ny = nx
				
				!print *, "nx = ", nx
				!print *, "nPhi = ", nPhi
				
				!Jk = 0.d0
				!do iPhi1 = 0,nPhi
				!	!print *, "iPhi1 = ",  iPhi1
				!	phi1=dble(iPhi1)*dPhi
				!	do iPhi2 = 0,nPhi
				!		phi2=dble(iPhi2)*dPhi
				!		do ix = -nx , nx
				!			x = dble(ix)*dx
				!			do iy = -ny , ny
				!				y = dble(iy)*dy								
				!				deltaR = sqrt((x-radius2*cos(phi2))**2+(y-radius1*cos(phi1))**2+(c2cDistance-radius1*sin(phi1)+radius2*sin(phi2))**2)
				!				Jk = Jk + cos(K1*x) * cos(K2*y) / deltaR 
				!			end do
				!		end do
				!	end do
				!end do
				
				Jk = (0.d0,0.d0)
				
				do iPhi1 = 0,nPhi
					phi1=dble(iPhi1)*dPhi
					do iPhi2 = 0,nPhi
						phi2=dble(iPhi2)*dPhi
						do iy = -ny , ny
							y = dble(iy)*dy
							tmpr = sqrt((y-radius1*cos(phi1))**2+(c2cDistance-radius1*sin(phi1)+radius2*sin(phi2))**2)
							Jk = Jk + exp(i1*K2*y) * exp(-i1*K1*radius2*cos(phi2)) * bessk0(abs(K1*tmpr))	
						end do
					end do
				end do
				
				Jk = Jk * dPhi * dPhi * dx * sqrt(2/pi) 
				
				matrixElementFinal = (0.d0,0.d0)
				matrixElementTemp = (0.d0,0.d0)
				
				do ikr1 = cnt1.ikr_low, cnt1.ikr_high
					do ikr2 = cnt2.ikr_low, cnt2.ikr_high					
						do is = 1,2
							do isp = 1,2
								matrixElementTemp = matrixElementTemp +	conjg(cnt1.Cc(1,ikr1+iKcm1,is))*cnt1.Cv(1,ikr1-iKcm1,is)*cnt2.Cc(1,ikr2+iKcm2,isp)*conjg(cnt2.Cv(1,ikr2-iKcm2,isp))
								matrixElementTemp = matrixElementTemp + conjg(cnt1.Cc(1,ikr1+iKcm1,is))*cnt1.Cv(1,ikr1-iKcm1,is)*cnt2.Cc(2,-ikr2+iKcm2,isp)*conjg(cnt2.Cv(2,-ikr2-iKcm2,isp))
								matrixElementTemp = matrixElementTemp + conjg(cnt1.Cc(2,-ikr1+iKcm1,is))*cnt1.Cv(2,-ikr1-iKcm1,is)*cnt2.Cc(1,ikr2+iKcm2,isp)*conjg(cnt2.Cv(1,ikr2-iKcm2,isp))
								matrixElementTemp = matrixElementTemp + conjg(cnt1.Cc(2,-ikr1+iKcm1,is))*cnt1.Cv(2,-ikr1-iKcm1,is)*cnt2.Cc(2,-ikr2+iKcm2,isp)*conjg(cnt2.Cv(2,-ikr2-iKcm2,isp))
              end do  
						end do
						matrixElementFinal = matrixElementFinal + matrixElementTemp	* conjg(cnt1.Psi0_A2(ikr1,ix1,iKcm1))*cnt2.Psi0_A2(ikr2,ix2,iKcm2) / (2.d0,0.d0)
						matrixElementTemp = (0.d0,0.d0)
					end do
				end do				
				
        matrixElementFinal = matrixElementFinal * dcmplx(q0*q0*abs(Jk) / (4.d0*pi*eps0*4.d0*pi**2*2.d0*pi/dk))
        
      end subroutine calculateMatrixElementPerpendicular
			
			!**************************************************************************************************************************
      ! calculate scattering rate from cnt1 to cnt2 when they are perpendicular
      !**************************************************************************************************************************
      subroutine calculateTransferRatePerpendicular(cnt1, cnt2)
				use physicalConstants, only : kb, pi, hb
				use inputParameters, only : Temperature, ppLen
        type(cnt), intent(in) :: cnt1, cnt2
				integer :: iC
				integer :: ix1, ix2, iKcm1, iKcm2
				integer :: nSameEnergy
				real*8 :: partitionFunction1, partitionFunction2
				real*8 :: dos1, dos2
				complex*16 :: matrixElement
				
        call calculatePartitionFunction(cnt1, partitionFunction1)
				call calculatePartitionFunction(cnt2, partitionFunction2)
				
				print *, "partitionFunction1 = ", partitionFunction1
				print *, "partitionFunction2 = ", partitionFunction2
				
				!pause
				!stop
				
				
				totalTransitionRate12 = 0.d0
				totalTransitionRate21 = 0.d0
				nSameEnergy = size(finalSameEnergy,1)
				
				print *, "nSameEnergy = ", nSameEnergy
				
				pause
				
				do iC = 1,nSameEnergy
					
					ix1 = finalSameEnergy(iC,1)
					ix2 = finalSameEnergy(iC,2)
					iKcm1 = finalSameEnergy(iC,3)
					iKcm2 = finalSameEnergy(iC,4)
					
					print *, "iC = ", iC
					
					if ((iKcm1 .ne. 0) .and. (iKcm2 .ne. 0)) then
					
						call calculateMatrixElementPerpendicular(cnt1,cnt2,ix1, ix2, iKcm1, iKcm2, matrixElement)
						call calculateDOS(cnt1,iKcm1,ix1,dos1)
						call calculateDOS(cnt2,iKcm2,ix2,dos2)

						totalTransitionRate12 = totalTransitionRate12 + exp(-cnt1.Ex0_A2(ix1,iKcm1)/kb/Temperature) * dble(conjg(matrixElement) * matrixElement) * dos2 * 2.d0 * pi / hb / ppLen/ partitionFunction1
						totalTransitionRate21 = totalTransitionRate21 + exp(-cnt2.Ex0_A2(ix1,iKcm1)/kb/Temperature) * dble(conjg(matrixElement) * matrixElement) * dos1 * 2.d0 * pi / hb / ppLen/ partitionFunction2
					end if
										
				end do
				
				write(*,10) "totalTransitionRate12 = ", totalTransitionRate12
				write(*,10) "totalTransitionRate21 = ", totalTransitionRate21
				
10			FORMAT (A,E16.8)
				
      end subroutine calculateTransferRatePerpendicular
      
    
end module forsterClass