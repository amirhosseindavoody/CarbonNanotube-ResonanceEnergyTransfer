module perpendicularForster_module
    use cntClass
    implicit none
    private
    
    public  :: calculatePerpendicularForsterRate
    
    contains
			
			!**************************************************************************************************************************
      ! calculate the matrix element for the crossing point number iC in two perpendicular tube
      !**************************************************************************************************************************
      subroutine calculateMatrixElement(cnt1,cnt2,ix1, ix2, iKcm1, iKcm2 ,matrixElementFinal)
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
				
				Jk = Jk * dPhi * dPhi * dx * 2.d0 
				
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
        
      end subroutine calculateMatrixElement
			
			!**************************************************************************************************************************
      ! calculate scattering rate from cnt1 to cnt2 when they are perpendicular
      !**************************************************************************************************************************
      subroutine calculatePerpendicularForsterRate(cnt1, cnt2)
				use physicalConstants, only : kb, pi, hb
				use inputParameters, only : Temperature, ppLen
				use prepareForster_module
        type(cnt), intent(in) :: cnt1, cnt2
				integer :: iC
				integer :: ix1, ix2, iKcm1, iKcm2
				integer :: nSameEnergy
				real*8 :: partitionFunction1, partitionFunction2
				real*8 :: dos1, dos2
				complex*16 :: matrixElement
				real*8 :: totalTransitionRate12, totalTransitionRate21
				
        call calculatePartitionFunction(cnt1, partitionFunction1)
				call calculatePartitionFunction(cnt2, partitionFunction2)
				
				totalTransitionRate12 = 0.d0
				totalTransitionRate21 = 0.d0
				nSameEnergy = size(sameEnergy,1)
				
				do iC = 1,nSameEnergy
					
					ix1 = sameEnergy(iC,1)
					ix2 = sameEnergy(iC,2)
					iKcm1 = sameEnergy(iC,3)
					iKcm2 = sameEnergy(iC,4)
					
					print *, "iC = ", iC, ", maximum iC = ", nSameEnergy
					
					if ((iKcm1 .ne. 0) .and. (iKcm2 .ne. 0)) then
					
						call calculateMatrixElement(cnt1,cnt2,ix1, ix2, iKcm1, iKcm2, matrixElement)
						call calculateDOS(cnt1,iKcm1,ix1,dos1)
						call calculateDOS(cnt2,iKcm2,ix2,dos2)

						totalTransitionRate12 = totalTransitionRate12 + exp(-cnt1.Ex0_A2(ix1,iKcm1)/kb/Temperature) * dble(conjg(matrixElement) * matrixElement) * dos2 / hb / ppLen/ partitionFunction1
						totalTransitionRate21 = totalTransitionRate21 + exp(-cnt2.Ex0_A2(ix1,iKcm1)/kb/Temperature) * dble(conjg(matrixElement) * matrixElement) * dos1 / hb / ppLen/ partitionFunction2
					end if
										
				end do
				
				write(*,*) "*** Forster transfer rate calculated for PERPENDICULAR geometry ***"
				write(*,10) " totalTransitionRate12 = ", totalTransitionRate12
				write(*,10) " totalTransitionRate21 = ", totalTransitionRate21
				write(*,*) "*******************************************************************"
				
10			FORMAT (A,E16.8)
				
      end subroutine calculatePerpendicularForsterRate
      
    
end module perpendicularForster_module