module parallelForster_module
    use cntClass
    implicit none
    private
    
    public  :: calculateParallelForsterRate
    
    contains
			!**************************************************************************************************************************
      ! calculate the matrix element for the crossing point number iC in two parallel tube
      !**************************************************************************************************************************
      subroutine calculateMatrixElement(cnt1,cnt2,ix1, ix2, iKcm1, iKcm2 ,matrixElementFinal)
        use inputParameters
        use physicalConstants, only : i1, pi, eps0, q0
				use smallFunctions
        type(cnt), intent(in) :: cnt1,cnt2
				complex*16, intent(out) :: matrixElementFinal
        integer, intent(in) :: ix1,ix2
        integer, intent(in) :: iKcm1, iKcm2
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
				
				if (iKcm1 .ne. iKcm2) then
					print *, "Momentum not conserved in parallel geometry rate calculation!"
					print *, "Press Enter to exit ..."
					pause
					stop
				end if
				
        iKcm = iKcm1
				
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
						Ck = Ck + 2.d0 * dPhi * dPhi * bessk0(abs(tmpr))
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
        
      end subroutine calculateMatrixElement
			
			!**************************************************************************************************************************
      ! calculate scattering rate from cnt1 to cnt2 when they are parallel
      !**************************************************************************************************************************
      subroutine calculateParallelForsterRate(cnt1, cnt2)
				use physicalConstants, only : kb, pi, hb
				use inputParameters, only : Temperature
				use prepareForster_module
        type(cnt), intent(in) :: cnt1, cnt2
				integer :: iC
				integer :: ix1, ix2, iKcm1, iKcm2
				integer :: nCrossing
				real*8 :: partitionFunction1, partitionFunction2
				real*8 :: dos1, dos2
				complex*16 :: matrixElement
				real*8 :: totalTransitionRate12, totalTransitionRate21
				
        call calculatePartitionFunction(cnt1, partitionFunction1)
				call calculatePartitionFunction(cnt2, partitionFunction2)
				
				!print *, "partitionFunction1 = ", partitionFunction1
				!print *, "partitionFunction2 = ", partitionFunction2
				
				!pause
				!stop
				
				
				totalTransitionRate12 = 0.d0
				totalTransitionRate21 = 0.d0
				nCrossing = size(crossingPoints,1)
				
				do iC = 1,nCrossing
					
					ix1 = crossingPoints(iC,1)
					ix2 = crossingPoints(iC,2)
					iKcm1 = crossingPoints(iC,3)
					iKcm2 = crossingPoints(iC,3)
					
					if ((iKcm1 .ne. 0) .and. (iKcm2 .ne. 0)) then
					
						call calculateMatrixElement(cnt1,cnt2,ix1, ix2, iKcm1, iKcm2 ,matrixElement)
						call calculateDOS(cnt1,iKcm1,ix1,dos1)
						call calculateDOS(cnt2,iKcm2,ix2,dos2)

						totalTransitionRate12 = totalTransitionRate12 + exp(-cnt1.Ex0_A2(ix1,iKcm1)/kb/Temperature) * dble(conjg(matrixElement) * matrixElement) * dos1 * 2.d0 * pi / hb / partitionFunction1
						totalTransitionRate21 = totalTransitionRate21 + exp(-cnt2.Ex0_A2(ix2,iKcm2)/kb/Temperature) * dble(conjg(matrixElement) * matrixElement) * dos2 * 2.d0 * pi / hb / partitionFunction2
					end if
										
				end do
				
				write(*,*)
				write(*,*) "*** Forster transfer rate calculated for PARALLEL geometry ***"
				write(*,10) " totalTransitionRate12 = ", totalTransitionRate12
				write(*,10) " totalTransitionRate21 = ", totalTransitionRate21
				write(*,*) "**************************************************************"
				
10			FORMAT (A,E16.8)
				
			end subroutine calculateParallelForsterRate
			
end module parallelForster_module