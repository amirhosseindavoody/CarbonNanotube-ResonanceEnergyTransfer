module a2a_matrix_element_mod
	implicit none
	private
    public  :: calculate_a2a_kSpaceMatrixElement, calculate_a2a_finitegeometricMatrixElement, calculate_a2a_infiniteGeometricMatrixElement_unparallel
    
contains
	
	!**************************************************************************************************************************
	! calculate the k-space part of matrix element for the crossing point number iC
	!**************************************************************************************************************************
	
	subroutine calculate_a2a_kSpaceMatrixElement(kSpaceMatrixElement)
		use comparams, only: cnt1, cnt2
		use physicalConstants, only: pi, eps0, q0, i1
		use transition_points_mod, only: sameEnergy
		use write_log_mod, only: writeLog

		complex*16, dimension(:), allocatable, intent(inout) :: kSpaceMatrixElement
		complex*16 :: tmpc, tmpc1, tmpc2, tmpc3, tmpc4
		integer :: ix1,ix2
		integer :: iKcm1, iKcm2
		integer :: ikr1, ikr2
		integer :: ikc1_p, ikc2_p, ikc1_m, ikc2_m, ikv1_p, ikv2_p, ikv1_m, ikv2_m
		integer :: is,isp
		integer :: iC
		integer :: nSameEnergy
		character(len=200) :: logInput
		real*8 :: K1, K2
		real*8 , dimension(2,2) :: ds1, ds2 ! this are relative displacement of carbon atoms in graphene unit cell

		nSameEnergy = size(sameEnergy,1)

		allocate(kSpaceMatrixElement(nSameEnergy))
		kSpaceMatrixElement = kSpaceMatrixElement * dcmplx(0.d0,0.d0)

		ds1(1,:) = 0.d0
		ds1(2,:) = cnt1%aCC_vec

		ds2(1,:) = 0.d0
		ds2(2,:) = cnt2%aCC_vec		
		
		do iC = 1,nSameEnergy
			if (mod(iC,100) .eq. 0) then
				write(logInput, '("Calculating k-space matrix element: iC = ", I6.6, "  nSameEnergy = ", I6.6)') iC, nSameEnergy
				call writeLog(logInput)
			end if

			ix1 = sameEnergy(iC,1)
			ix2 = sameEnergy(iC,2)
			iKcm1 = sameEnergy(iC,3)
			iKcm2 = sameEnergy(iC,4)
			kSpaceMatrixElement(iC) = (0.d0,0.d0)

			K1 = dble(iKcm1) * cnt1%dkx
			K2 = dble(iKcm2) * cnt2%dkx

			tmpc = (0.d0,0.d0)
			do ikr1 = cnt1%ikr_low, cnt1%ikr_high
				ikc1_p = +ikr1 * cnt1%dk_dkx_ratio + iKcm1
				ikv1_p = +ikr1 * cnt1%dk_dkx_ratio - iKcm1
				ikc1_m = -ikr1 * cnt1%dk_dkx_ratio + iKcm1
				ikv1_m = -ikr1 * cnt1%dk_dkx_ratio - iKcm1
				do ikr2 = cnt2%ikr_low, cnt2%ikr_high
					ikc2_p = +ikr2 * cnt2%dk_dkx_ratio + iKcm2
					ikv2_p = +ikr2 * cnt2%dk_dkx_ratio - iKcm2
					ikc2_m = -ikr2 * cnt2%dk_dkx_ratio + iKcm2
					ikv2_m = -ikr2 * cnt2%dk_dkx_ratio - iKcm2
					do is = 1,2
						do isp = 1,2
							tmpc1 = conjg(cnt1%Cc(1,ikc1_p,is))*cnt1%Cv(1,ikv1_p,is)*cnt2%Cc(1,ikc2_p,isp)*conjg(cnt2%Cv(1,ikv2_p,isp))
							tmpc2 = conjg(cnt1%Cc(1,ikc1_p,is))*cnt1%Cv(1,ikv1_p,is)*cnt2%Cc(2,ikc2_m,isp)*conjg(cnt2%Cv(2,ikv2_m,isp))*dcmplx(cnt2%ex_symmetry)
							tmpc3 = conjg(cnt1%Cc(2,ikc1_m,is))*cnt1%Cv(2,ikv1_m,is)*cnt2%Cc(1,ikc2_p,isp)*conjg(cnt2%Cv(1,ikv2_p,isp))*dcmplx(cnt1%ex_symmetry)
							tmpc4 = conjg(cnt1%Cc(2,ikc1_m,is))*cnt1%Cv(2,ikv1_m,is)*cnt2%Cc(2,ikc2_m,isp)*conjg(cnt2%Cv(2,ikv2_m,isp))*dcmplx(cnt1%ex_symmetry*cnt2%ex_symmetry)
							tmpc = tmpc + (tmpc1 + tmpc2 + tmpc3 + tmpc4)*exp(dcmplx(-2.d0*K1*ds1(is,2)+2.d0*K2*ds2(isp,2))*i1)
						end do  
					end do
					kSpaceMatrixElement(iC) = kSpaceMatrixElement(iC) + tmpc*conjg(cnt1%Psi_t(ikr1,ix1,iKcm1))*cnt2%Psi_t(ikr2,ix2,iKcm2)/dcmplx(2.d0,0.d0)
					tmpc = (0.d0,0.d0)
				end do
			end do
		end do

		kSpaceMatrixElement = kSpaceMatrixElement * dcmplx(q0**2/(4.d0*pi*eps0*4.d0*(pi*pi)*sqrt(2.d0*pi/cnt1%dk * 2.d0*pi/cnt2%dk)))
			
		return		
	end subroutine calculate_a2a_kSpaceMatrixElement

	!**************************************************************************************************************************
	! calculate the geometric part of matrix element for two finite tubes forming angle theta
	!**************************************************************************************************************************
	
	subroutine calculate_a2a_finiteGeometricMatrixElement(iKcm1, iKcm2, theta, c2cDistance, geometricMatrixElement)
		use comparams, only: cnt1, cnt2
		use physicalConstants, only: i1, pi, A_u

		real*8, intent(in) :: theta
		real*8, intent(in) :: c2cDistance
		integer, intent(in) :: iKcm1, iKcm2
		complex*16, intent(out) :: geometricMatrixElement

		real*8 :: K1, K2
		integer :: iu

		K2 = dble(iKcm2)*cnt2%dkx
		K1 = dble(iKcm1)*cnt1%dkx
		geometricMatrixElement = (0.d0, 0.d0)
		do iu = lbound(cnt1%r_posA3,1), ubound(cnt1%r_posA3,1)
			geometricMatrixElement = geometricMatrixElement + sum(exp(-i1*dcmplx(2.d0*K1*cnt1%ur_posA3(iu,1)))*exp(i1*dcmplx(2.d0*K2*cnt2%ur_posA3(:,1)))/dcmplx(sqrt((cnt2%r_posA3(:,1)-cnt1%r_posA3(iu,1))**2+(cnt2%r_posA3(:,2)-cnt1%r_posA3(iu,2))**2+(cnt2%r_posA3(:,3)-cnt1%r_posA3(iu,3))**2)))
		end do

		geometricMatrixElement = geometricMatrixElement * dcmplx(A_u*A_u/cnt1%radius/cnt2%radius)
						
		return		
	end subroutine calculate_a2a_finiteGeometricMatrixElement


	!**************************************************************************************************************************
	! calculate the geometric part of matrix element for two infinite tubes forming angle theta
	!**************************************************************************************************************************
	
	subroutine calculate_a2a_infiniteGeometricMatrixElement_unparallel(iKcm1, iKcm2, theta, c2cDistance, geometricMatrixElement)
		use comparams, only: cnt1, cnt2
		use physicalConstants, only: i1, pi

		real*8, intent(in) :: theta
		real*8, intent(in) :: c2cDistance
		integer, intent(in) :: iKcm1, iKcm2
		complex*16, intent(out) :: geometricMatrixElement

! 		integer :: iKcm1, iKcm2

		real*8 :: K1, K2
		integer :: iPhi1, iPhi2
		integer :: nPhi1, nPhi2
		real*8 :: dPhi1, dPhi2
		real*8, dimension(:), allocatable :: phi1, phi2

		real*8 :: arg1, arg2, arg3

		real*8 :: radius1, radius2

		radius1 = cnt1%radius
		radius2	= cnt2%radius

		nPhi1 = 8
		dPhi1 = 2.d0*pi/nPhi1
		
		nPhi2 = nPhi1
		dPhi2 = 2.d0*pi/nPhi2

		allocate(phi1(nPhi1))
		do iPhi1 = 1, nPhi1
			phi1(iPhi1) = dble(iPhi1)*dPhi1
		end do

		allocate(phi2(nPhi2))
		do iPhi2 = 1, nPhi2
			phi2(iPhi2) = dble(iPhi2)*dPhi2
		end do

		K2 = dble(iKcm2)*cnt2%dkx
		K1 = dble(iKcm1)*cnt1%dkx

		geometricMatrixElement = (0.d0, 0.d0)
		if ((iKcm1 .ne. 0) .and. (iKcm2 .ne. 0)) then
			arg1 = sqrt(K1**2+K2**2-2.d0*K1*K2*cos(theta))
			do iPhi1 = 1,nPhi1
				do iPhi2 = 1,nPhi2
					arg2 = 2.d0 * (K1*(radius2*cos(phi2(iPhi2))-radius1*cos(phi1(iPhi1))*cos(theta))+K2*(radius1*cos(phi1(iPhi1))-radius2*cos(phi2(iPhi2))*cos(theta))) / (sin(theta))
					arg3 = 2.d0 * abs((c2cDistance+radius2*sin(phi2(iPhi2))-radius1*sin(phi1(iPhi1)))/(sin(theta)))
					geometricMatrixElement = geometricMatrixElement + exp(i1*dcmplx(arg2)) * dcmplx(exp(-arg3 * arg1))
				end do
			end do
			geometricMatrixElement = geometricMatrixElement * dcmplx(dPhi1*dPhi2*pi/arg1)
		end if
						
		return		
	end subroutine calculate_a2a_infiniteGeometricMatrixElement_unparallel

				
end module a2a_matrix_element_mod