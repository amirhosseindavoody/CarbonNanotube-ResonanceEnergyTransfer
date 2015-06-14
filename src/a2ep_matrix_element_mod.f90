module a2ep_matrix_element_mod
	implicit none
	private
    public  :: calculate_a2ep_kSpaceMatrixElement, calculate_a2ep_finitegeometricMatrixElement, calculate_a2ep_infiniteGeometricMatrixElement_unparallel
    
contains
	
	!**************************************************************************************************************************
	! calculate the k-space part of matrix element for the crossing point number iC
	!**************************************************************************************************************************
	
	subroutine calculate_a2ep_kSpaceMatrixElement(kSpaceMatrixElement)
		use comparams, only: cnt1, cnt2
		use physicalConstants, only: pi, eps0, q0, i1
		use transition_points_mod, only: sameEnergy
		use write_log_mod, only: writeLog

		complex*16, dimension(:), allocatable, intent(inout) :: kSpaceMatrixElement
		complex*16 :: tmpc, tmpc1, tmpc2
		integer :: ix1,ix2
		integer :: iKcm1, iKcm2
		integer :: ikr1, ikr2
		integer :: ikc1_p, ikc1_m, ikv1_p, ikv1_m
		integer :: ikc2, ikv2
		integer :: is,isp
		integer :: iC
		integer :: nSameEnergy
		character(len=200) :: logInput
		real*8, dimension(2) :: Kcm1, Kcm2
		real*8 , dimension(2,2) :: ds1, ds2 ! this are relative displacement of carbon atoms in graphene unit cell

! 		write(*,*) "Calculating A to Ep k-space matrix element"
! 		read(*,*)

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

			Kcm1 = dble(cnt1%mu_cm) * cnt1%K1 + dble(iKcm1) * cnt1%dkx * cnt1%K2
			Kcm2 = dble(cnt2%mu_cm) * cnt2%K1 + dble(iKcm2) * cnt2%dkx * cnt2%K2

			tmpc = (0.d0,0.d0)
			do ikr1 = cnt1%ikr_low, cnt1%ikr_high
				ikc1_p = +ikr1 * cnt1%dk_dkx_ratio + iKcm1
				ikv1_p = +ikr1 * cnt1%dk_dkx_ratio - iKcm1
				ikc1_m = -ikr1 * cnt1%dk_dkx_ratio + iKcm1
				ikv1_m = -ikr1 * cnt1%dk_dkx_ratio - iKcm1
				do ikr2 = cnt2%ikr_low, cnt2%ikr_high
					ikc2 = +ikr2 * cnt2%dk_dkx_ratio + iKcm2
					ikv2 = +ikr2 * cnt2%dk_dkx_ratio - iKcm2
					do is = 1,2
						do isp = 1,2
							tmpc1 = conjg(cnt1%Cc(1,ikc1_p,is))*cnt1%Cv(1,ikv1_p,is)*cnt2%Cc(1,ikc2,isp)*conjg(cnt2%Cv(2,ikv2,isp))
							tmpc2 = conjg(cnt1%Cc(2,ikc1_m,is))*cnt1%Cv(2,ikv1_m,is)*cnt2%Cc(1,ikc2,isp)*conjg(cnt2%Cv(2,ikv2,isp))*dcmplx(cnt1%ex_symmetry)
							tmpc = tmpc + (tmpc1 + tmpc2)*exp(i1*dcmplx(-2.d0* dot_product(Kcm1,ds1(is,:)) + 2.d0*dot_product(Kcm2,ds2(isp,:))))
						end do  
					end do
					kSpaceMatrixElement(iC) = kSpaceMatrixElement(iC) + tmpc*conjg(cnt1%Psi_t(ikr1,ix1,iKcm1))*cnt2%Psi_t(ikr2,ix2,iKcm2)/dcmplx(sqrt(2.d0),0.d0)
					tmpc = (0.d0,0.d0)
				end do
			end do
		end do

		kSpaceMatrixElement = kSpaceMatrixElement * dcmplx(q0**2/(4.d0*pi*eps0*4.d0*(pi*pi)*sqrt(2.d0*pi/cnt1%dk * 2.d0*pi/cnt2%dk)))
			
		return		
	end subroutine calculate_a2ep_kSpaceMatrixElement

	!**************************************************************************************************************************
	! calculate the geometric part of matrix element for two finite tubes forming angle theta
	!**************************************************************************************************************************
	
	subroutine calculate_a2ep_finiteGeometricMatrixElement(iKcm1, iKcm2, theta, c2cDistance, geometricMatrixElement)
		use comparams, only: cnt1, cnt2
		use physicalConstants, only: i1, pi

		real*8, intent(in) :: theta
		real*8, intent(in) :: c2cDistance
		integer, intent(in) :: iKcm1, iKcm2
		complex*16, intent(out) :: geometricMatrixElement
		real*8 :: dx
		real*8 :: x1, x2, xp1, xp2
		real*8 :: K1, K2
		integer :: ix1, ix2
		integer :: nx1, nx2
		integer :: iPhi1, iPhi2
		integer :: nPhi1, nPhi2
		real*8 :: dPhi1, dPhi2
		real*8, dimension(:), allocatable :: xvec1, xvec2
		real*8, dimension(:), allocatable :: phi1, phi2
		complex*16, dimension(:), allocatable :: integrand

		real*8 :: radius1, radius2

! 		write(*,*) "Calculating A to Ep finite length geometric matrix element"
! 		read(*,*)

		radius1 = cnt1%radius
		radius2	= cnt2%radius

		nPhi1 = 20
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



		dx = 5.0d-10

		x1 = cnt1%center_position - cnt1%length/2.d0
		x2 = cnt1%center_position + cnt1%length/2.d0
		xp1 = cnt2%center_position - cnt2%length/2.d0
		xp2 = cnt2%center_position + cnt2%length/2.d0

		nx1 = nint((x2-x1)/dx)
		nx2 = nint((xp2-xp1)/dx)

		allocate(xvec1(nx1))
		do ix1 = 1, nx1
			xvec1(ix1) = x1+dble(ix1-1)*dx
		end do

		allocate(xvec2(nx2))
		do ix2 = 1, nx2
			xvec2(ix2) = xp1+dble(ix2-1)*dx
		end do

		allocate(integrand(nx2))

		K2 = dble(iKcm2)*cnt2%dkx
		K1 = dble(iKcm1)*cnt1%dkx
		geometricMatrixElement = (0.d0, 0.d0)
		do ix1 = 1, nx1
			integrand = dcmplx(0.d0, 0.d0)
			do iPhi1 = 1, nPhi1
				do iPhi2 = 1, nPhi2
					integrand = integrand + exp(-i1*dcmplx(2.d0*(K1*xvec1(ix1)+dble(cnt1%mu_cm)*phi1(iPhi1))))*exp(i1*dcmplx(2.d0*(K2*xvec2+dble(cnt2%mu_cm)*phi2(iPhi2))))/dcmplx(sqrt((xvec2*cos(theta)-radius2*cos(phi2(iPhi2))*sin(theta)-xvec1(ix1))**2+(xvec2*sin(theta)+radius2*cos(phi2(iPhi2))*cos(theta)-radius1*cos(phi1(iPhi1)))**2+(c2cDistance+radius2*sin(phi2(iPhi2))-radius1*sin(phi1(iPhi1)))**2))
				enddo
			enddo
			geometricMatrixElement = geometricMatrixElement + sum(integrand)
		end do
		geometricMatrixElement = geometricMatrixElement * dcmplx(dPhi1*dPhi2*(dx**2))
						
		return		
	end subroutine calculate_a2ep_finiteGeometricMatrixElement


	!**************************************************************************************************************************
	! calculate the geometric part of matrix element for two infinite tubes forming angle theta
	!**************************************************************************************************************************
	
	subroutine calculate_a2ep_infiniteGeometricMatrixElement_unparallel(iKcm1, iKcm2, theta, c2cDistance, geometricMatrixElement)
		use comparams, only: cnt1, cnt2
		use physicalConstants, only: i1, pi

		real*8, intent(in) :: theta
		real*8, intent(in) :: c2cDistance
		integer, intent(in) :: iKcm1, iKcm2
		complex*16, intent(out) :: geometricMatrixElement

		real*8 :: K1, K2
		integer :: iPhi1, iPhi2
		integer :: nPhi1, nPhi2
		real*8 :: dPhi1, dPhi2
		real*8, dimension(:), allocatable :: phi1, phi2

		real*8 :: arg1, arg2, arg3

		real*8 :: radius1, radius2


		radius1 = cnt1%radius
		radius2	= cnt2%radius

		nPhi1 = 7
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

! 		write(*,*) "Calculating A to Ep infinite length geometric matrix element"
! 		write(*,*) "cnt1%mu_cm = ", cnt1%mu_cm
! 		write(*,*) "cnt2%mu_cm = ", cnt2%mu_cm
! 		write(*,*) "radius1 = ", radius1
! 		write(*,*) "radius2 = ", radius2
! 		write(*,*) "c2cDistance = ", c2cDistance
! 		read(*,*)

		geometricMatrixElement = (0.d0, 0.d0)
		if ((iKcm1 .ne. 0) .and. (iKcm2 .ne. 0)) then
			arg1 = sqrt(K1**2+K2**2-2.d0*K1*K2*cos(theta))
			do iPhi1 = 1,nPhi1
				do iPhi2 = 1, nPhi2
					arg2 = 2.d0 * (K1*(radius2*cos(phi2(iPhi2))-radius1*cos(phi1(iPhi1))*cos(theta))+K2*(radius1*cos(phi1(iPhi1))-radius2*cos(phi2(iPhi2))*cos(theta))) / (sin(theta))
					arg3 = 2.d0 * abs((c2cDistance+radius2*sin(phi2(iPhi2))-radius1*sin(phi1(iPhi1)))/(sin(theta)))
					geometricMatrixElement = geometricMatrixElement + exp(i1*dcmplx(-2.d0*dble(cnt1%mu_cm)*phi1(iPhi1)+2.d0*dble(cnt2%mu_cm)*phi2(iPhi2))) * exp(i1*dcmplx(arg2)) * dcmplx(exp(-arg3 * arg1))
! 					geometricMatrixElement = geometricMatrixElement + exp(i1*dcmplx(arg2)) * dcmplx(exp(-arg3 * arg1))
				end do
			end do
			geometricMatrixElement = geometricMatrixElement * dcmplx(dPhi1*dPhi2*pi/arg1)
		end if
						
		return		
	end subroutine calculate_a2ep_infiniteGeometricMatrixElement_unparallel

				
end module a2ep_matrix_element_mod