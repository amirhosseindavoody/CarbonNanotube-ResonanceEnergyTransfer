module matrix_element_mod
	implicit none
	private
    public  :: calculate_kSpaceMatrixElement, calculate_finitegeometricMatrixElement, calculate_infinitegeometricMatrixElement

	complex*16, dimension(:), allocatable, public :: kSpaceMatrixElement
	complex*16, dimension(:,:), allocatable, public :: geometricMatrixElement
    
contains
	
	!**************************************************************************************************************************
	! calculate the k-space part of matrix element for the crossing point number iC
	!**************************************************************************************************************************
	
	subroutine calculate_kSpaceMatrixElement()
		use comparams, only: cnt1, cnt2
		use transition_points_mod, only: sameEnergy
		use write_log_mod, only: writeLog

		complex*16 :: tmpc
		integer :: ix1,ix2
		integer :: iKcm1, iKcm2
		integer :: ikr1, ikr2
		integer :: is,isp
		integer :: iC
		integer :: nSameEnergy
		character(len=100) :: logInput

		nSameEnergy = size(sameEnergy)

		allocate(kSpaceMatrixElement(nSameEnergy))
		kSpaceMatrixElement = kSpaceMatrixElement * dcmplx(0.d0,0.d0)
		
		do iC = 1,nSameEnergy
			if (mod(iC,100) .eq. 0) then
				write(logInput, '("Calculating k-space matrix element: iC = ", I5.5, "  nSameEnergy = ", I5.5)') iC, nSameEnergy
				call writeLog(logInput)
			end if

			ix1 = sameEnergy(iC,1)
			ix2 = sameEnergy(iC,2)
			iKcm1 = sameEnergy(iC,3)
			iKcm2 = sameEnergy(iC,4)
			kSpaceMatrixElement(iC) = (0.d0,0.d0)

			if ((iKcm1 .ne. 0) .and. (iKcm2 .ne. 0)) then
				tmpc = (0.d0,0.d0)
				do ikr1 = cnt1%ikr_low, cnt1%ikr_high
					do ikr2 = cnt2%ikr_low, cnt2%ikr_high					
						do is = 1,2
							do isp = 1,2
								tmpc = tmpc + conjg(cnt1%Cc(1,ikr1+iKcm1,is))*cnt1%Cv(1,ikr1-iKcm1,is)*cnt2%Cc(1,ikr2+iKcm2,isp)*conjg(cnt2%Cv(1,ikr2-iKcm2,isp))
								tmpc = tmpc + conjg(cnt1%Cc(1,ikr1+iKcm1,is))*cnt1%Cv(1,ikr1-iKcm1,is)*cnt2%Cc(2,-ikr2+iKcm2,isp)*conjg(cnt2%Cv(2,-ikr2-iKcm2,isp))
								tmpc = tmpc + conjg(cnt1%Cc(2,-ikr1+iKcm1,is))*cnt1%Cv(2,-ikr1-iKcm1,is)*cnt2%Cc(1,ikr2+iKcm2,isp)*conjg(cnt2%Cv(1,ikr2-iKcm2,isp))
								tmpc = tmpc + conjg(cnt1%Cc(2,-ikr1+iKcm1,is))*cnt1%Cv(2,-ikr1-iKcm1,is)*cnt2%Cc(2,-ikr2+iKcm2,isp)*conjg(cnt2%Cv(2,-ikr2-iKcm2,isp))
							end do  
						end do
						kSpaceMatrixElement(iC) = kSpaceMatrixElement(iC) + tmpc*conjg(cnt1%Psi_t(ikr1,ix1,iKcm1))*cnt2%Psi_t(ikr2,ix2,iKcm2)/dcmplx(2.d0,0.d0)
						tmpc = (0.d0,0.d0)
					end do
				end do
			end if
		end do
			
		return		
	end subroutine calculate_kSpaceMatrixElement

	!**************************************************************************************************************************
	! calculate the geometric part of matrix element for two finite tubes forming angle theta
	!**************************************************************************************************************************
	
	subroutine calculate_finiteGeometricMatrixElement(theta, c2cDistance)
		use comparams, only: cnt1, cnt2
		use physicalConstants, only: i1, pi

		real*8, intent(in) :: theta
		real*8, intent(in) :: c2cDistance

		integer :: iKcm1, iKcm2

		real*8 :: dx
		real*8 :: x1, x2, xp1, xp2
		real*8 :: K1, K2
		integer :: ix1, ix2
		integer :: nx1, nx2
		real*8, dimension(:), allocatable :: xvec1, xvec2

		real*8 :: radius1, radius2

		if (.not. allocated(geometricMatrixElement)) allocate(geometricMatrixElement(cnt1%iKcm_min:cnt1%iKcm_max,cnt2%iKcm_min:cnt2%iKcm_max))

		radius1 = cnt1%radius
		radius2	= cnt2%radius

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

! 		if (theta .eq. 0.d0*pi/180.d0) open(unit=100,file='geometricMatrixElement_0.dat',status="unknown")
! 		if (theta .eq. 45.d0*pi/180.d0) open(unit=100,file='geometricMatrixElement_45.dat',status="unknown")
! 		if (theta .eq. 90.d0*pi/180.d0) open(unit=100,file='geometricMatrixElement_90.dat',status="unknown")

		do iKcm2 = cnt2%iKcm_min, cnt2%iKcm_max
			K2 = dble(iKcm2)*cnt2%dk
			do iKcm1 = cnt1%iKcm_min, cnt1%iKcm_max
				K1 = dble(iKcm1)*cnt1%dk
				geometricMatrixElement(iKcm1, iKcm2) = (0.d0, 0.d0)
				do ix1 = 1, nx1
					geometricMatrixElement(iKcm1, iKcm2) = geometricMatrixElement(iKcm1, iKcm2) + sum(exp(-i1*dcmplx(2.d0*K1*xvec1(ix1)))*exp(i1*dcmplx(2.d0*K2*xvec2))/dcmplx(sqrt((xvec2*cos(theta)-xvec1(ix1))**2+(xvec2*sin(theta))**2+c2cDistance**2)))
				end do
				geometricMatrixElement(iKcm1, iKcm2) = geometricMatrixElement(iKcm1, iKcm2) * dcmplx((2.d0*pi*dx)**2)

! 				write(100,'(SP , E16.3 )', advance='no') abs(geometricMatrixElement(iKcm1, iKcm2))
			end do
! 			write(100, *)
		end do
! 		close(100)
						
		return		
	end subroutine calculate_finiteGeometricMatrixElement


	!**************************************************************************************************************************
	! calculate the geometric part of matrix element for two infinite tubes forming angle theta
	!**************************************************************************************************************************
	
	subroutine calculate_infiniteGeometricMatrixElement(theta, c2cDistance)
		use comparams, only: cnt1, cnt2
		use physicalConstants, only: i1, pi

		real*8, intent(in) :: theta
		real*8, intent(in) :: c2cDistance

		integer :: iKcm1, iKcm2

		real*8 :: dx
		real*8 :: x1, x2, xp1, xp2
		real*8 :: K1, K2
		integer :: ix1, ix2
		integer :: nx1, nx2
		real*8, dimension(:), allocatable :: xvec1, xvec2

		real*8 :: radius1, radius2

		if (.not. allocated(geometricMatrixElement)) allocate(geometricMatrixElement(cnt1%iKcm_min:cnt1%iKcm_max,cnt2%iKcm_min:cnt2%iKcm_max))

		radius1 = cnt1%radius
		radius2	= cnt2%radius

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

! 		if (theta .eq. 0.d0*pi/180.d0) open(unit=100,file='geometricMatrixElement_0.dat',status="unknown")
! 		if (theta .eq. 45.d0*pi/180.d0) open(unit=100,file='geometricMatrixElement_45.dat',status="unknown")
! 		if (theta .eq. 90.d0*pi/180.d0) open(unit=100,file='geometricMatrixElement_90.dat',status="unknown")

		do iKcm2 = cnt2%iKcm_min, cnt2%iKcm_max
			K2 = dble(iKcm2)*cnt2%dk
			do iKcm1 = cnt1%iKcm_min, cnt1%iKcm_max
				K1 = dble(iKcm1)*cnt1%dk
				geometricMatrixElement(iKcm1, iKcm2) = (0.d0, 0.d0)
				do ix1 = 1, nx1
					geometricMatrixElement(iKcm1, iKcm2) = geometricMatrixElement(iKcm1, iKcm2) + sum(exp(-i1*dcmplx(2.d0*K1*xvec1(ix1)))*exp(i1*dcmplx(2.d0*K2*xvec2))/dcmplx(sqrt((xvec2*cos(theta)-xvec1(ix1))**2+(xvec2*sin(theta))**2+c2cDistance**2)))
				end do
				geometricMatrixElement(iKcm1, iKcm2) = geometricMatrixElement(iKcm1, iKcm2) * dcmplx((2.d0*pi*dx)**2)

! 				write(100,'(SP , E16.3 )', advance='no') abs(geometricMatrixElement(iKcm1, iKcm2))
			end do
! 			write(100, *)
		end do
! 		close(100)
						
		return		
	end subroutine calculate_infiniteGeometricMatrixElement

				
end module matrix_element_mod