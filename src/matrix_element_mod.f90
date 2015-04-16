module matrix_element_mod
	implicit none
	private
    public  :: calculate_kSpaceMatrixElement, calculate_geometricMatrixElement

	complex*16, dimension(:), allocatable, public :: kSpaceMatrixElement
	complex*16, dimension(:,:,:), allocatable, public :: geometricMatrixElement
    
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
	! calculate the geometric part of matrix element for two tubes forming angle theta
	!**************************************************************************************************************************
	
	subroutine calculate_geometricMatrixElement(thetaMin, thetaMax, nTheta, c2cDistance)
		use comparams, only: cnt1, cnt2
		use write_log_mod, only: writeLog
		use physicalConstants, only: i1, pi

		character(len=100) :: logInput

		real*8, intent(in) :: thetaMin, thetaMax
		integer, intent(in) :: nTheta
		real*8, intent(in) :: c2cDistance

		complex*16 :: tmpc

		integer :: iTheta
		integer :: iKcm1, iKcm2
		real*8 :: dTheta, theta

		real*8 :: dx
		real*8 :: x1, x2, xp1, xp2
		real*8 :: K1, K2
		integer :: ix1, ix2
		integer :: nx1, nx2
		real*8, dimension(:), allocatable :: xvec1, xvec2

		real*8 :: radius1, radius2

		call writeLog("Calculating geometric part of matrix element")

		allocate(geometricMatrixElement(cnt1%iKcm_min:cnt1%iKcm_max,cnt2%iKcm_min:cnt2%iKcm_max,nTheta))

		! set orientation properties
		if (nTheta .ne. 1) then
			dTheta = (thetaMax-thetaMin)/dble(nTheta-1)
		else
			dTheta = 0.d0
		end if

		radius1 = cnt1%radius
		radius2	= cnt2%radius

		dx = 1.0d-10

		x1 = cnt1%center_position - cnt1%length/2.d0
		x2 = cnt1%center_position + cnt1%length/2.d0
		xp1 = cnt2%center_position - cnt2%length/2.d0
		xp2 = cnt2%center_position + cnt2%length/2.d0

		nx1 = nint((x1-x2)/dx)
		nx2 = nint((xp1-xp2)/dx)

		allocate(xvec1(nx1))
		do ix1 = 1, nx1
			xvec1(ix1) = x1+dble(ix1-1)*dx
		end do

		allocate(xvec2(nx2))
		do ix2 = 1, nx2
			xvec2(ix2) = xp1+dble(ix2-1)*dx
		end do

		do iTheta = 1, nTheta
			
			theta = thetaMin + dble(iTheta-1)*dTheta

			do iKcm2 = cnt2%iKcm_min, cnt2%iKcm_max
				K2 = dble(iKcm2)*cnt2%dk
				do iKcm1 = cnt1%iKcm_min, cnt1%iKcm_max
					K1 = dble(iKcm1)*cnt1%dk
					geometricMatrixElement(iKcm1, iKcm2, iTheta) = (0.d0, 0.d0)
! 					tmpc = (0.d0, 0.d0)
					do ix1 = 1, nx1
						geometricMatrixElement(iKcm1, iKcm2, iTheta) = geometricMatrixElement(iKcm1, iKcm2, iTheta) + sum(exp(-i1*dcmplx(2.d0*K1*xvec1(ix1)))*exp(i1*dcmplx(2.d0*K2*xvec2))/dcmplx(sqrt((xvec2*cos(theta)-xvec1(ix1))**2+(xvec2*sin(theta))**2+c2cDistance**2)))
					end do
				end do
			end do
			write(logInput, '("Calculating geometric matrix elemetn: theta = ", I3.3, ", matrix elemetn = (", E8.3, ",",E8.3,")" )') nint(theta*180/pi), real(geometricMatrixElement(iKcm1, iKcm2, iTheta)), aimag(geometricMatrixElement(iKcm1, iKcm2, iTheta))
			call writeLog(logInput)
		end do
			
		return		
	end subroutine calculate_geometricMatrixElement

				
end module matrix_element_mod