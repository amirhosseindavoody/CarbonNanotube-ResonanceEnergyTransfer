module geometric_matrix_element_mod
	implicit none
	private
    public  :: calculate_finite_geometric_matrix_element, calculate_infinite_geometric_matrix_element
    
contains

	!**************************************************************************************************************************
	! calculate the geometric part of matrix element for two finite tubes forming angle theta
	!**************************************************************************************************************************
	
	subroutine calculate_finite_geometric_matrix_element(iKcm1, iKcm2, theta, c2cDistance, geometricMatrixElement)
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
						
		return		
	end subroutine calculate_finite_geometric_matrix_element


	!**************************************************************************************************************************
	! calculate the geometric part of matrix element for two infinite tubes forming angle theta
	!**************************************************************************************************************************
	
	subroutine calculate_infinite_geometric_matrix_element(iKcm1, iKcm2, theta, c2cDistance, geometricMatrixElement)
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
	end subroutine calculate_infinite_geometric_matrix_element

				
end module geometric_matrix_element_mod