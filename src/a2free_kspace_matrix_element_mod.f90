module a2free_kspace_matrix_element_mod
	implicit none
	private
    public  :: calculate_a2free_kSpaceMatrixElement
    
contains
	
	!**************************************************************************************************************************
	! calculate the k-space part of matrix element for the crossing point number iT
	!**************************************************************************************************************************
	
	subroutine calculate_a2free_kSpaceMatrixElement(nTransitionPoints, transitionPoints, kSpaceMatrixElement)
		use comparams, only: cnt1, cnt2
		use physicalConstants, only: pi, eps0, q0, i1
		use write_log_mod, only: writeLog

		integer, intent(in) :: nTransitionPoints
		integer, dimension(nTransitionPoints,4), intent(in) :: transitionPoints
		complex*16, dimension(:), allocatable, intent(inout) :: kSpaceMatrixElement
		complex*16 :: tmpc, tmpc1, tmpc2
		integer :: ix1
		integer :: iKcm1, iKcm2
		integer :: ikr1, ikr2
		integer :: ikc1_p, ikc1_m, ikv1_p, ikv1_m
		integer :: ikc2, ikv2
		integer :: is,isp
		integer :: iT
		character(len=200) :: logInput
		real*8 :: K1, K2
		real*8 , dimension(2,2) :: ds1, ds2 ! this are relative displacement of carbon atoms in graphene unit cell

		allocate(kSpaceMatrixElement(nTransitionPoints))
		kSpaceMatrixElement = kSpaceMatrixElement * dcmplx(0.d0,0.d0)

		ds1(1,:) = 0.d0
		ds1(2,:) = cnt1%aCC_vec

		ds2(1,:) = 0.d0
		ds2(2,:) = cnt2%aCC_vec		
		
		do iT = 1,nTransitionPoints
			if (mod(iT,100) .eq. 0) then
				write(logInput, '("Calculating k-space matrix element: iT = ", I6.6, "  nTransitionPoints = ", I6.6)') iT, nTransitionPoints
				call writeLog(logInput)
			end if

			ix1 = transitionPoints(iT,1)
			ikr2 = transitionPoints(iT,2)
			iKcm1 = transitionPoints(iT,3)
			iKcm2 = transitionPoints(iT,4)
			kSpaceMatrixElement(iT) = (0.d0,0.d0)

			K1 = dble(iKcm1) * cnt1%dkx
			K2 = dble(iKcm2) * cnt2%dkx

			tmpc = (0.d0,0.d0)
			do ikr1 = cnt1%ikr_low, cnt1%ikr_high
				ikc1_p = +ikr1 * cnt1%dk_dkx_ratio + iKcm1
				ikv1_p = +ikr1 * cnt1%dk_dkx_ratio - iKcm1
				ikc1_m = -ikr1 * cnt1%dk_dkx_ratio + iKcm1
				ikv1_m = -ikr1 * cnt1%dk_dkx_ratio - iKcm1
				
				ikc2 = +ikr2 * cnt2%dk_dkx_ratio + iKcm2
				ikv2 = +ikr2 * cnt2%dk_dkx_ratio - iKcm2

				do is = 1,2
					do isp = 1,2
						tmpc1 = conjg(cnt1%Cc(1,ikc1_p,is))*cnt1%Cv(1,ikv1_p,is)*cnt2%Cc(1,ikc2,isp)*conjg(cnt2%Cv(1,ikv2,isp))
						tmpc2 = conjg(cnt1%Cc(1,ikc1_p,is))*cnt1%Cv(1,ikv1_p,is)*cnt2%Cc(1,ikc2,isp)*conjg(cnt2%Cv(1,ikv2,isp))
						tmpc = tmpc + (tmpc1 + tmpc2)*exp(dcmplx(-2.d0*K1*ds1(is,2)+2.d0*K2*ds2(isp,2))*i1)
					end do  
				end do
				kSpaceMatrixElement(iT) = kSpaceMatrixElement(iT) + tmpc*conjg(cnt1%Psi_t(ikr1,ix1,iKcm1))/dcmplx(sqrt(2.d0),0.d0)
				tmpc = (0.d0,0.d0)
			end do
		end do

		kSpaceMatrixElement = kSpaceMatrixElement * dcmplx(q0**2/(4.d0*pi*eps0*sqrt(2.d0*pi/cnt1%dk)))
			
		return		
	end subroutine calculate_a2free_kSpaceMatrixElement
	
end module a2free_kspace_matrix_element_mod