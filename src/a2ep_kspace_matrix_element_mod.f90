module a2ep_kspace_matrix_element_mod
	implicit none
	private
    public  :: calculate_a2ep_kSpaceMatrixElement
    
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

		kSpaceMatrixElement = kSpaceMatrixElement * dcmplx(q0**2/(4.d0*pi*eps0*sqrt(2.d0*pi/cnt1%dk * 2.d0*pi/cnt2%dk)))
			
		return		
	end subroutine calculate_a2ep_kSpaceMatrixElement

				
end module a2ep_kspace_matrix_element_mod