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
	
	subroutine calculate_geometricMatrixElement()
		use comparams, only: cnt1, cnt2
		use write_log_mod, only: writeLog

		integer :: iTheta, nTheta = 1
		integer :: iKcm1, iKcm2

		call writeLog("Calculating geometric part of matrix element")

		allocate(geometricMatrixElement(cnt1%iKcm_min:cnt1%iKcm_max,cnt2%iKcm_min:cnt2%iKcm_max,nTheta))

		do iTheta = 1, nTheta
			do iKcm2 = cnt2%iKcm_min, cnt2%iKcm_max
				do iKcm1 = cnt1%iKcm_min, cnt1%iKcm_max
					geometricMatrixElement(iKcm1, iKcm2, iTheta) = (0.d0, 0.d0)
				end do
			end do
		end do
			
		return		
	end subroutine calculate_geometricMatrixElement

				
end module matrix_element_mod