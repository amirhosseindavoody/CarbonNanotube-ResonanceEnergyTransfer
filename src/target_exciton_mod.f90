module target_exciton_mod
	implicit none
	private
	public :: set_target_exciton
	
contains		
	subroutine set_target_exciton(currcnt, exciton_type)
		use cnt_class, only: cnt
		use write_log_mod, only: writeLog

		type(cnt), intent(inout) :: currcnt
		character(len=*), intent(in) :: exciton_type

		if (allocated(currcnt%Ex_t)) deallocate(currcnt%Ex_t)
		if (allocated(currcnt%Psi_t)) deallocate(currcnt%Psi_t)

		! set the information of the target exciton type
		select case (trim(exciton_type))
		case ('Ex_A1')
			call writeLog("currcnt target exciton: Ex_A1")
			currcnt%targetExcitonType = 'Ex_A1'
			allocate(currcnt%Ex_t(1:currcnt%nX_a,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
			allocate(currcnt%Psi_t(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%nX_a,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
			currcnt%Ex_t = currcnt%Ex_A1
			currcnt%Psi_t = currcnt%Psi_A1
			currcnt%ex_symmetry = -1.d0
		case ('Ex0_A2')
			call writeLog("currcnt target exciton: Ex0_A2")
			currcnt%targetExcitonType = 'Ex0_A2'
			allocate(currcnt%Ex_t(1:currcnt%nX_a,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
			allocate(currcnt%Psi_t(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%nX_a,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
			currcnt%Ex_t = currcnt%Ex0_A2
			currcnt%Psi_t = currcnt%Psi0_A2
			currcnt%ex_symmetry = +1.d0
		case ('Ex1_A2')
			call writeLog("currcnt target exciton: Ex1_A2")
			currcnt%targetExcitonType = 'Ex1_A2'
			allocate(currcnt%Ex_t(1:currcnt%nX_a,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
			allocate(currcnt%Psi_t(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%nX_a,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
			currcnt%Ex_t = currcnt%Ex1_A2
			currcnt%Psi_t = currcnt%Psi1_A2
			currcnt%ex_symmetry = +1.d0

		case ('Ex0_Ep')
			call writeLog("currcnt target exciton: Ex0_Ep")
			currcnt%targetExcitonType = 'Ex0_Ep'
			allocate(currcnt%Ex_t(1:currcnt%nX_a,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
			allocate(currcnt%Psi_t(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%nX_a,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
			currcnt%Ex_t = currcnt%Ex0_Ep
			currcnt%Psi_t = currcnt%Psi0_Ep
			currcnt%ex_symmetry = 0.d0
		case ('Ex0_Em')
			call writeLog("currcnt target exciton: Ex0_Em")
			currcnt%targetExcitonType = 'Ex0_Em'
			allocate(currcnt%Ex_t(1:currcnt%nX_a,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
			allocate(currcnt%Psi_t(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%nX_a,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
			currcnt%Ex_t = currcnt%Ex0_Em
			currcnt%Psi_t = currcnt%Psi0_Em
			currcnt%ex_symmetry = 0.d0
		case ('Ex1_Ep')
			call writeLog("currcnt target exciton: Ex1_Ep")
			currcnt%targetExcitonType = 'Ex1_Ep'
			allocate(currcnt%Ex_t(1:currcnt%nX_a,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
			allocate(currcnt%Psi_t(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%nX_a,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
			currcnt%Ex_t = currcnt%Ex1_Ep
			currcnt%Psi_t = currcnt%Psi1_Ep
			currcnt%ex_symmetry = 0.d0
		case ('Ex1_Em')
			call writeLog("currcnt target exciton: Ex1_Em")
			currcnt%targetExcitonType = 'Ex1_Em'
			allocate(currcnt%Ex_t(1:currcnt%nX_a,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
			allocate(currcnt%Psi_t(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%nX_a,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
			currcnt%Ex_t = currcnt%Ex1_Em
			currcnt%Psi_t = currcnt%Psi1_Em
			currcnt%ex_symmetry = 0.d0
		case default
			write(*,*) "Could not recognize exciton type!!!!"
			call exit()
		end select

		select case (trim(currcnt%targetExcitonType))
		case('Ex_A1', 'Ex0_A2', 'Ex1_A2')
			currcnt%mu_cm = 0
		case('Ex0_Ep', 'Ex1_Ep')
			currcnt%mu_cm = +1 * currcnt%min_sub(currcnt%i_sub)
		case('Ex0_Em', 'Ex1_Em')
			currcnt%mu_cm = -1 * currcnt%min_sub(currcnt%i_sub)
		case default
			write(*,*) "ERROR: undetermined target exciton type!!!!"
			call exit()
		end select

		currcnt%nX_t = size(currcnt%Ex_t,1)

		return

	end subroutine set_target_exciton

end module target_exciton_mod