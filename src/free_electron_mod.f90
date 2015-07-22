!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Declaration of input parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
module free_electron_mod
	implicit none
	private
	public :: calculate_free_electron_transition, save_free_electron_transition

contains	
	!**************************************************************************************************************************
	! calculate the free-electron transition energies between subbands of CNTs
	!**************************************************************************************************************************

	subroutine calculate_free_electron_transition(currcnt)
		use cnt_class, only: cnt
		use write_log_mod, only: writeLog

		type(cnt), intent(inout) :: currcnt
		integer :: iKcm, ikr

		call writeLog("Calculating free-electron free-hole transition energies")

		allocate(currcnt%E_free_eh(2, 2, currcnt%ikr_low:currcnt%ikr_high, currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))

		do iKcm=currcnt%iKcm_min_fine,currcnt%iKcm_max_fine
			do ikr=currcnt%ikr_low,currcnt%ikr_high
				currcnt%E_free_eh(1,1,ikr,iKcm) = (currcnt%Ek(1,ikr*currcnt%dk_dkx_ratio+iKcm,1)+currcnt%Sk(1,ikr*currcnt%dk_dkx_ratio+iKcm,1)) - (currcnt%Ek(1,ikr*currcnt%dk_dkx_ratio-iKcm,2)+currcnt%Sk(1,ikr*currcnt%dk_dkx_ratio-iKcm,2))
				currcnt%E_free_eh(1,2,ikr,iKcm) = (currcnt%Ek(1,ikr*currcnt%dk_dkx_ratio+iKcm,1)+currcnt%Sk(1,ikr*currcnt%dk_dkx_ratio+iKcm,1)) - (currcnt%Ek(2,ikr*currcnt%dk_dkx_ratio-iKcm,2)+currcnt%Sk(2,ikr*currcnt%dk_dkx_ratio-iKcm,2))
				currcnt%E_free_eh(2,1,ikr,iKcm) = (currcnt%Ek(2,ikr*currcnt%dk_dkx_ratio+iKcm,1)+currcnt%Sk(2,ikr*currcnt%dk_dkx_ratio+iKcm,1)) - (currcnt%Ek(1,ikr*currcnt%dk_dkx_ratio-iKcm,2)+currcnt%Sk(1,ikr*currcnt%dk_dkx_ratio-iKcm,2))
				currcnt%E_free_eh(2,2,ikr,iKcm) = (currcnt%Ek(2,ikr*currcnt%dk_dkx_ratio+iKcm,1)+currcnt%Sk(2,ikr*currcnt%dk_dkx_ratio+iKcm,1)) - (currcnt%Ek(2,ikr*currcnt%dk_dkx_ratio-iKcm,2)+currcnt%Sk(2,ikr*currcnt%dk_dkx_ratio-iKcm,2))
			enddo
		enddo
		
		return
	end subroutine calculate_free_electron_transition


	!**************************************************************************************************************************
	! save the free-electron transition energies
	!**************************************************************************************************************************
	
	subroutine save_free_electron_transition (currcnt)
		use cnt_class, only: cnt
		use write_log_mod, only: writeLog

		type(cnt), intent(inout) :: currcnt
		integer :: iKcm, ikr

		open(unit=100,file='free_electron_transition_pp.dat',status="unknown")
		open(unit=101,file='free_electron_transition_pm.dat',status="unknown")
		open(unit=102,file='free_electron_transition_mp.dat',status="unknown")
		open(unit=103,file='free_electron_transition_mm.dat',status="unknown")


		do iKcm=currcnt%iKcm_min_fine,currcnt%iKcm_max_fine
			do ikr=currcnt%ikr_low,currcnt%ikr_high
				write(100,'(E16.8)', advance='no') currcnt%E_free_eh(1,1,ikr,iKcm)
				write(101,'(E16.8)', advance='no') currcnt%E_free_eh(1,2,ikr,iKcm)
				write(102,'(E16.8)', advance='no') currcnt%E_free_eh(2,1,ikr,iKcm)
				write(103,'(E16.8)', advance='no') currcnt%E_free_eh(2,2,ikr,iKcm)
			enddo
			write(100,*)
			write(101,*)
			write(102,*)
			write(103,*)
		enddo

		close(100)
		close(101)
		close(102)
		close(103)
		
	end subroutine save_free_electron_transition

	
	
end module free_electron_mod