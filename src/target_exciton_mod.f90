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
		real*8 :: tmpr
		integer :: iKcm, ikr, ix

		if (allocated(currcnt%Ex_t)) deallocate(currcnt%Ex_t)
		if (allocated(currcnt%Psi_t)) deallocate(currcnt%Psi_t)

		! set the information of the target exciton type
		select case (trim(exciton_type))
		case ('Ex_A1')
			call writeLog("currcnt target exciton: Ex_A1")
			currcnt%targetExcitonType = 'Ex_A1'
			
			allocate(currcnt%Ex_t(1:currcnt%nX_a,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
			allocate(currcnt%Psi_t(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%nX_a,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
			
			open(unit=100,file=trim(currcnt%directory)//'Ex_A1.dat',status="old")
			open(unit=101,file=trim(currcnt%directory)//'Psi_A1.dat',status="old")
			do iKcm=currcnt%iKcm_min_fine,currcnt%iKcm_max_fine
				do iX=1,currcnt%nX_a
					read(100,'(E16.8)', advance='no') currcnt%Ex_t(iX,iKcm)
					do ikr=currcnt%ikr_low,currcnt%ikr_high
						read(101,'(E16.8,E16.8)', advance='no') currcnt%Psi_t(ikr,iX,iKcm)
					enddo
				enddo
				
				read(100,'(E16.8)')
				read(101,'(E16.8)')
			enddo
			close(100)
			close(101)

			currcnt%ex_symmetry = -1.d0
		case ('Ex0_A2')
			call writeLog("currcnt target exciton: Ex0_A2")
			currcnt%targetExcitonType = 'Ex0_A2'
			
			allocate(currcnt%Ex_t(1:currcnt%nX_a,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
			allocate(currcnt%Psi_t(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%nX_a,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
			
			open(unit=100,file=trim(currcnt%directory)//'Ex0_A2.dat',status="old")
			open(unit=101,file=trim(currcnt%directory)//'Psi0_A2.dat',status="old")
			do iKcm=currcnt%iKcm_min_fine,currcnt%iKcm_max_fine
				do iX=1,currcnt%nX_a
					read(100,'(E16.8)', advance='no') currcnt%Ex_t(iX,iKcm)
					do ikr=currcnt%ikr_low,currcnt%ikr_high
						read(101,'(E16.8,E16.8)', advance='no') currcnt%Psi_t(ikr,iX,iKcm)
					enddo
				enddo
				
				read(100,'(E16.8)')
				read(101,'(E16.8)')
			enddo
			close(100)
			close(101)

			currcnt%ex_symmetry = +1.d0

		case ('Ex1_A2')
			call writeLog("currcnt target exciton: Ex1_A2")
			currcnt%targetExcitonType = 'Ex1_A2'
			
			allocate(currcnt%Ex_t(1:currcnt%nX_a,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
			allocate(currcnt%Psi_t(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%nX_a,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))

			open(unit=100,file=trim(currcnt%directory)//'Ex1_A2.dat',status="old")
			open(unit=101,file=trim(currcnt%directory)//'Psi1_A2.dat',status="old")
			do iKcm=currcnt%iKcm_min_fine,currcnt%iKcm_max_fine
				do iX=1,currcnt%nX_a
					read(100,'(E16.8)', advance='no') currcnt%Ex_t(iX,iKcm)
					do ikr=currcnt%ikr_low,currcnt%ikr_high
						read(101,'(E16.8,E16.8)', advance='no') currcnt%Psi_t(ikr,iX,iKcm)
					enddo
				enddo
				
				read(100,'(E16.8)')
				read(101,'(E16.8)')
			enddo
			close(100)
			close(101)

			currcnt%ex_symmetry = +1.d0

		case ('Ex0_Ep')
			call writeLog("currcnt target exciton: Ex0_Ep")
			currcnt%targetExcitonType = 'Ex0_Ep'
			
			allocate(currcnt%Ex_t(1:currcnt%nX_a,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
			allocate(currcnt%Psi_t(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%nX_a,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))

			open(unit=100,file=trim(currcnt%directory)//'Ex0_Ep.dat',status="old")
			open(unit=101,file=trim(currcnt%directory)//'Psi0_Ep.dat',status="old")
			do iKcm=currcnt%iKcm_min_fine,currcnt%iKcm_max_fine
				do iX=1,currcnt%nX_e
					read(100,'(E16.8)', advance='no') currcnt%Ex_t(iX,iKcm)
					do ikr=currcnt%ikr_low,currcnt%ikr_high
						read(101,'(E16.8,E16.8)', advance='no') currcnt%Psi_t(ikr,iX,iKcm)
					enddo
				enddo
				
				read(100,'(E16.8)')
				read(101,'(E16.8)')
			enddo
			close(100)
			close(101)

			currcnt%ex_symmetry = 0.d0
		case ('Ex0_Em')
			call writeLog("currcnt target exciton: Ex0_Em")
			currcnt%targetExcitonType = 'Ex0_Em'
			
			allocate(currcnt%Ex_t(1:currcnt%nX_e,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
			allocate(currcnt%Psi_t(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%nX_e,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
			
			open(unit=100,file=trim(currcnt%directory)//'Ex0_Em.dat',status="old")
			open(unit=101,file=trim(currcnt%directory)//'Psi0_Em.dat',status="old")
			do iKcm=currcnt%iKcm_min_fine,currcnt%iKcm_max_fine
				do iX=1,currcnt%nX_e
					read(100,'(E16.8)', advance='no') currcnt%Ex_t(iX,iKcm)
					do ikr=currcnt%ikr_low,currcnt%ikr_high
						read(101,'(E16.8,E16.8)', advance='no') currcnt%Psi_t(ikr,iX,iKcm)
					enddo
				enddo
				
				read(100,'(E16.8)')
				read(101,'(E16.8)')
			enddo
			close(100)
			close(101)

			currcnt%ex_symmetry = 0.d0
		case ('Ex1_Ep')
			call writeLog("currcnt target exciton: Ex1_Ep")
			currcnt%targetExcitonType = 'Ex1_Ep'
			
			allocate(currcnt%Ex_t(1:currcnt%nX_e,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
			allocate(currcnt%Psi_t(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%nX_e,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
			
			open(unit=100,file=trim(currcnt%directory)//'Ex1_Ep.dat',status="old")
			open(unit=101,file=trim(currcnt%directory)//'Psi1_Ep.dat',status="old")
			do iKcm=currcnt%iKcm_min_fine,currcnt%iKcm_max_fine
				do iX=1,currcnt%nX_e
					read(100,'(E16.8)', advance='no') currcnt%Ex_t(iX,iKcm)
					do ikr=currcnt%ikr_low,currcnt%ikr_high
						read(101,'(E16.8,E16.8)', advance='no') currcnt%Psi_t(ikr,iX,iKcm)
					enddo
				enddo
				
				read(100,'(E16.8)')
				read(101,'(E16.8)')
			enddo
			close(100)
			close(101)

			currcnt%ex_symmetry = 0.d0
		case ('Ex1_Em')
			call writeLog("currcnt target exciton: Ex1_Em")
			currcnt%targetExcitonType = 'Ex1_Em'
			
			allocate(currcnt%Ex_t(1:currcnt%nX_e,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
			allocate(currcnt%Psi_t(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%nX_e,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
			
			open(unit=100,file=trim(currcnt%directory)//'Ex1_Em.dat',status="old")
			open(unit=101,file=trim(currcnt%directory)//'Psi1_Em.dat',status="old")
			do iKcm=currcnt%iKcm_min_fine,currcnt%iKcm_max_fine
				do iX=1,currcnt%nX_e
					read(100,'(E16.8)', advance='no') currcnt%Ex_t(iX,iKcm)
					do ikr=currcnt%ikr_low,currcnt%ikr_high
						read(101,'(E16.8,E16.8)', advance='no') currcnt%Psi_t(ikr,iX,iKcm)
					enddo
				enddo
				
				read(100,'(E16.8)')
				read(101,'(E16.8)')
			enddo
			close(100)
			close(101)

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

		!make sure the exciton wavefunctions are normalized
		do iX=1,currcnt%nX_t
			do iKcm=currcnt%iKcm_min_fine,currcnt%iKcm_max_fine
				tmpr = 0.d0
				do ikr=currcnt%ikr_low,currcnt%ikr_high
					tmpr = tmpr + abs(currcnt%Psi_t(ikr,iX,iKcm))
				enddo
				currcnt%Psi_t(:,iX,iKcm) = currcnt%Psi_t(:,iX,iKcm) / dcmplx(sqrt(tmpr))
			enddo
		enddo

		return

	end subroutine set_target_exciton

end module target_exciton_mod