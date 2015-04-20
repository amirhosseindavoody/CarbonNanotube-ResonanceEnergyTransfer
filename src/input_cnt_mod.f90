!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Declaration of input parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
module input_cnt_mod
	implicit none
	private
	public :: input_cnt

contains	
	!**************************************************************************************************************************
	! read the CNT information from the log file of exciton calculation which is done previously
	!**************************************************************************************************************************

	subroutine input_cnt(currcnt)
		use physicalConstants, only: eV
		use cnt_class, only: cnt, cnt_geometry, cnt_band
		use write_log_mod, only: writeLog

		type(cnt) :: currcnt
		character(len=200) :: buffer, label
		integer :: ios=0
		integer :: pos=0
		integer, parameter :: nparam=10
		integer :: iparam=0
		character(len=100) :: logInput

		open(unit=100, file= trim(currcnt%directory)//'log.dat', status="old")

		iparam=0
		ios=0
		do while ((ios == 0) .and. (iparam .lt. nparam))
			read (100,'(A)',iostat=ios) buffer
			if (ios == 0) then
				pos = scan(buffer,'=')
				label = buffer(1:pos-1)
				buffer = buffer(pos+1:)
				label = adjustl(label)
				buffer = adjustl(buffer)

				select case (trim(label))
				case ('n_ch')
					read(buffer, *, iostat=ios) currcnt%n_ch
					write(logInput,"('n_ch = ', I2.2)") currcnt%n_ch
					call writeLog(logInput)
					iparam = iparam+1
				case ('m_ch')
					read(buffer, *, iostat=ios) currcnt%m_ch
					write(logInput,"('m_ch = ', I2.2)") currcnt%m_ch
					call writeLog(logInput)
					iparam = iparam+1
				case ('nkg')
					read(buffer, *, iostat=ios) currcnt%nkg
					write(logInput,"('nkg = ', I4.4)") currcnt%nkg
					call writeLog(logInput)
					iparam = iparam+1
				case ('nr')
					read(buffer, *, iostat=ios) currcnt%nr
					write(logInput,"('nr = ', I4.4)") currcnt%nr
					call writeLog(logInput)
					iparam = iparam+1
				case ('E_th')
					read(buffer, *, iostat=ios) currcnt%E_th
					write(logInput,"('E_th[eV] = ', f3.1)") currcnt%E_th
					call writeLog(logInput)
					iparam = iparam+1
					currcnt%E_th=currcnt%E_th*eV
				case ('Kcm_max')
					read(buffer, *, iostat=ios) currcnt%Kcm_max
					write(logInput,"('Kcm_max[1/nm] = ', f3.1)") currcnt%Kcm_max
					call writeLog(logInput)
					currcnt%Kcm_max = currcnt%Kcm_max*1.d9
					iparam = iparam+1
				case ('i_sub')
					read(buffer, *, iostat=ios) currcnt%i_sub
					write(logInput,"('i_sub = ', I1.1)") currcnt%i_sub
					call writeLog(logInput)
					iparam = iparam+1
				case ('Ckappa')
					read(buffer, *, iostat=ios) currcnt%Ckappa
					write(logInput,"('Ckappa = ', f3.1)") currcnt%Ckappa
					call writeLog(logInput)
					iparam = iparam+1
				case ('kappa')
					read(buffer, *, iostat=ios) currcnt%kappa
					write(logInput,"('kappa = ', f3.1)") currcnt%kappa
					call writeLog(logInput)
					iparam = iparam+1
				case ('nX')
					read(buffer, *, iostat=ios) currcnt%nX
					write(logInput,"('nX = ', I3.3)") currcnt%nX
					call writeLog(logInput)
					iparam = iparam+1
				end select
			end if
		end do

		close(100)

		if (iparam .lt. nparam) then
			call writeLog("Error in reading parameters!")
			call exit()
		end if

		! create the cnt object and calculate bands and load exciton wavefunction
		call cnt_geometry(currcnt)
		call cnt_band(currcnt)
		call input_exciton(currcnt)

		return
	end subroutine input_cnt

	!**************************************************************************************************************************
	! load the exciton wavefunction and energies from the ExcitonEnergy calculation
	!**************************************************************************************************************************
	
	subroutine input_exciton(currcnt)
		use cnt_class, only: cnt
		use write_log_mod, only: writeLog

		type(cnt), intent(inout) :: currcnt
		integer :: iX, iKcm, ikr
		real*8 :: tmpr
		
		select case (trim(currcnt%targetExcitonType))
		case ('Ex_A1')
			call writeLog("Target exciton: Ex_A1")
			open(unit=100,file=trim(currcnt%directory)//'Ex_A1.dat',status="old")
			open(unit=101,file=trim(currcnt%directory)//'Psi_A1.dat',status="old")
			currcnt%ex_symmetry = -1.d0
		case ('Ex0_A2')
			call writeLog("Target exciton: Ex0_A2")
			open(unit=100,file=trim(currcnt%directory)//'Ex0_A2.dat',status="old")
			open(unit=101,file=trim(currcnt%directory)//'Psi0_A2.dat',status="old")
			currcnt%ex_symmetry = +1.d0
		case ('Ex1_A2')
			call writeLog("Target exciton: Ex1_A2")
			open(unit=100,file=trim(currcnt%directory)//'Ex1_A2.dat',status="old")
			open(unit=101,file=trim(currcnt%directory)//'Psi1_A2.dat',status="old")
			currcnt%ex_symmetry = +1.d0
		case default
			call writeLog("Could not recognize target exciton type!!!!")
			call exit()
		end select

		allocate(currcnt%Ex_t(1:currcnt%nX,currcnt%iKcm_min:currcnt%iKcm_max))
		allocate(currcnt%Psi_t(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%nX,currcnt%iKcm_min:currcnt%iKcm_max))
		
		do iKcm=currcnt%iKcm_min,currcnt%iKcm_max
			do iX=1,currcnt%nX
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
				
		!make sure the exciton wavefunctions are normalized
		do iX=1,currcnt%nX
			do iKcm=currcnt%iKcm_min,currcnt%iKcm_max
				tmpr = 0.d0
				do ikr=currcnt%ikr_low,currcnt%ikr_high
					tmpr = tmpr + abs(currcnt%Psi_t(ikr,iX,iKcm))
				enddo
				currcnt%Psi_t(:,iX,iKcm) = currcnt%Psi_t(:,iX,iKcm) / dcmplx(sqrt(tmpr))
			enddo
		enddo

		return
	end subroutine input_exciton

	
end module input_cnt_mod