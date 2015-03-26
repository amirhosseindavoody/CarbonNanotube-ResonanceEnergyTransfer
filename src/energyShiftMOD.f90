!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate transition table
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module energyShiftMOD
	implicit none	
	private
	public :: shiftEnergy
	
contains
	!**************************************************************************************************************************
	! calculate transition table
	!**************************************************************************************************************************
	subroutine shiftEnergy (currcnt)
		use output_module, only: writeLog, logInput
		use physicalConstants, only: eV
		use cnt_class, only: cnt
		
		type(cnt), intent(inout) :: currcnt
		real*8 :: tmpr
		
		if ((currcnt%n_ch .eq. 7) .and. (currcnt%m_ch .eq. 5) .and. (currcnt%i_sub .eq. 2)) then
			tmpr = minval(currcnt%Ex0_A2)
			write(logInput,*) "min Ex0_A2 = ", tmpr/eV
			call writeLog()
			currcnt%Ex0_A2 = currcnt%Ex0_A2 - tmpr + 1.921 * eV
		end if
		
		if ((currcnt%n_ch .eq. 7) .and. (currcnt%m_ch .eq. 6) .and. (currcnt%i_sub .eq. 2)) then
			tmpr = minval(currcnt%Ex0_A2)
			write (logInput,*)"min Ex0_A2 = ", tmpr/eV
			call writeLog()
			currcnt%Ex0_A2 = currcnt%Ex0_A2 - tmpr + 1.914 * eV
		end if
		
		if ((currcnt%n_ch .eq. 8) .and. (currcnt%m_ch .eq. 6) .and. (currcnt%i_sub .eq. 2)) then
			tmpr = minval(currcnt%Ex0_A2)
			write(logInput,*) "min Ex0_A2 = ", tmpr/eV
			call writeLog()
			currcnt%Ex0_A2 = currcnt%Ex0_A2 - tmpr + 1.727 * eV
		end if
		
		if ((currcnt%n_ch .eq. 8) .and. (currcnt%m_ch .eq. 7) .and. (currcnt%i_sub .eq. 2)) then
			tmpr = minval(currcnt%Ex0_A2)
			write(logInput,*) "min Ex0_A2 = ", tmpr/eV
			call writeLog()
			currcnt%Ex0_A2 = currcnt%Ex0_A2 - tmpr + 1.702 * eV
		end if
		return
	end subroutine shiftEnergy
	
end module energyShiftMOD