!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate transition table
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
module energyCorrection_module
	use cntClass
  implicit none
	
	private
	
	public :: correctEnergy
    
	contains
		!**************************************************************************************************************************
		! calculate transition table
    !**************************************************************************************************************************
		subroutine correctEnergy (currcnt)
			
			use physicalConstants	
			
			type(cnt), intent(inout) :: currcnt
			real*8 :: tmpr
			
			if ((currcnt.n_ch .eq. 7) .and. (currcnt.m_ch .eq. 5) .and. (currcnt.i_sub .eq. 2)) then
				tmpr = minval(currcnt.Ex0_A2)
				print *, "min Ex0_A2 = ", tmpr/eV
				currcnt.Ex0_A2 = currcnt.Ex0_A2 - tmpr + 1.921 * eV
			end if
			
			if ((currcnt.n_ch .eq. 7) .and. (currcnt.m_ch .eq. 6) .and. (currcnt.i_sub .eq. 2)) then
				tmpr = minval(currcnt.Ex0_A2)
				print *, "min Ex0_A2 = ", tmpr/eV
				currcnt.Ex0_A2 = currcnt.Ex0_A2 - tmpr + 1.914 * eV
			end if
			
			if ((currcnt.n_ch .eq. 8) .and. (currcnt.m_ch .eq. 6) .and. (currcnt.i_sub .eq. 2)) then
				tmpr = minval(currcnt.Ex0_A2)
				print *, "min Ex0_A2 = ", tmpr/eV
				currcnt.Ex0_A2 = currcnt.Ex0_A2 - tmpr + 1.727 * eV
			end if
			
			if ((currcnt.n_ch .eq. 8) .and. (currcnt.m_ch .eq. 7) .and. (currcnt.i_sub .eq. 2)) then
				tmpr = minval(currcnt.Ex0_A2)
				print *, "min Ex0_A2 = ", tmpr/eV
				currcnt.Ex0_A2 = currcnt.Ex0_A2 - tmpr + 1.702 * eV
			end if
        
		end subroutine correctEnergy
	
end module energyCorrection_module