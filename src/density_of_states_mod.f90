module density_of_states_mod
	implicit none
	private
    public  :: calculate_dos, save_dos
    
contains

	!**************************************************************************************************************************
	! calculate the density of states at a given point
	!**************************************************************************************************************************
	
	subroutine calculate_dos(currcnt,iKcm,iX,dos)
		use cnt_class, only: cnt

		type(cnt), intent(in) :: currcnt
		integer, intent(in) :: iKcm,iX
		real*8, intent(out) :: dos
		
		if (iKcm .le. (currcnt%iKcm_min_fine+1)) then
			dos = 2.d0*currcnt%dkx/abs(-3.d0*currcnt%Ex_t(iX,iKcm)+4.d0*currcnt%Ex_t(iX,iKcm+1)-currcnt%Ex_t(iX,iKcm+2))
		else if(iKcm .ge. (currcnt%iKcm_max_fine-1)) then
			dos = 2.d0*currcnt%dkx/abs(3.d0*currcnt%Ex_t(iX,iKcm)-4.d0*currcnt%Ex_t(iX,iKcm-1)+currcnt%Ex_t(iX,iKcm-2))
		else if(iKcm == 0) then
			dos = 2.d0*currcnt%dkx/abs(3.d0*currcnt%Ex_t(iX,iKcm)-4.d0*currcnt%Ex_t(iX,iKcm-1)+currcnt%Ex_t(iX,iKcm-2))
		else
			dos = 12.d0*currcnt%dkx / abs(currcnt%Ex_t(iX,iKcm-2)-8.d0*currcnt%Ex_t(iX,iKcm-1)+8.d0*currcnt%Ex_t(iX,iKcm+1)-currcnt%Ex_t(iX,iKcm+2))
		end if
		return
	end subroutine calculate_dos

	!**************************************************************************************************************************
	! save total exciton density of states for a given cnt
	!**************************************************************************************************************************
	
	subroutine save_dos(currcnt)
		use cnt_class, only: cnt

		type(cnt), intent(in) :: currcnt
		real*8 :: dos
		integer :: iKcm, iX
		character(len=100) :: filename

		write(filename,"( 'dos-cnt(', I2.2, ',', I2.2, ')-i_sub(', I1.1, ')-Ckappa(', F3.1, ').dat' )") currcnt%n_ch, currcnt%m_ch, currcnt%i_sub, currcnt%Ckappa

		!write currcnt dos
		open(unit=100,file=filename,status="unknown")
		do iKcm=currcnt%iKcm_min_fine,currcnt%iKcm_max_fine
			do iX=1,currcnt%nX_t
				call calculate_dos(currcnt,iKcm,iX,DOS)
				write(100,'(E16.8)', advance='no') dos
			enddo
			write(100,'(E16.8)')
		enddo
		close(100)

		return    
	end subroutine save_dos
			
end module density_of_states_mod