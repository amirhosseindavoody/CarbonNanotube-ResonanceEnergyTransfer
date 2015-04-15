module prepareForster_module
	implicit none
	private
    public  :: calculateDOS, calculatePartitionFunction, saveDOS
    
contains
	!**************************************************************************************************************************
	! calculate the partition function for a given carbon nanotube
	!**************************************************************************************************************************
	
	subroutine calculatePartitionFunction(currcnt, partitionFunction)
		use physicalConstants, only : kb
		use comparams, only : Temperature
		use cnt_class, only: cnt

		type(cnt), intent(in) :: currcnt
		real*8, intent(out) :: partitionFunction
		integer :: ix, iKcm
        
		partitionFunction = 0.d0
        
		do ix = 1,currcnt%nX
			do iKcm = currcnt%iKcm_min,currcnt%iKcm_max
				partitionFunction = partitionFunction + currcnt%dk * exp(-currcnt%Ex_t(ix,iKcm)/kb/Temperature)    
			end do
		end do
        
	end subroutine calculatePartitionFunction

	!**************************************************************************************************************************
	! calculate the density of states at a given point
	!**************************************************************************************************************************
	
	subroutine calculateDOS(currcnt,iKcm,iX,dos)
		use cnt_class, only: cnt

		type(cnt), intent(in) :: currcnt
		integer, intent(in) :: iKcm,iX
		real*8, intent(out) :: dos
		
		if (iKcm .le. (currcnt%iKcm_min+1)) then
			dos = 2.d0*currcnt%dk/abs(-3.d0*currcnt%Ex_t(iX,iKcm)+4.d0*currcnt%Ex_t(iX,iKcm+1)-currcnt%Ex_t(iX,iKcm+2))
		else if(iKcm .ge. (currcnt%iKcm_max-1)) then
			dos = 2.d0*currcnt%dk/abs(3.d0*currcnt%Ex_t(iX,iKcm)-4.d0*currcnt%Ex_t(iX,iKcm-1)+currcnt%Ex_t(iX,iKcm-2))
		else if(iKcm == 0) then
			dos = 2.d0*currcnt%dk/abs(3.d0*currcnt%Ex_t(iX,iKcm)-4.d0*currcnt%Ex_t(iX,iKcm-1)+currcnt%Ex_t(iX,iKcm-2))
		else
			dos = 12.d0*currcnt%dk / abs(currcnt%Ex_t(iX,iKcm-2)-8.d0*currcnt%Ex_t(iX,iKcm-1)+8.d0*currcnt%Ex_t(iX,iKcm+1)-currcnt%Ex_t(iX,iKcm+2))
		end if
		return
	end subroutine calculateDOS

	!**************************************************************************************************************************
	! save total exciton density of states for a given cnt
	!**************************************************************************************************************************
	
	subroutine saveDOS(currcnt)
		use cnt_class, only: cnt

		type(cnt), intent(in) :: currcnt
		real*8 :: dos
		integer :: iKcm, iX
		character(len=100) :: filename

		write(filename,"( 'dos-cnt(', I2.2, ',', I2.2, ')-i_sub(', I1.1, ')-Ckappa(', F3.1, ').dat' )") currcnt%n_ch, currcnt%m_ch, currcnt%i_sub, currcnt%Ckappa

		!write currcnt dos
		open(unit=100,file=filename,status="unknown")
		do iKcm=currcnt%iKcm_min,currcnt%iKcm_max
			do iX=1,currcnt%nX
				call calculateDOS(currcnt,iKcm,iX,DOS)
				write(100,'(E16.8)', advance='no') dos
			enddo
			write(100,'(E16.8)')
		enddo
		close(100)

		return    
	end subroutine saveDOS
			
end module prepareForster_module