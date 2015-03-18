module data_class
    use cnt_class
    implicit none
    
    private

    public  :: saveTransitionPoints, saveDOS
    
contains
      
    !**************************************************************************************************************************
    ! save CNT dispersions and the crossing points and the same energy points
    !**************************************************************************************************************************
    subroutine saveTransitionPoints(cnt1,cnt2)
        use prepareForster_module
		use comparams
        use smallFunctions
        type(cnt), intent(in) :: cnt1,cnt2
        integer :: iKcm, iX, i
        
        !write carbon nanotube 1 k_vector
        open(unit=100,file='cnt1_kvec.dat',status="unknown")
        do iKcm=cnt1%iKcm_min,cnt1%iKcm_max
          write(100,10, advance='no') dble(iKcm)*cnt1%dk
        enddo
        close(100)
        
        !write carbon nanotube 2 k_vector
        open(unit=100,file='cnt2_kvec.dat',status="unknown")
        do iKcm=cnt2%iKcm_min,cnt2%iKcm_max
          write(100,10, advance='no') dble(iKcm)*cnt2%dk
        enddo
        close(100)
        
        !write carbon nanotube 1 Ex0_A2 dispersion
        open(unit=100,file='cnt1_Ex0_A2.dat',status="unknown")
        do iKcm=cnt1%iKcm_min,cnt1%iKcm_max
          do iX=1,cnt1%nX
            write(100,10, advance='no') cnt1%Ex0_A2(iX,iKcm)
          enddo
          write(100,10)
        enddo
        close(100)
    
        !write carbon nanotube 2 Ex0_A2 dispersion
        open(unit=100,file='cnt2_Ex0_A2.dat',status="unknown")
        do iKcm=cnt2%iKcm_min,cnt2%iKcm_max
          do iX=1,cnt2%nX
            write(100,10, advance='no') cnt2%Ex0_A2(iX,iKcm)
          enddo
          write(100,10)
        enddo
        close(100) 
        
        !write crossing points indexes
        open(unit=100,file='crossingPoints.dat',status="unknown")
        do i=1,ubound(crossingPoints,1)
          write(100,11) crossingPoints(i,1), crossingPoints(i,2), crossingPoints(i,3)
        enddo
        close(100) 
        
				!write crossing points indexes
        open(unit=100,file='sameEnergy.dat',status="unknown")
        do i=1,ubound(sameEnergy,1)
          write(100,12) sameEnergy(i,1), sameEnergy(i,2), sameEnergy(i,3), sameEnergy(i,4)
        enddo
        close(100) 
				
10		FORMAT (E16.8)
11		FORMAT (4I8, 4I8, 4I8)
12		FORMAT (4I8, 4I8, 4I8, 4I8)
					 
        return    
	end subroutine saveTransitionPoints
			
	!**************************************************************************************************************************
    ! save total exciton density of states for a given cnt
    !**************************************************************************************************************************
      subroutine saveDOS(cnt1, cnt2)
        use prepareForster_module
		use comparams
        use smallFunctions
        type(cnt), intent(in) :: cnt1, cnt2
        real*8 :: dos
        integer :: iKcm, iX
        
        !write cnt1 Ex0_A2 dispersion
        open(unit=100,file='cnt1_DOS.dat',status="unknown")
        do iKcm=cnt1%iKcm_min,cnt1%iKcm_max
          do iX=1,cnt1%nX
						call calculateDOS(cnt1,iKcm,iX,DOS)
            write(100,10, advance='no') dos
          enddo
          write(100,10)
        enddo
        close(100)
        
				!write cnt2 Ex0_A2 dispersion
        open(unit=100,file='cnt2_DOS.dat',status="unknown")
        do iKcm=cnt2%iKcm_min,cnt2%iKcm_max
          do iX=1,cnt2%nX
						call calculateDOS(cnt2,iKcm,iX,DOS)
            write(100,10, advance='no') dos
          enddo
          write(100,10)
        enddo
        close(100)
				
        10 FORMAT (E16.8)
        
        return    
	end subroutine saveDOS
end module data_class