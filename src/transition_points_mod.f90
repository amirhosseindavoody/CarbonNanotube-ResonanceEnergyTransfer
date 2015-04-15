module transition_points_mod
	implicit none
	private
    public  :: findCrossings, findSameEnergy

	integer, dimension(:,:), allocatable, public :: crossingPoints, sameEnergy
    
contains
	!**************************************************************************************************************************
	! find the points that the bands cross each other
	!**************************************************************************************************************************
	
	subroutine findCrossings(cnt1,cnt2)
		use cnt_class, only: cnt

		type(cnt), intent(in) :: cnt1,cnt2
		integer :: ix1,ix2
		integer :: iKcm
		integer :: tmp,i, nCrossing
		integer, dimension(:,:), allocatable :: tempCrossingPoints
		real*8 :: rtmp1, rtmp2
        
		tmp = cnt1%nX * cnt2%nX * (2*cnt1%iKcm_max+1) * (2*cnt2%iKcm_max+1)
		allocate(tempCrossingPoints(tmp,3))
		do i = 1,tmp
			tempCrossingPoints(i,1) = 0    
			tempCrossingPoints(i,2) = 0
			tempCrossingPoints(i,3) = 0
		end do
        
		nCrossing = 0
		do ix1 = 1,cnt1%nX
			do ix2 = 1,cnt2%nX
				do iKcm = cnt1%iKcm_min+1 , cnt1%iKcm_max
					rtmp1 = (cnt1%Ex_t(ix1,iKcm)-cnt2%Ex_t(ix2,iKcm))
					rtmp2 = (cnt1%Ex_t(ix1,iKcm-1)-cnt2%Ex_t(ix2,iKcm-1))
					if ((rtmp1 * rtmp2) .le. 0.d0) then
						if ((abs(rtmp1) .le. abs(rtmp2)) .or. (abs(rtmp1) .eq. 0.d0)) then
							nCrossing = nCrossing+1
							tempCrossingPoints(nCrossing,1) = ix1
							tempCrossingPoints(nCrossing,2) = ix2
							tempCrossingPoints(nCrossing,3) = iKcm
						else
							nCrossing = nCrossing+1
							tempCrossingPoints(nCrossing,1) = ix1
							tempCrossingPoints(nCrossing,2) = ix2
							tempCrossingPoints(nCrossing,3) = iKcm-1
						endif
					end if
				end do
			end do
		end do
        
		if(allocated(crossingPoints))	deallocate(crossingPoints)
		allocate(crossingPoints(nCrossing ,3))
		crossingPoints(:,:) = tempCrossingPoints(1:nCrossing,:)
        return
	end subroutine findCrossings
			
	!**************************************************************************************************************************
	! find the points that the bands cross each other
	!**************************************************************************************************************************
	
	subroutine findSameEnergy(cnt1,cnt2)
		use cnt_class, only: cnt
		use write_log_mod, only: writeLog

		type(cnt), intent(in) :: cnt1,cnt2
		integer :: ix1,ix2
		integer :: iKcm1, iKcm2
		integer :: tmp,i, nSameEnergy
		integer, dimension(:,:), allocatable :: tempSameEnergy
		real*8 :: rtmp1, rtmp2

		tmp = cnt1%nX * cnt2%nX * (2*cnt1%iKcm_max+1) * (2*cnt2%iKcm_max+1)
		allocate(tempSameEnergy(tmp,4))
		do i = 1,tmp
			tempSameEnergy(i,1) = 0    
			tempSameEnergy(i,2) = 0
			tempSameEnergy(i,3) = 0
			tempSameEnergy(i,4) = 0
		end do
        
		nSameEnergy = 0
		do ix1 = 1,cnt1%nX
			do iKcm1 = cnt1%iKcm_min , cnt1%iKcm_max
				do ix2 = 1,cnt2%nX
					
					if (cnt1%Ex_t(ix1,iKcm1) .eq. cnt2%Ex_t(ix2,cnt2%iKcm_min)) then
							nSameEnergy = nSameEnergy+1
							tempSameEnergy(nSameEnergy,1) = ix1
							tempSameEnergy(nSameEnergy,2) = ix2
							tempSameEnergy(nSameEnergy,3) = iKcm1
							tempSameEnergy(nSameEnergy,4) = cnt2%iKcm_min
					endif

					do iKcm2 = cnt2%iKcm_min+1 , cnt2%iKcm_max
						rtmp1 = (cnt1%Ex_t(ix1,iKcm1)-cnt2%Ex_t(ix2,iKcm2))
						rtmp2 = (cnt1%Ex_t(ix1,iKcm1)-cnt2%Ex_t(ix2,iKcm2-1))
						if (rtmp1 .eq. 0.d0) then
							nSameEnergy = nSameEnergy+1
							tempSameEnergy(nSameEnergy,1) = ix1
							tempSameEnergy(nSameEnergy,2) = ix2
							tempSameEnergy(nSameEnergy,3) = iKcm1
							tempSameEnergy(nSameEnergy,4) = iKcm2
						elseif ((rtmp1 * rtmp2) .lt. 0.d0) then
							nSameEnergy = nSameEnergy+1
							if (abs(rtmp1) .lt. abs(rtmp2)) then
								tempSameEnergy(nSameEnergy,1) = ix1
								tempSameEnergy(nSameEnergy,2) = ix2
								tempSameEnergy(nSameEnergy,3) = iKcm1
								tempSameEnergy(nSameEnergy,4) = iKcm2
							else
								tempSameEnergy(nSameEnergy,1) = ix1
								tempSameEnergy(nSameEnergy,2) = ix2
								tempSameEnergy(nSameEnergy,3) = iKcm1
								tempSameEnergy(nSameEnergy,4) = iKcm2-1
							endif
						endif
					end do
				end do
			end do
		end do	
        
		if(allocated(sameEnergy))	deallocate(sameEnergy)
        allocate(sameEnergy(nSameEnergy ,4))
		sameEnergy(:,:) = tempSameEnergy(1:nSameEnergy,:)

		call writeLog(new_line('A')//"Same energy table calculated!!!"//new_line('A'))

		call saveTransitionPoints(cnt1,cnt2)

        return
	end subroutine findSameEnergy

	!**************************************************************************************************************************
	! save CNT dispersions and the crossing points and the same energy points
	!**************************************************************************************************************************
	
	subroutine saveTransitionPoints(cnt1,cnt2)
		use cnt_class, only: cnt

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

		!write carbon nanotube 1 Ex_t dispersion
		open(unit=100,file='cnt1_Ex_t.dat',status="unknown")
		do iKcm=cnt1%iKcm_min,cnt1%iKcm_max
			do iX=1,cnt1%nX
				write(100,10, advance='no') cnt1%Ex_t(iX,iKcm)
			enddo
			write(100,10)
		enddo
		close(100)

		!write carbon nanotube 2 Ex_t dispersion
		open(unit=100,file='cnt2_Ex_t.dat',status="unknown")
		do iKcm=cnt2%iKcm_min,cnt2%iKcm_max
			do iX=1,cnt2%nX
				write(100,10, advance='no') cnt2%Ex_t(iX,iKcm)
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
			
end module transition_points_mod