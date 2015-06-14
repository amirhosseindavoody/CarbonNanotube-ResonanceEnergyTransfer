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
        
		tmp = cnt1%nX_t * cnt2%nX_t * (2*cnt1%iKcm_max_fine+1) * (2*cnt2%iKcm_max_fine+1)
		allocate(tempCrossingPoints(tmp,3))
		do i = 1,tmp
			tempCrossingPoints(i,1) = 0    
			tempCrossingPoints(i,2) = 0
			tempCrossingPoints(i,3) = 0
		end do
        
		nCrossing = 0
		do ix1 = 1,cnt1%nX_t
			do ix2 = 1,cnt2%nX_t
				do iKcm = cnt1%iKcm_min_fine+1 , cnt1%iKcm_max_fine
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
		use comparams, only: Temperature
		use math_functions_mod, only: bisect_root
		use physicalConstants, only: kb
		use write_log_mod, only: writeLog

		type(cnt), intent(in) :: cnt1,cnt2
		integer :: ix1,ix2
		integer :: iKcm1, iKcm2
		integer :: tmp, nSameEnergy
		integer :: n, iKcm_raw
		integer, dimension(:,:), allocatable :: tempSameEnergy
		real*8 :: min_energy, deltaE

		tmp = cnt1%nX_t * cnt2%nX_t * (2*cnt1%iKcm_max_fine+1) * (2*cnt2%iKcm_max_fine+1)
		allocate(tempSameEnergy(tmp,4))

		! calculate relevant same energy points for transition from cnt1 to cnt2
		tempSameEnergy = 0 * tempSameEnergy
		deltaE = (-1.d0) * log(1.d-5) * kb*Temperature

		min_energy = max(minval(cnt1%Ex_t),minval(cnt2%Ex_t))

		n = cnt2%iKcm_max_fine

		nSameEnergy = 0
		do ix1 = 1,cnt1%nX_t
			do iKcm1 = cnt1%iKcm_min_fine , -1
				if ((cnt1%Ex_t(ix1,iKcm1) - min_energy) .lt. deltaE) then
					do ix2 = 1, cnt2%nX_t
						call bisect_root(n, cnt2%Ex_t(ix2,cnt2%iKcm_min_fine:-1), cnt1%Ex_t(ix1,iKcm1), iKcm_raw)
						if (iKcm_raw .gt. 0) then
							iKcm2 = cnt2%iKcm_min_fine + iKcm_raw - 1

							nSameEnergy = nSameEnergy + 1
							tempSameEnergy(nSameEnergy, 1) = ix1
							tempSameEnergy(nSameEnergy, 2) = ix2
							tempSameEnergy(nSameEnergy, 3) = +iKcm1
							tempSameEnergy(nSameEnergy, 4) = +iKcm2

							nSameEnergy = nSameEnergy + 1
							tempSameEnergy(nSameEnergy, 1) = ix1
							tempSameEnergy(nSameEnergy, 2) = ix2
							tempSameEnergy(nSameEnergy, 3) = +iKcm1
							tempSameEnergy(nSameEnergy, 4) = -iKcm2

							nSameEnergy = nSameEnergy + 1
							tempSameEnergy(nSameEnergy, 1) = ix1
							tempSameEnergy(nSameEnergy, 2) = ix2
							tempSameEnergy(nSameEnergy, 3) = -iKcm1
							tempSameEnergy(nSameEnergy, 4) = +iKcm2

							nSameEnergy = nSameEnergy + 1
							tempSameEnergy(nSameEnergy, 1) = ix1
							tempSameEnergy(nSameEnergy, 2) = ix2
							tempSameEnergy(nSameEnergy, 3) = -iKcm1
							tempSameEnergy(nSameEnergy, 4) = -iKcm2

						endif
					enddo
				endif
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
		do iKcm=cnt1%iKcm_min_fine,cnt1%iKcm_max_fine
			write(100,'(E16.8)', advance='no') dble(iKcm)*cnt1%dkx
		enddo
		close(100)

		!write carbon nanotube 2 k_vector
		open(unit=100,file='cnt2_kvec.dat',status="unknown")
		do iKcm=cnt2%iKcm_min_fine,cnt2%iKcm_max_fine
			write(100,'(E16.8)', advance='no') dble(iKcm)*cnt2%dkx
		enddo
		close(100)

		!write carbon nanotube 1 Ex_t dispersion
		open(unit=100,file='cnt1_Ex_t.dat',status="unknown")
		do iKcm=cnt1%iKcm_min_fine,cnt1%iKcm_max_fine
			do iX=1,cnt1%nX_t
				write(100,'(E16.8)', advance='no') cnt1%Ex_t(iX,iKcm)
			enddo
			write(100,*)
		enddo
		close(100)

		!write carbon nanotube 2 Ex_t dispersion
		open(unit=100,file='cnt2_Ex_t.dat',status="unknown")
		do iKcm=cnt2%iKcm_min_fine,cnt2%iKcm_max_fine
			do iX=1,cnt2%nX_t
				write(100,'(E16.8)', advance='no') cnt2%Ex_t(iX,iKcm)
			enddo
			write(100,*)
		enddo
		close(100) 

		!write crossing points indexes
		open(unit=100,file='crossingPoints.dat',status="unknown")
		do i=1,ubound(crossingPoints,1)
			write(100,'(4I8, 4I8, 4I8)') crossingPoints(i,1), crossingPoints(i,2), crossingPoints(i,3)
		enddo
		close(100) 

		!write same energy points indexes for transition from cnt1 to cnt2
		open(unit=100,file='sameEnergy.dat',status="unknown")
		do i=1,ubound(sameEnergy,1)
			write(100,'(4I8, 4I8, 4I8, 4I8)') sameEnergy(i,1), sameEnergy(i,2), sameEnergy(i,3), sameEnergy(i,4)
		enddo
		close(100) 

		return    
	end subroutine saveTransitionPoints
			
end module transition_points_mod
