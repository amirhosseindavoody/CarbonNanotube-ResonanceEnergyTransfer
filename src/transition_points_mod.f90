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
		use write_log_mod, only: writeLog

		type(cnt), intent(in) :: cnt1,cnt2
		integer :: ix1
		integer :: ikr2
		integer :: iKcm
		integer :: nCrossing
		real*8 :: rtmp1, rtmp2
		character(len=200) :: logInput

		nCrossing = 0
		do ix1 = 1,cnt1%nX_t
			do ikr2 = cnt2%ikr_low,cnt2%ikr_high
				do iKcm = cnt1%iKcm_min_fine , -1
					rtmp1 = (cnt1%Ex_t(ix1,iKcm)-cnt2%E_free_eh(1,1,ikr2,iKcm))
					rtmp2 = (cnt1%Ex_t(ix1,iKcm+1)-cnt2%E_free_eh(1,1,ikr2,iKcm+1))
					if ((rtmp1 * rtmp2) .le. 0.d0) then
						nCrossing = nCrossing+2
					end if
				end do
			end do
		end do

		write(logInput,*) "Number of crossing points = ", nCrossing
		call writeLog(logInput)

		if(allocated(crossingPoints))	deallocate(crossingPoints)
		allocate(crossingPoints(nCrossing ,4))

		nCrossing = 0
		do ix1 = 1,cnt1%nX_t
			do ikr2 = cnt2%ikr_low,cnt2%ikr_high
				do iKcm = cnt1%iKcm_min_fine , -1
					rtmp1 = (cnt1%Ex_t(ix1,iKcm)-cnt2%E_free_eh(1,1,ikr2,iKcm))
					rtmp2 = (cnt1%Ex_t(ix1,iKcm+1)-cnt2%E_free_eh(1,1,ikr2,iKcm+1))
					if ((rtmp1 * rtmp2) .le. 0.d0) then
						nCrossing = nCrossing+1
						crossingPoints(nCrossing,1) = ix1
						crossingPoints(nCrossing,2) = ikr2
						crossingPoints(nCrossing,3) = iKcm
						crossingPoints(nCrossing,4) = iKcm

						nCrossing = nCrossing+1
						crossingPoints(nCrossing,1) = ix1
						crossingPoints(nCrossing,2) = ikr2
						crossingPoints(nCrossing,3) = -iKcm
						crossingPoints(nCrossing,4) = -iKcm
					end if
				end do
			end do
		end do

		call writeLog(new_line('A')//"Crossing points table calculated!!!"//new_line('A'))

		call saveTransitionPoints(cnt1,cnt2)

        return
	end subroutine findCrossings
			
	!**************************************************************************************************************************
	! find the points that the bands have equal energy
	!**************************************************************************************************************************
	
	subroutine findSameEnergy(cnt1,cnt2)
		use cnt_class, only: cnt
		use comparams, only: Temperature
		use math_functions_mod, only: bisect_root
		use physicalConstants, only: kb
		use write_log_mod, only: writeLog
		character(len=200) :: logInput

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
		deltaE = (-1.d0) * log(1.d-3) * kb*Temperature

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
		write(logInput,*) "Number of same energy points = ", nSameEnergy
		call writeLog(logInput)

! 		call saveTransitionPoints(cnt1,cnt2)

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

		!write crossing points indexes
		open(unit=100,file='crossingPoints.dat',status="unknown")
		do i=lbound(crossingPoints,1),ubound(crossingPoints,1)
			write(100,'(4I8, 4I8, 4I8, 4I8)') crossingPoints(i,1), crossingPoints(i,2)-cnt2%ikr_low+1, crossingPoints(i,3)-cnt1%iKcm_min_fine+1, crossingPoints(i,4)-cnt2%iKcm_min_fine+1
		enddo
		close(100) 

! 		!write same energy points indexes for transition from cnt1 to cnt2
! 		open(unit=100,file='sameEnergy.dat',status="unknown")
! 		do i=lbound(sameEnergy,1),ubound(sameEnergy,1)
! 			write(100,'(4I8, 4I8, 4I8, 4I8)') sameEnergy(i,1), sameEnergy(i,2), sameEnergy(i,3), sameEnergy(i,4)
! 		enddo
! 		close(100) 

		return    
	end subroutine saveTransitionPoints
			
end module transition_points_mod
