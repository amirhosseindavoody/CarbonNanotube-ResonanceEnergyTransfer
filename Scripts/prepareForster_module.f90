module prepareForster_module
	use cntClass
	implicit none
	private
    
	integer, dimension(:,:), allocatable, public :: crossingPoints, sameEnergy
	complex*16, dimension(:), allocatable, public :: kSpaceMatrixElement
		
	public  :: findCrossings, findSameEnergy, calculateDOS, calculatePartitionFunction
    
	contains
	!**************************************************************************************************************************
	! find the points that the bands cross each other
	!**************************************************************************************************************************
	subroutine findCrossings(cnt1,cnt2)
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
					rtmp1 = (cnt1%Ex0_A2(ix1,iKcm)-cnt2%Ex0_A2(ix2,iKcm))
					rtmp2 = (cnt1%Ex0_A2(ix1,iKcm-1)-cnt2%Ex0_A2(ix2,iKcm-1))
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
        
	end subroutine findCrossings
			
	!**************************************************************************************************************************
	! find the points that the bands cross each other
	!**************************************************************************************************************************
	subroutine findSameEnergy(cnt1,cnt2)
		type(cnt), intent(in) :: cnt1,cnt2
		integer :: ix1,ix2
		integer :: iKcm1, iKcm2
		integer :: tmp,i, nSameEnergy
		integer, dimension(:,:), allocatable :: tempSameEnergy
		integer :: iC
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
					do iKcm2 = cnt2%iKcm_min+1 , cnt2%iKcm_max
						rtmp1 = (cnt1%Ex0_A2(ix1,iKcm1)-cnt2%Ex0_A2(ix2,iKcm2))
						rtmp2 = (cnt1%Ex0_A2(ix1,iKcm1)-cnt2%Ex0_A2(ix2,iKcm2-1))
						if ((rtmp1 * rtmp2) .le. 0.d0) then
							if ((iKcm2 .eq. iKcm1) .and. (rtmp2 .ne. 0.d0)) then
								nSameEnergy = nSameEnergy+1
								tempSameEnergy(nSameEnergy,1) = ix1
								tempSameEnergy(nSameEnergy,2) = ix2
								tempSameEnergy(nSameEnergy,3) = iKcm1
								tempSameEnergy(nSameEnergy,4) = iKcm2
							else if (((iKcm2-1) .eq. iKcm1) .and. (rtmp1 .ne. 0.d0)) then
								nSameEnergy = nSameEnergy+1
								tempSameEnergy(nSameEnergy,1) = ix1
								tempSameEnergy(nSameEnergy,2) = ix2
								tempSameEnergy(nSameEnergy,3) = iKcm1
								tempSameEnergy(nSameEnergy,4) = iKcm2-1
							else if ((abs(rtmp1) .le. abs(rtmp2)) .or. (abs(rtmp1) .eq. 0.d0)) then
								nSameEnergy = nSameEnergy+1
								tempSameEnergy(nSameEnergy,1) = ix1
								tempSameEnergy(nSameEnergy,2) = ix2
								tempSameEnergy(nSameEnergy,3) = iKcm1
								tempSameEnergy(nSameEnergy,4) = iKcm2
							else
								nSameEnergy = nSameEnergy+1
								tempSameEnergy(nSameEnergy,1) = ix1
								tempSameEnergy(nSameEnergy,2) = ix2
								tempSameEnergy(nSameEnergy,3) = iKcm1
								tempSameEnergy(nSameEnergy,4) = iKcm2-1
							endif
						end if
					end do
				end do
			end do
		end do
				
		!nSameEnergy = 0
		!do ix1 = 1,cnt1%nX
		!	do iKcm1 = cnt1%iKcm_min , cnt1%iKcm_max
		!		do ix2 = 1,cnt2%nX
		!			do iKcm2 = cnt2%iKcm_min+1 , cnt2%iKcm_max
		!				rtmp1 = (cnt1%Ex0_A2(ix1,iKcm1)-cnt2%Ex0_A2(ix2,iKcm2))
		!				rtmp2 = (cnt1%Ex0_A2(ix1,iKcm1)-cnt2%Ex0_A2(ix2,iKcm2-1))
		!				if ((rtmp1 * rtmp2) .le. 0.d0) then
		!					if ((abs(iKcm2-iKcm1) .le. abs(iKcm2-1-iKcm1)) .or. (abs(iKcm2-iKcm1) .eq. 0)) then
		!						nSameEnergy = nSameEnergy+1
		!						tempSameEnergy(nSameEnergy,1) = ix1
		!						tempSameEnergy(nSameEnergy,2) = ix2
		!						tempSameEnergy(nSameEnergy,3) = iKcm1
		!						tempSameEnergy(nSameEnergy,4) = iKcm2
		!					else
		!						nSameEnergy = nSameEnergy+1
		!						tempSameEnergy(nSameEnergy,1) = ix1
		!						tempSameEnergy(nSameEnergy,2) = ix2
		!						tempSameEnergy(nSameEnergy,3) = iKcm1
		!						tempSameEnergy(nSameEnergy,4) = iKcm2-1
		!					endif
		!				end if
		!			end do
		!	  end do
		!	end do
		!end do
				
        
		if(allocated(sameEnergy))	deallocate(sameEnergy)
        allocate(sameEnergy(nSameEnergy ,4))
		sameEnergy(:,:) = tempSameEnergy(1:nSameEnergy,:)
				
		if(allocated(kSpaceMatrixElement))	deallocate(kSpaceMatrixElement)
		allocate(kSpaceMatrixElement(nSameEnergy))
		kSpaceMatrixElement = kSpaceMatrixElement * 0.d0
				
		do iC = 1,nSameEnergy
			ix1 = sameEnergy(iC,1)
			ix2 = sameEnergy(iC,2)
			iKcm1 = sameEnergy(iC,3)
			iKcm2 = sameEnergy(iC,4)
			if ((iKcm1 .ne. 0) .and. (iKcm2 .ne. 0)) then
				call calculateKSpaceMatrixElement(cnt1,cnt2,ix1, ix2, iKcm1, iKcm2, kSpaceMatrixElement(iC))
			end if
		end do
        
	end subroutine findSameEnergy
      
      
	!**************************************************************************************************************************
	! calculate the partition function for a given carbon nanotube
	!**************************************************************************************************************************
	subroutine calculatePartitionFunction(currcnt, partitionFunction)
		use physicalConstants, only : kb
		use inputParameters, only : Temperature
		type(cnt), intent(in) :: currcnt
		real*8, intent(out) :: partitionFunction
		integer :: ix, iKcm
        
		partitionFunction = 0.d0
        
		do ix = 1,currcnt%nX
			do iKcm = currcnt%iKcm_min,currcnt%iKcm_max
				partitionFunction = partitionFunction + currcnt%dk * exp(-currcnt%Ex0_A2(ix,iKcm)/kb/Temperature)    
			end do
		end do
        
	end subroutine calculatePartitionFunction

	!**************************************************************************************************************************
	! calculate the density of states at a given point
	!**************************************************************************************************************************
	subroutine calculateDOS(currcnt,iKcm,iX,DOS)
        type(cnt), intent(in) :: currcnt
        integer, intent(in) :: iKcm,iX
        real*8, intent(out) :: DOS
        
		if (iKcm == currcnt%iKcm_min) then
			DOS = currcnt%dk/abs(currcnt%Ex0_A2(iX,iKcm)-currcnt%Ex0_A2(iX,iKcm+1))
		else if(iKcm == currcnt%iKcm_max) then
			DOS = currcnt%dk/abs(currcnt%Ex0_A2(iX,iKcm-1)-currcnt%Ex0_A2(iX,iKcm))
		else if(iKcm == 0) then
			DOS = currcnt%dk/abs(currcnt%Ex0_A2(iX,iKcm-1)-currcnt%Ex0_A2(iX,iKcm))
		else
			DOS = 2.d0*currcnt%dk/abs(currcnt%Ex0_A2(iX,iKcm-1)-currcnt%Ex0_A2(iX,iKcm+1))
		end if
	end subroutine calculateDOS
			
	!**************************************************************************************************************************
	! calculate the matrix element for the crossing point number iC in two unparallel tube
	!**************************************************************************************************************************
	subroutine calculateKSpaceMatrixElement(cnt1,cnt2,ix1, ix2, iKcm1, iKcm2 ,kSpaceMatrixElementTemp)
		use inputParameters
		use physicalConstants, only : i1, pi, eps0, q0
		type(cnt), intent(in) :: cnt1,cnt2
		complex*16, intent(out) :: kSpaceMatrixElementTemp
		complex*16 :: tmpc
		integer, intent(in) :: ix1,ix2
		integer, intent(in) :: iKcm1, iKcm2
		integer :: ikr1, ikr2
		integer :: is,isp
		
		kSpaceMatrixElementTemp = (0.d0,0.d0)
		tmpc = (0.d0,0.d0)
		
		do ikr1 = cnt1%ikr_low, cnt1%ikr_high
			do ikr2 = cnt2%ikr_low, cnt2%ikr_high					
				do is = 1,2
					do isp = 1,2
						tmpc = tmpc + conjg(cnt1%Cc(1,ikr1+iKcm1,is))*cnt1%Cv(1,ikr1-iKcm1,is)*cnt2%Cc(1,ikr2+iKcm2,isp)*conjg(cnt2%Cv(1,ikr2-iKcm2,isp))
						tmpc = tmpc + conjg(cnt1%Cc(1,ikr1+iKcm1,is))*cnt1%Cv(1,ikr1-iKcm1,is)*cnt2%Cc(2,-ikr2+iKcm2,isp)*conjg(cnt2%Cv(2,-ikr2-iKcm2,isp))
						tmpc = tmpc + conjg(cnt1%Cc(2,-ikr1+iKcm1,is))*cnt1%Cv(2,-ikr1-iKcm1,is)*cnt2%Cc(1,ikr2+iKcm2,isp)*conjg(cnt2%Cv(1,ikr2-iKcm2,isp))
						tmpc = tmpc + conjg(cnt1%Cc(2,-ikr1+iKcm1,is))*cnt1%Cv(2,-ikr1-iKcm1,is)*cnt2%Cc(2,-ikr2+iKcm2,isp)*conjg(cnt2%Cv(2,-ikr2-iKcm2,isp))
					end do  
				end do
				kSpaceMatrixElementTemp = kSpaceMatrixElementTemp + tmpc*conjg(cnt1%Psi0_A2(ikr1,ix1,iKcm1))*cnt2%Psi0_A2(ikr2,ix2,iKcm2)/(2.d0,0.d0)
				tmpc = (0.d0,0.d0)
			end do
		end do				
	end subroutine calculateKSpaceMatrixElement
			
end module prepareForster_module