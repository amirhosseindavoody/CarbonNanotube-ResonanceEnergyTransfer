module forsterClass
    use cntClass
    implicit none
    integer, dimension(:,:), allocatable :: finalCrossingPoints
    private
    
    public  :: findCrossings, finalCrossingPoints
    
    contains
      !**************************************************************************************************************************
      ! find the points that the bands cross each other
      !**************************************************************************************************************************
      subroutine findCrossings(cnt1,cnt2)
        use physicalConstants, only : eV
        type(cnt), intent(in) :: cnt1,cnt2
        integer :: ix1,ix2
        integer :: iKcm
        integer :: tmp,i, nCrossing
        integer, dimension(:,:), allocatable :: crossingPoints
        real*8 :: rtmp1, rtmp2
        
        tmp = cnt1.nX * cnt2.nX * (2*cnt1.iKcm_max+1) * (2*cnt2.iKcm_max+1)
        allocate(crossingPoints(tmp,3))
        do i = 1,tmp
          crossingPoints(i,1) = 0    
          crossingPoints(i,2) = 0
          crossingPoints(i,3) = 0
        end do
        
        nCrossing = 0
        do ix1 = 1,cnt1.nX
          do ix2 = 1,cnt2.nX
            do iKcm = cnt1.iKcm_min+1 , cnt1.iKcm_max
              rtmp1 = (cnt1.Ex0_A2(ix1,iKcm)-cnt2.Ex0_A2(ix2,iKcm))
              rtmp2 = (cnt1.Ex0_A2(ix1,iKcm-1)-cnt2.Ex0_A2(ix2,iKcm-1))
              if ((rtmp1 * rtmp2) .le. 0.d0) then
                if ((abs(rtmp1) .le. abs(rtmp2)) .or. (abs(rtmp1) .eq. 0.d0)) then
                  nCrossing = nCrossing+1
                  crossingPoints(nCrossing,1) = ix1
                  crossingPoints(nCrossing,2) = ix2
                  crossingPoints(nCrossing,3) = iKcm
                else
                  nCrossing = nCrossing+1
                  crossingPoints(nCrossing,1) = ix1
                  crossingPoints(nCrossing,2) = ix2
                  crossingPoints(nCrossing,3) = iKcm-1
                endif
              end if
            end do
          end do
        end do
        
        allocate(finalCrossingPoints(nCrossing ,3))
        finalCrossingPoints(:,:) = crossingPoints(1:nCrossing,:)
        print *, "nCrossing = ", nCrossing
        
      end subroutine findCrossings
      
      
    
end module forsterClass