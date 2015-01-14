module forsterClass
    use cntClass
    implicit none
    private
    
    integer, dimension(:,:), allocatable, public :: finalCrossingPoints
    real*8, public :: partitionFunction
    
    public  :: findCrossings, calculatePartitionFunction
    
    contains
      !**************************************************************************************************************************
      ! find the points that the bands cross each other
      !**************************************************************************************************************************
      subroutine findCrossings(cnt1,cnt2)
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
      
      
      !**************************************************************************************************************************
      ! calculate the partition function for a given carbon nanotube
      !**************************************************************************************************************************
      subroutine calculatePartitionFunction(currcnt)
        use physicalConstants, only : kb
        use inputParameters, only : Temperature
        type(cnt), intent(in) :: currcnt
        integer :: ix, iKcm
        
        partitionFunction = 0.d0
        
        do ix = 1,currcnt.nX
          do iKcm = currcnt.iKcm_min,currcnt.iKcm_max
            partitionFunction = partitionFunction + currcnt.dk * exp(-currcnt.Ex0_A2(ix,iKcm)/kb/Temperature)    
          end do
        end do
        
        print *, "partitionFunction = ", partitionFunction
        
      end subroutine calculatePartitionFunction
      
      !**************************************************************************************************************************
      ! calculate the matrix element for the crossing point number iC
      !**************************************************************************************************************************
      subroutine calculateMatrixElement(cnt1,cnt2,iC)
        use inputParameters
        use physicalConstants, only : i1
        type(cnt), intent(in) :: cnt1,cnt2
        integer, intent(in) :: iC
        integer :: ix1,ix2
        integer :: iKcm
        integer :: mu1,mu2
        integer :: ikr1, ikr2
        integer :: iv, ivp
        integer :: imu1,imu2
        integer :: is,isp
        complex*16 :: matrixElement
        
        matrixElement = (0.d0,0.d0)
        
        ix1 = finalCrossingPoints(iC,1)
        ix2 = finalCrossingPoints(iC,2)
        iKcm = finalCrossingPoints(iC,3)
        
        mu1 = cnt1.min_sub(i_sub1)
        mu2 = cnt2.min_sub(i_sub2)
        
        ! Calculate matrix element for mu1 = + and mu2 = +
        imu1 = 1
        imu2 = 1
        
        do ikr1 = cnt1.ikr_low, cnt1.ikr_high
          do ikr2 = cnt2.ikr_low, cnt2.ikr_high
              do iv = 1, cnt1.Nu
                do ivp = 1,cnt2.Nu
                  is = 1
                  isp = 1
                  matrixElement = matrixElement + conjg(cnt1.Cc(imu1,ikr1+iKcm,is))*cnt1.Cv(imu1,ikr1-iKcm,is)*cnt2.Cc(imu2,ikr2+iKcm,isp)*conjg(cnt2.Cv(imu2,ikr2-iKcm,isp))* &
                      exp(i1*dcmplx(2.d0*(-dot_product(cnt1.K1,cnt1.posA(iv,:))+dot_product(cnt2.K1,cnt1.posA(iv,:)))))
                
                end do  
              end do
          end do
        end do
        
      end subroutine calculateMatrixElement
      
    
end module forsterClass