!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Declaration of input parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
module inputParameters
  use physicalConstants, only: eV
  implicit none
  
  real*8 :: E_th = 1.5d0 * eV !Threshold energy
  real*8 :: Kcm_max = 1.5d9 !Maximum value of center of mass k-vector
  integer :: nkg = 501, nr = 200
  
  integer :: n_ch1 = 7, m_ch1 = 6, i_sub1 = 1
  integer :: n_ch2 = 7, m_ch2 = 5, i_sub2 = 1
  
end module inputParameters