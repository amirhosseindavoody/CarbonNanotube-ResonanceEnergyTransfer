!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Declaration of input parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
module inputParameters
  use physicalConstants, only: eV
  implicit none
  
  real*8 :: Temperature = 300 !Temperature of the system in Kelvin units
  real*8 :: E_th = 1.5d0 * eV !Threshold energy that is used as a cutoff on how large the limits of k-points included in the ExcitonEnergy program are.
  real*8 :: Kcm_max = 1.5d9 !Maximum value of center of mass k-vector
  integer :: nkg = 501 !number of k points along b1 vector in graphene Brillouin zone.
  integer :: nr = 200 !number of cnt unit cells in real space that are considered in calculating the overlap integrals in ExcitonEnergy program.
  
  
  integer :: n_ch1 = 7, m_ch1 = 6, i_sub1 = 1
  integer :: n_ch2 = 7, m_ch2 = 5, i_sub2 = 1
  
end module inputParameters