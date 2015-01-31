!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Declaration of input parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
module inputParameters
  use physicalConstants, only: eV, pi
  implicit none
  
  real*8 :: Temperature = 300 !Temperature of the system in Kelvin units
  real*8 :: E_th = 1.5d0 * eV !Threshold energy that is used as a cutoff on how large the limits of k-points included in the ExcitonEnergy program are.
  real*8 :: Kcm_max = 1.5d9 !Maximum value of center of mass k-vector
  integer :: nkg = 501 !number of k points along b1 vector in graphene Brillouin zone.
  integer :: nr = 200 !number of cnt unit cells in real space that are considered in calculating the overlap integrals in ExcitonEnergy program.
  
	real*8 :: ppLen = 20.d-10 !length per perpendicular tubes
	
	real*8 :: c2cDistance !center to center distance between parallel carbon nanotubes
	real*8 :: c2cMin = 1.2d-10
	real*8 :: c2cMax = 10.d-10
	integer :: nc2c = 100
	integer :: ic2c
	real*8 :: dc2c
	
	real*8 :: theta
	real*8 :: thetaMax = pi/2.d0
	integer :: nTheta = 200
	integer :: iTheta
	real*8 :: dTheta
	
	
	real*8, dimension(:,:,:), allocatable :: transitionRate
  
  
  integer :: n_ch1 = 7, m_ch1 = 6, i_sub1 = 1
  integer :: n_ch2 = 8, m_ch2 = 6, i_sub2 = 1
  
end module inputParameters