!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Declaration of global variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module physicalConstants
  implicit none
  
  real*8, parameter :: pi=3.141592d0
  complex*16, parameter :: i1=(0.d0,1.d0)
  
  !Physical constants
  real*8, parameter :: eV=1.6d-19 ![Joules]
  real*8, parameter :: hb=6.5d-16*eV ![eV.s]
  
  real*8, parameter :: a_cc=1.42d-10 !carbon-carbon distance [meters] 
  real*8, parameter :: a_l=dsqrt(3.d0)*a_cc !lattice constants
  
  real*8, parameter :: e2p = 0.d0 !tight binding constants
  real*8, parameter :: t0 = 2.7d0*eV  !tight binding constants
  real*8, parameter :: s0 = 0.d0 !tight binding constants
  
  real*8, parameter :: Upp = 11.3d0*eV !constant used in the Ohno potential
  real*8, parameter :: eps0 = 8.85d-12 !permittivity of free space
  real*8, parameter :: q0 = 1.6d-19 !charge of electron
  real*8, parameter :: kappa = 2.d0 !dielectric constant due to core electrons in CNT
  
end module physicalConstants