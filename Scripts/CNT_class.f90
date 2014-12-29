module cnt_class
    implicit none
    
    private
    
    public  :: cnt
    
    type cnt
      integer :: n_ch,m_ch !chiral vector parameters
      integer :: nkg !reciprocal space mesh size in graphene
      integer :: nr !number of CNT unit cells in real space
      integer :: i_sub !subband index used in exciton energy calculation
    
      !Geometrical properties
      real*8, dimension(2) :: a1,a2,b1,b2,ch_vec,t_vec,aCC_vec
      real*8 :: len_ch,radius
      integer :: dR,Nu,MC
      integer :: t1,t2
      real*8, dimension(:,:), allocatable :: posA,posB,posAA,posBB,posAB,posBA
    
      !Reciprocal lattice properties
      real*8 :: dk
      real*8, dimension(2) :: K1, K2
      
      !CNT band structure properties
      integer, dimension(:), allocatable :: min_sub
      integer :: ikc_max, ikc_min, ik_max, ik_min, iKcm_max, iKcm_min, ik_high, ik_low, ikr_high, ikr_low, iq_max, iq_min
      
      !CNT self energy and tight binding coefficients
      real*8, dimension(:,:,:), allocatable:: Ek  !Tight-binding energy
      real*8, dimension(:,:,:), allocatable:: Sk  !Self-energy 
      complex*16, dimension(:,:,:), allocatable:: Cc,Cv !Tight-binding wavefunction coefficients
      
    end type cnt
    
end module cnt_class