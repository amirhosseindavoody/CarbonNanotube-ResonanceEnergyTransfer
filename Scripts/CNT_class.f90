module cntClass
    use physicalConstants
    implicit none
    
    private
    
    public  :: cnt
    
    type cnt
      integer, public :: n_ch,m_ch !chiral vector parameters
      integer, public :: i_sub !subband index used in exciton energy calculation
    
      !Geometrical properties
      real*8, dimension(2) :: a1,a2,b1,b2,ch_vec,t_vec,aCC_vec
      real*8 :: len_ch,radius
      integer :: Nu
      real*8, dimension(:,:), allocatable, public :: posA,posB,posAA,posBB,posAB,posBA
      real*8, dimension(:,:), allocatable, public :: posA3, posB3
			real*8, dimension(:,:,:), allocatable, public :: pos2d, pos3d
    
      !Reciprocal lattice properties
      real*8, public :: dk
      real*8, dimension(2) :: K1, K2
      
      !CNT band structure properties
      integer, dimension(:), allocatable :: min_sub
      integer, public :: ikc_max, ikc_min, ik_max, ik_min, iKcm_max, iKcm_min, ik_high, ik_low, ikr_high, ikr_low, iq_max, iq_min
      
      !CNT self energy and tight binding coefficients
      real*8, dimension(:,:,:), allocatable, public :: Ek  !Tight-binding energy
      real*8, dimension(:,:,:), allocatable:: Sk  !Self-energy 
      complex*16, dimension(:,:,:), allocatable:: Cc,Cv !Tight-binding wavefunction coefficients
      
      !Exciton wavefunction and energies
      real*8, dimension(:,:), allocatable, public :: Ex_A1, Ex0_A2, Ex1_A2 !the first index is subband, the second index is iKcm
      complex*16, dimension(:,:,:), allocatable, public :: Psi_A1, Psi0_A2, Psi1_A2 !the first index is ikr, the scond index is the subband, the third index is iKcm
      integer, public :: nX
      
    contains
        procedure :: calculateBands
        procedure :: printProperties
    end type cnt
    
    interface cnt
      module procedure init_cnt
    end interface
      
    contains
      
      !**************************************************************************************************************************
      ! initialize CNT by calculating its geometrical properties
      !**************************************************************************************************************************
      type(cnt) function init_cnt(n_ch,m_ch,nkg)
        use smallFunctions
        integer, intent(in) :: n_ch, m_ch
        integer, intent(in) :: nkg
        integer :: dR = 1
        integer :: t1,t2
        real*8 :: cosTh, sinTh
        real*8, dimension(2,2) :: Rot
        integer :: i,j,k
        
        init_cnt.n_ch = n_ch
        init_cnt.m_ch = m_ch
        
        ! unit vectors and reciprocal lattice vectors.
        init_cnt.a1=(/dsqrt(3.d0)/2.d0*a_l, 1.d0/2.d0*a_l/)
        init_cnt.a2=(/dsqrt(3.d0)/2.d0*a_l, -1.d0/2.d0*a_l/)
        init_cnt.b1=(/1.d0/dsqrt(3.d0)*2.d0*pi/a_l, +1.d0*2.d0*pi/a_l/)
        init_cnt.b2=(/1.d0/dsqrt(3.d0)*2.d0*pi/a_l, -1.d0*2.d0*pi/a_l/)
        init_cnt.aCC_vec=1.d0/3.d0*(init_cnt.a1+init_cnt.a2)
        
        ! calculate chirality and translational vectors of CNT unit cell.
        init_cnt.ch_vec=dble(n_ch)*init_cnt.a1+dble(m_ch)*init_cnt.a2
        init_cnt.len_ch=a_l*dsqrt(dble(n_ch)**2+dble(m_ch)**2+dble(n_ch)*dble(m_ch))
        init_cnt.radius=init_cnt.len_ch/2.d0/pi
  
        call gcd(dR,2*n_ch+m_ch,2*m_ch+n_ch)
  
        t1=(2*m_ch+n_ch)/dR
        t2=-(2*n_ch+m_ch)/dR
  
        init_cnt.t_vec=dble(t1)*init_cnt.a1+dble(t2)*init_cnt.a2
  
        init_cnt.Nu=2*(n_ch**2+m_ch**2+n_ch*m_ch)/dR
 
 
        ! rotate basis vectors so that ch_vec is along x-axis.
        cosTh=init_cnt.ch_vec(1)/norm2(init_cnt.ch_vec)
        sinTh=init_cnt.ch_vec(2)/norm2(init_cnt.ch_vec)
        Rot=reshape((/ cosTh, -sinTh , sinTh, cosTh /), (/2,2/))
        init_cnt.ch_vec=matmul(Rot,init_cnt.ch_vec)
        init_cnt.t_vec=matmul(Rot,init_cnt.t_vec)
        init_cnt.a1=matmul(Rot,init_cnt.a1)
        init_cnt.a2=matmul(Rot,init_cnt.a2)
        init_cnt.b1=matmul(Rot,init_cnt.b1)
        init_cnt.b2=matmul(Rot,init_cnt.b2)
        init_cnt.aCC_vec=matmul(Rot,init_cnt.aCC_vec)
  
        ! calculate reciprocal lattice of CNT.
        init_cnt.dk=norm2(init_cnt.b1)/(dble(nkg)-1.d0)
        init_cnt.K1=(-t2*init_cnt.b1+t1*init_cnt.b2)/(dble(init_cnt.Nu))
        init_cnt.K2=(dble(m_ch)*init_cnt.b1-dble(n_ch)*init_cnt.b2)/dble(init_cnt.Nu)
        init_cnt.K2=init_cnt.K2/norm2(init_cnt.K2)
  
        ! calculate coordinates of atoms in the unwarped CNT unit cell.
        allocate(init_cnt.posA(init_cnt.Nu,2))
        allocate(init_cnt.posB(init_cnt.Nu,2))
    
        k=0
        do i=0,t1+n_ch
          do j=t2,m_ch
            if ((dble(t2)/dble(t1)*i .le. j) .and. (dble(m_ch)/dble(n_ch)*i .ge. j) .and. (dble(t2)/dble(t1)*(i-n_ch) .gt. (j-m_ch)) .and. (dble(m_ch)/dble(n_ch)*(i-t1) .lt. (j-t2))) then
              k=k+1
              init_cnt.posA(k,1)=dble(i)*init_cnt.a1(1)+dble(j)*init_cnt.a2(1)
              init_cnt.posA(k,2)=dble(i)*init_cnt.a1(2)+dble(j)*init_cnt.a2(2)
              init_cnt.posB(k,1)=init_cnt.posA(k,1)+init_cnt.aCC_vec(1)
              init_cnt.posB(k,2)=init_cnt.posA(k,2)+init_cnt.aCC_vec(2)
        
              if (init_cnt.posA(k,1) .gt. init_cnt.ch_vec(1)) init_cnt.posA(k,1)=init_cnt.posA(k,1)-init_cnt.ch_vec(1);
              if (init_cnt.posA(k,1) .lt. 0) init_cnt.posA(k,1)=init_cnt.posA(k,1)+init_cnt.ch_vec(1);
              if (init_cnt.posA(k,2) .gt. init_cnt.t_vec(2)) init_cnt.posA(k,2)=init_cnt.posA(k,2)-init_cnt.t_vec(2);
              if (init_cnt.posA(k,2) .lt. 0) init_cnt.posA(k,2)=init_cnt.posA(k,2)+init_cnt.t_vec(2);
          
              if (init_cnt.posB(k,1) .gt. init_cnt.ch_vec(1)) init_cnt.posB(k,1)=init_cnt.posB(k,1)-init_cnt.ch_vec(1);
              if (init_cnt.posB(k,1) .lt. 0) init_cnt.posB(k,1)=init_cnt.posB(k,1)+init_cnt.ch_vec(1);
              if (init_cnt.posB(k,2) .gt. init_cnt.t_vec(2)) init_cnt.posB(k,2)=init_cnt.posB(k,2)-init_cnt.t_vec(2);
              if (init_cnt.posB(k,2) .lt. 0) init_cnt.posB(k,2)=init_cnt.posB(k,2)+init_cnt.t_vec(2);
              
            endif
          enddo
				enddo
				
				allocate(init_cnt.pos2d(2,init_cnt.Nu,2))
				init_cnt.pos2d(1,:,:) = init_cnt.posA(:,:)
				init_cnt.pos2d(2,:,:) = init_cnt.posB(:,:)
    
        if (k .ne. init_cnt.Nu) stop "*** Error in calculating atom positions ***"
  
        ! calculate distances between atoms in a warped CNT unit cell.
        allocate(init_cnt.posAA(init_cnt.Nu,2))
        allocate(init_cnt.posAB(init_cnt.Nu,2))
        allocate(init_cnt.posBA(init_cnt.Nu,2))
        allocate(init_cnt.posBB(init_cnt.Nu,2))
  
        do i=1,init_cnt.Nu
          init_cnt.posAA(i,:)=init_cnt.posA(i,:)-init_cnt.posA(1,:)
          init_cnt.posAB(i,:)=init_cnt.posA(i,:)-init_cnt.posB(1,:)
          init_cnt.posBA(i,:)=init_cnt.posB(i,:)-init_cnt.posA(1,:)
          init_cnt.posBB(i,:)=init_cnt.posB(i,:)-init_cnt.posB(1,:)
          if (init_cnt.posAA(i,1) .gt. init_cnt.ch_vec(1)/2.d0) init_cnt.posAA(i,1)=init_cnt.posAA(i,1)-init_cnt.ch_vec(1)
          if (init_cnt.posAB(i,1) .gt. init_cnt.ch_vec(1)/2.d0) init_cnt.posAB(i,1)=init_cnt.posAB(i,1)-init_cnt.ch_vec(1)
          if (init_cnt.posBA(i,1) .gt. init_cnt.ch_vec(1)/2.d0) init_cnt.posBA(i,1)=init_cnt.posBA(i,1)-init_cnt.ch_vec(1)
          if (init_cnt.posBB(i,1) .gt. init_cnt.ch_vec(1)/2.d0) init_cnt.posBB(i,1)=init_cnt.posBB(i,1)-init_cnt.ch_vec(1)
        end do
        
        ! calculate coordinates of atoms in 3D unit cell
        allocate(init_cnt.posA3(init_cnt.Nu,3))
        allocate(init_cnt.posB3(init_cnt.Nu,3))
        
        do i=1,init_cnt.Nu
          init_cnt.posA3(i,1) = init_cnt.radius*sin(2*pi*init_cnt.posA(i,1)/init_cnt.len_ch)
          init_cnt.posA3(i,2) = init_cnt.posA(i,2)
          init_cnt.posA3(i,3) = -init_cnt.radius*cos(2*pi*init_cnt.posA(i,1)/init_cnt.len_ch)
        end do
        
        do i=1,init_cnt.Nu
          init_cnt.posB3(i,1) = init_cnt.radius*sin(2*pi*init_cnt.posB(i,1)/init_cnt.len_ch)
          init_cnt.posB3(i,2) = init_cnt.posB(i,2)
          init_cnt.posB3(i,3) = -init_cnt.radius*cos(2*pi*init_cnt.posB(i,1)/init_cnt.len_ch)
				end do
				
				allocate(init_cnt.pos3d(2,init_cnt.Nu,3))
				init_cnt.pos3d(1,:,:) = init_cnt.posA3(:,:)
				init_cnt.pos3d(2,:,:) = init_cnt.posB3(:,:)
        
      end function init_cnt
    
    
      !**************************************************************************************************************************
      ! calculate band structure of the CNT
      !**************************************************************************************************************************
      subroutine calculateBands(self, i_sub, E_th, Kcm_max)
        class(cnt), intent(inout) :: self
        integer, intent(in) :: i_sub
        real*8, intent(in) :: E_th, Kcm_max
        
        integer :: nkc, imin_sub
        integer :: i,j,mu,ik,tmpi
        integer, dimension(:), allocatable :: min_loc
        real*8 :: tmpr
        real*8, dimension(2) :: k,E1_tmp,E2_tmp
        real*8, dimension(:), allocatable :: k_vec,min_energy
        real*8, dimension(:,:,:), allocatable :: E_k
        complex*16, dimension(:,:,:), allocatable :: Cc_k,Cv_k
        complex*16, dimension(2) :: Cc_tmp,Cv_tmp
  
        ! set the value of i_sub
        self.i_sub = i_sub
        
        ! calculate CNT energy dispersion.
        self.ikc_max=floor(pi/norm2(self.t_vec)/self.dk)
        self.ikc_min=-self.ikc_max
        nkc=2*self.ikc_max+1
  
        allocate(k_vec(self.ikc_min:self.ikc_max))
        allocate(E_k(1-self.Nu/2:self.Nu/2,self.ikc_min:self.ikc_max,2))
        allocate(Cc_k(1-self.Nu/2:self.Nu/2,self.ikc_min:self.ikc_max,2))
        allocate(Cv_k(1-self.Nu/2:self.Nu/2,self.ikc_min:self.ikc_max,2))
        allocate(min_loc(0:self.Nu/2))
  
        do ik=self.ikc_min,self.ikc_max
          k_vec(ik)=dble(ik)*self.dk
        end do
  
        do mu=1-self.Nu/2,self.Nu/2
          do ik=self.ikc_min,self.ikc_max
            k = dble(mu) * self.K1 + dble(ik) * self.dk * self.K2
            call grapheneEnergy(self, E_k(mu,ik,:), Cc_k(mu,ik,:), Cv_k(mu,ik,:), k)
          enddo
        enddo
  
  
        ! find the subbands with a minimum energy.
        min_loc=minloc(E_k(0:self.Nu/2,:,1),2)
        imin_sub=count((min_loc .lt. nkc) .and. (min_loc .gt. 1))
        allocate(self.min_sub(imin_sub))
        allocate(min_energy(imin_sub))
  
        ! store the value of mu for subbands with minimums in the variable min_sub
        i=1
        do mu=0,self.Nu/2
          if ((min_loc(mu) .gt. 1) .and. (min_loc(mu) .lt. nkc)) then
            self.min_sub(i)=mu
            min_energy(i)=minval(E_k(mu,:,1))
            i=i+1
          end if
        end do
  
        ! sort the subbands
        do i=imin_sub,2,-1
          do j=i-1,1,-1
            if (min_energy(i) .lt. min_energy(j)) then
              tmpr=min_energy(i)
              tmpi=self.min_sub(i)
              min_energy(i)=min_energy(j)
              self.min_sub(i)=self.min_sub(j)
              min_energy(j)=tmpr
              self.min_sub(j)=tmpi
            end if    
          end do
        end do
  
        ! find the max k-index that energy is below threshold energy (E_th).
        ik=0
        E1_tmp=(/ min_energy(self.i_sub),0 /)
        E2_tmp=(/ min_energy(self.i_sub),0 /)
        do while ((min(E1_tmp(1),E2_tmp(1))-min_energy(self.i_sub)) .le. E_th )
          k=dble(self.min_sub(self.i_sub))*self.K1+dble(ik)*self.dk*self.K2
          call grapheneEnergy(self,E1_tmp,Cc_tmp,Cv_tmp,k)
          k=dble(self.min_sub(self.i_sub))*self.K1-dble(ik)*self.dk*self.K2
          call grapheneEnergy(self,E2_tmp(:),Cc_tmp,Cv_tmp,k)
          ik=ik+1
        enddo
  
        ! set the index boundaries for some arrays and kernels.
        self.ik_max=ik                              !the higher limit of k-vector that is below E_th
        self.ik_min=-ik                             !the lower limit of k-vector that is below E_th
        self.iKcm_max=floor(Kcm_max/self.dk)        !the higher limit of center of mass wave vector that we calculate
        self.iKcm_min = - self.iKcm_max             !the lower limit of center of mass wave vector that we calculate
        self.ikr_high=self.iKcm_max-self.ik_min     !the maximum index that the relative wavenumber in the entire simulation.
        self.ikr_low=-self.ikr_high                 !the minimum index that the relative wavenumber in the entire simulation.
        self.ik_high=self.ikr_high+self.iKcm_max    !the maximum index that the wavenumber in the entire simulation.
        self.ik_low=-self.ik_high                   !the minimum index that the wavenumber in the entire simulation.
        self.iq_max=2*self.ikr_high                 !the higher limit of the index in v_FT and esp_q
        self.iq_min=-self.iq_max                    !the lower limit of the index in v_FT and esp_q
        
        ! calculate the tight-binding energies and coefficients.
        allocate(self.Ek(2,self.ik_low:self.ik_high,2))        
        allocate(self.Cc(2,self.ik_low:self.ik_high,2))
        allocate(self.Cv(2,self.ik_low:self.ik_high,2))
  
        
        do ik=self.ik_low,self.ik_high
          mu=self.min_sub(i_sub) !first band
          k=dble(mu)*self.K1+dble(ik)*self.dk*self.K2
          call grapheneEnergy(self,self.Ek(1,ik,:),self.Cc(1,ik,:),self.Cv(1,ik,:),k)

          mu=-self.min_sub(i_sub) !second band
          k=dble(mu)*self.K1+dble(ik)*self.dk*self.K2
          call grapheneEnergy(self,self.Ek(2,ik,:),self.Cc(2,ik,:),self.Cv(2,ik,:),k)
        enddo
      
      end subroutine calculateBands
    
      !**************************************************************************************************************************
      ! private subroutine to calculate Bloch functions and energy in graphene
      !**************************************************************************************************************************
      subroutine grapheneEnergy(currCNT,E,Cc,Cv,k)
        type(cnt), intent(in) :: currCNT
        complex*16 :: f_k
        real*8, dimension(2), intent(in) :: k
        real*8, dimension(2), intent(out) :: E
        complex*16, dimension(2), intent(out) :: Cv
        complex*16, dimension(2), intent(out) :: Cc
  
        f_k=exp(i1*dot_product(k,(currCNT.a1+currCNT.a2)/3.d0))+exp(i1*dot_product(k,(currCNT.a1-2.d0*currCNT.a2)/3.d0))+exp(i1*dot_product(k,(currCNT.a2-2.d0*currCNT.a1)/3.d0))  
  
        E(1)=+t0*abs(f_k)
        E(2)=-t0*abs(f_k)
  
        Cc(1)=+1.d0/sqrt(2.d0)
        Cc(2)=+1.d0/sqrt(2.d0)*conjg(f_k)/abs(f_k)
        Cv(1)=+1.d0/sqrt(2.d0)
        Cv(2)=-1.d0/sqrt(2.d0)*conjg(f_k)/abs(f_k)
      end subroutine grapheneEnergy

      !**************************************************************************************************************************
      ! print properties of CNT to the console
      !**************************************************************************************************************************
      subroutine printProperties(self)
        class(cnt), intent(in) :: self
        print *, "chirality = (", self.n_ch, ",", self.m_ch,")"
        print *, "radius = ", self.radius
        print *, "t_vec = ", self.t_vec
      end subroutine printProperties
    
end module cntClass