!*******************************************************************************
! This program calculates the excitonic band structure of single wall carbon nanotubes through simple tight-binding method.
! Amirhossein Davoody
! Last modified: 3/17/2014
!*******************************************************************************

program cntForsterEnergyTransfer
  use cnt_class
  use inputParameters
  implicit none
  
  type (cnt) :: firstCNT, secondCNT
  
  
  
  firstCNT = cnt( n_ch1, m_ch1, nkg1)
  secondCNT = cnt (n_ch2, m_ch2, nkg2)
  
  call firstCNT.calculateBands(i_sub1, E_th, Kcm_max)
  call firstCNT.printProperties()
  
  call secondCNT.calculateBands(i_sub2, E_th, Kcm_max)
  call secondCNT.printProperties()
  

  print *,'Finish!!!!'
  pause
  
end program cntForsterEnergyTransfer

