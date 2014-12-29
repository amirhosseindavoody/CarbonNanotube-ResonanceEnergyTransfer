!*******************************************************************************
! This subroutines sets the physical constants of the simulation
!*******************************************************************************
subroutine fnPhysConst()
  use comparams
  implicit none
  
  E_th=E_th*eV
  
  
  write(fh1,10) "TIGHTBINDING PARAMETERS"
  write(fh1,11) "e2p [eV]=",e2p/eV
  write(fh1,11) "t0 [eV]=",t0/eV
  write(fh1,11) "s0 [eV]=",s0
  
10 FORMAT (A100)
11 FORMAT (A10,E16.8)   
12 FORMAT (A10,I3) 
  return
end