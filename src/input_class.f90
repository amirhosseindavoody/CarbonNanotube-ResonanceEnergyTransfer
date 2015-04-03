!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Declaration of input parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
module input_class
	implicit none
	private
	public :: inputCNT

contains	
	!**************************************************************************************************************************
	! read the CNT information from the log file of exciton calculation which is done previously
	!**************************************************************************************************************************

	subroutine inputCNT(currcnt,dirname)
		use physicalConstants, only: eV
		use cnt_class, only: cnt, init_cnt, calculateBands
		use output_module, only: writeLog, logInput

		type(cnt) :: currcnt
		character(len=100) :: dirname
		character(len=200) :: buffer, label
		integer :: istat=0
		integer :: ios=0
		integer :: pos=0
		integer, parameter :: nparam=10
		integer :: iparam=0
		integer :: n_ch, m_ch, nkg, nr, i_sub, nX
		real*8 :: E_th, Kcm_max, kappa, Ckappa

		istat=chdir(dirname)
		if (istat .ne. 0) then
			write(*,*) "Directory does not exist:", dirname
			call exit()
		end if

		open(unit=100,file='log.dat',status="old")

		iparam=0
		ios=0
		do while ((ios == 0) .and. (iparam .lt. nparam))
			read (100,'(A)',iostat=ios) buffer
			if (ios == 0) then
				pos = scan(buffer,'=')
				label = buffer(1:pos-1)
				buffer = buffer(pos+1:)
				label = adjustl(label)
				buffer = adjustl(buffer)

				select case (trim(label))
				case ('n_ch')
					read(buffer, *, iostat=ios) n_ch
! 					write(logInput,*) 'Read n_ch: ', n_ch
! 					call writeLog()
					iparam = iparam+1
				case ('m_ch')
					read(buffer, *, iostat=ios) m_ch
! 					write(logInput,*) 'Read m_ch: ', m_ch
! 					call writeLog()
					iparam = iparam+1
				case ('nkg')
					read(buffer, *, iostat=ios) nkg
! 					write(logInput,*) 'Read nkg: ', nkg
! 					call writeLog()
					iparam = iparam+1
				case ('nr')
					read(buffer, *, iostat=ios) nr
! 					write(logInput,*) 'Read nr: ', nr
! 					call writeLog()
					iparam = iparam+1
				case ('E_th')
					read(buffer, *, iostat=ios) E_th
! 					write(logInput,*) 'Read E_th: ', E_th
! 					call writeLog()
					iparam = iparam+1
					E_th=E_th*eV
				case ('Kcm_max')
					read(buffer, *, iostat=ios) Kcm_max
! 					write(logInput,*) 'Read Kcm_max: ', Kcm_max
! 					call writeLog()
					Kcm_max = Kcm_max*1.d9
					iparam = iparam+1
				case ('i_sub')
					read(buffer, *, iostat=ios) i_sub
! 					write(logInput,*) 'Read i_sub: ', i_sub
! 					call writeLog()
					iparam = iparam+1
				case ('Ckappa')
					read(buffer, *, iostat=ios) Ckappa
! 					write(logInput,*) 'Read Ckappa: ', Ckappa
! 					call writeLog()
					iparam = iparam+1
				case ('kappa')
					read(buffer, *, iostat=ios) kappa
! 					write(logInput,*) 'Read kappa: ', kappa
! 					call writeLog()
					iparam = iparam+1
				case ('nX')
					read(buffer, *, iostat=ios) nX
! 					write(logInput,*) 'Read nX: ', nX
! 					call writeLog()
					iparam = iparam+1
				end select
			end if
		end do

		close(100)

		if (iparam .lt. nparam) then
			write(*,*) "Error in reading parameters!"
			call exit()
		end if

		! create the cnt object and calculate bands and load exciton wavefunction
		currcnt = init_cnt( n_ch, m_ch, nkg, nX, kappa, Ckappa)
		call calculateBands(currcnt,i_sub, E_th, Kcm_max)
		call inputExciton(currcnt)
		
		!change the directory back
		istat=chdir('..')
		if (istat .ne. 0) then
			write(*,*) "Directory not changed to root!"
			call exit()
		end if

		return
	end subroutine inputCNT

	!**************************************************************************************************************************
	! load the exciton wavefunction and energies from the ExcitonEnergy calculation
	!**************************************************************************************************************************
	
	subroutine inputExciton(currcnt)
		use cnt_class, only: cnt
		type(cnt), intent(inout) :: currcnt
		integer :: iX, iKcm, ikr
		real*8 :: tmpr1, tmpr2, tmpr3
		
		open(unit=100,file='Ex_A1.dat',status="old")
		open(unit=101,file='Ex0_A2.dat',status="old")
		open(unit=102,file='Ex1_A2.dat',status="old")
		open(unit=103,file='Psi_A1.dat',status="old")
		open(unit=104,file='Psi0_A2.dat',status="old")
		open(unit=105,file='Psi1_A2.dat',status="old")
		
		allocate(currcnt%Ex_A1(1:currcnt%nX,currcnt%iKcm_min:currcnt%iKcm_max))
		allocate(currcnt%Ex0_A2(1:currcnt%nX,currcnt%iKcm_min:currcnt%iKcm_max))
		allocate(currcnt%Ex1_A2(1:currcnt%nX,currcnt%iKcm_min:currcnt%iKcm_max))
		allocate(currcnt%Psi_A1(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%nX,currcnt%iKcm_min:currcnt%iKcm_max))
		allocate(currcnt%Psi0_A2(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%nX,currcnt%iKcm_min:currcnt%iKcm_max))
		allocate(currcnt%Psi1_A2(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%nX,currcnt%iKcm_min:currcnt%iKcm_max))

		allocate(currcnt%Ex_t(1:currcnt%nX,currcnt%iKcm_min:currcnt%iKcm_max))
		allocate(currcnt%Psi_t(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%nX,currcnt%iKcm_min:currcnt%iKcm_max))
		
		do iKcm=currcnt%iKcm_min,currcnt%iKcm_max
			do iX=1,currcnt%nX
				read(100,'(E16.8)', advance='no') currcnt%Ex_A1(iX,iKcm)
				read(101,'(E16.8)', advance='no') currcnt%Ex0_A2(iX,iKcm)
				read(102,'(E16.8)', advance='no') currcnt%Ex1_A2(iX,iKcm)
				do ikr=currcnt%ikr_low,currcnt%ikr_high
					read(103,'(E16.8,E16.8)', advance='no') currcnt%Psi_A1(ikr,iX,iKcm)
					read(104,'(E16.8,E16.8)', advance='no') currcnt%Psi0_A2(ikr,iX,iKcm)
					read(105,'(E16.8,E16.8)', advance='no') currcnt%Psi1_A2(ikr,iX,iKcm)
				enddo
			enddo
			
			read(100,'(E16.8)')
			read(101,'(E16.8)')
			read(102,'(E16.8)')
			read(103,'(E16.8)')
			read(104,'(E16.8)')
			read(105,'(E16.8)')
		enddo
		
		
		close(100)
		close(101)
		close(102)
		close(103)
		close(104)
		close(105)
				
		!make sure the exciton wavefunctions are normalized
		do iX=1,currcnt%nX
			do iKcm=currcnt%iKcm_min,currcnt%iKcm_max
				tmpr1 = 0.d0
				tmpr2 = 0.d0
				tmpr3 = 0.d0
				do ikr=currcnt%ikr_low,currcnt%ikr_high
					tmpr1 = tmpr1 + dble(conjg(currcnt%Psi_A1(ikr,iX,iKcm))*currcnt%Psi_A1(ikr,iX,iKcm))
					tmpr2 = tmpr2 + dble(conjg(currcnt%Psi0_A2(ikr,iX,iKcm))*currcnt%Psi0_A2(ikr,iX,iKcm))
					tmpr3 = tmpr3 + dble(conjg(currcnt%Psi1_A2(ikr,iX,iKcm))*currcnt%Psi1_A2(ikr,iX,iKcm))
				enddo
				currcnt%Psi_A1(:,iX,iKcm) = currcnt%Psi_A1(:,iX,iKcm) / dcmplx(sqrt(tmpr1))
				currcnt%Psi0_A2(:,iX,iKcm) = currcnt%Psi0_A2(:,iX,iKcm) / dcmplx(sqrt(tmpr2))
				currcnt%Psi1_A2(:,iX,iKcm) = currcnt%Psi1_A2(:,iX,iKcm) / dcmplx(sqrt(tmpr3))
			enddo
		enddo

		return
	end subroutine inputExciton

	
end module input_class