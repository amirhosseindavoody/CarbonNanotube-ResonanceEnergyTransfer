module dataClass
    use cntClass
    implicit none
    
    private

    public  :: saveCNTProperties, loadExcitonWavefunction, saveTransitionPoints, saveDOS
    
    contains
      !**************************************************************************************************************************
      ! save the properties of carbon nanotubes that are calculated in cnt_class
      !**************************************************************************************************************************
      subroutine saveCNTProperties(currCNT)
        use inputParameters
        use smallFunctions
        type(cnt), intent(in) :: currCNT
        character(len=200) :: command
        integer(4) :: istat
        integer :: i
        
        !change the directory to that of the CNT
				write(command,'("mkdir ",A100)') currcnt%excitonDirectory
        call execute_command_line(command,wait=.true.,exitstat=i)
    
        istat=chdir(currcnt%excitonDirectory)
        if (istat .ne. 0) then
            write(logInput,*) "Directory did not changed!!!"
            call writeLog()
            call exit()
        end if
        
        ! save the position of A-type and B-type carbon atoms
        open(unit=100, file='posA.dat',status='unknown')
        open(unit=101, file='posB.dat',status='unknown')
        do i=1,size(currcnt%posA,1)
          write(100,14) currcnt%posA(i,1), currcnt%posA(i,2)
          write(101,14) currcnt%posB(i,1), currcnt%posB(i,2)
        end do
        close(100)
        close(101)
        
        open(unit=100, file="CondBand_Sub.dat", status="unknown")
        open(unit=101, file="ValeBand_Sub.dat", status="unknown")
        
        do i=currcnt%ik_low,currcnt%ik_high
          write(100,10, advance='no') dble(i)*currcnt%dk
          write(101,10, advance='no') dble(i)*currcnt%dk
        end do
        write(100,10)
        write(101,10)
        
        do i=currcnt%ik_low,currcnt%ik_high
          write(100,10, advance='no') currcnt%Ek(1,i,1)
          write(101,10, advance='no') currcnt%Ek(1,i,2)
        end do
        write(100,10)
        write(101,10)
        
        do i=currcnt%ik_low,currcnt%ik_high
          write(100,10, advance='no') currcnt%Ek(2,i,1)
          write(101,10, advance='no') currcnt%Ek(2,i,2)
        end do
        
        close(100)
        close(101)
        
        !change the directory back
        istat=chdir('..')
        if (istat .ne. 0) then
            write(logInput,*) "Directory did not changed!!!"
            call writeLog()
            call exit()
        end if
        
        10 FORMAT (E16.8)
        14 FORMAT (E16.8,E16.8) 
      end subroutine saveCNTProperties
      
      !**************************************************************************************************************************
      ! load the exciton wavefunction and energies from the ExcitonEnergy calculation
      !**************************************************************************************************************************
      subroutine loadExcitonWavefunction(currCNT)
        use inputParameters
        use smallFunctions
        type(cnt), intent(inout) :: currCNT
        integer(4) :: istat
        integer :: iX, iKcm, ikr
				real*8 :: tmpr1, tmpr2, tmpr3
				
        call loadMiscData(currcnt)

        !change the directory to that of the CNT
        istat=chdir(currcnt%excitonDirectory)
        if (istat .ne. 0) then
            write(logInput,*) "Directory did not changed!!!"
            call writeLog()
            call exit()
        end if
        
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
        
        do iKcm=currcnt%iKcm_min,currcnt%iKcm_max
          do iX=1,currcnt%nX
            read(100,10, advance='no') currcnt%Ex_A1(iX,iKcm)
            read(101,10, advance='no') currcnt%Ex0_A2(iX,iKcm)
            read(102,10, advance='no') currcnt%Ex1_A2(iX,iKcm)
            do ikr=currcnt%ikr_low,currcnt%ikr_high
              read(103,11, advance='no') currcnt%Psi_A1(ikr,iX,iKcm)
              read(104,11, advance='no') currcnt%Psi0_A2(ikr,iX,iKcm)
              read(105,11, advance='no') currcnt%Psi1_A2(ikr,iX,iKcm)
            enddo
          enddo
          
          read(100,10)
          read(101,10)
          read(102,10)
          read(103,10)
          read(104,10)
          read(105,10)
        enddo
        
        
        close(100)
        close(101)
        close(102)
        close(103)
        close(104)
        close(105)
    
        10 FORMAT (E16.8)  
        11 FORMAT (E16.8,E16.8) 
        
        !change the directory back
        istat=chdir('..')
        if (istat .ne. 0) then
            write(logInput,*) "Directory did not changed!!!"
            call writeLog()
            call exit()
				end if
				
				!make sure the exciton wavefunctions are normalized
				do iX=1,currcnt%nX
					do iKcm=currcnt%iKcm_min,currcnt%iKcm_max
						tmpr1 = 0.d0
						tmpr2 = 0.d0
						tmpr3 = 0.d0
						do ikr=currcnt%ikr_low,currcnt%ikr_high
              tmpr1 = tmpr1 + dble(conjg(currcnt%Psi_A1(ikr,iX,iKcm)*currcnt%Psi_A1(ikr,iX,iKcm)))
							tmpr2 = tmpr2 + dble(conjg(currcnt%Psi0_A2(ikr,iX,iKcm)*currcnt%Psi0_A2(ikr,iX,iKcm)))
							tmpr3 = tmpr3 + dble(conjg(currcnt%Psi1_A2(ikr,iX,iKcm)*currcnt%Psi1_A2(ikr,iX,iKcm)))
						enddo
						currcnt%Psi_A1(:,iX,iKcm) = currcnt%Psi_A1(:,iX,iKcm) / dcmplx(tmpr1)
						currcnt%Psi0_A2(:,iX,iKcm) = currcnt%Psi0_A2(:,iX,iKcm) / dcmplx(tmpr2)
						currcnt%Psi1_A2(:,iX,iKcm) = currcnt%Psi1_A2(:,iX,iKcm) / dcmplx(tmpr3)
          enddo
				enddo
        
      end subroutine loadExcitonWavefunction
      
      !**************************************************************************************************************************
      ! load the data from miscellaneous.dat file from ExcitonEnergy calculation
      !**************************************************************************************************************************
      subroutine loadMiscData(currCNT)
        use inputParameters
        use smallFunctions
        type(cnt), intent(inout) :: currCNT
        integer(4) :: istat
        integer :: min_sub, ik_max, iKcm_max, ikr_high, ik_high, iq_max, nX
        real*8 :: dk
        
        !change the directory to that of the CNT
        istat=chdir(currcnt%excitonDirectory)
        if (istat .ne. 0) then
            write(logInput,*) "Directory did not changed!!!"
            call writeLog()
            call exit()
				end if
				
        open(unit=100,file='miscellaneous.dat',status="old")
        read(100,10) min_sub
        read(100,10) ik_max
        read(100,10) iKcm_max
        read(100,10) ikr_high
        read(100,10) ik_high
        read(100,10) iq_max
        read(100,10) nX
        read(100,11) dk
        
        close(100)
  
        if (min_sub .ne. currcnt%min_sub(currcnt%i_sub)) then
          write(logInput,*) min_sub
          call writeLog()
          write(logInput,*) currcnt%min_sub(currcnt%i_sub)
          call writeLog()
          write(logInput,*) "Error in matching min_sub variables!"
          call writeLog()
          call exit()
        end if

        if (ik_max .ne. currcnt%ik_max) then
          write(logInput,*) "Error in matching ik_max variables!"
          call writeLog()
          call exit()
        end if
        
        if (iKcm_max .ne. currcnt%iKcm_max) then
          write(logInput,*) "Error in matching iKcm_max variables!"
          call writeLog()
          call exit()
        end if

        if (ikr_high .ne. currcnt%ikr_high) then
          write(logInput,*) "Error in matching ikr_high variables!"
          call writeLog()
          call exit()
        end if
        
        if (ik_high .ne. currcnt%ik_high) then
          write(logInput,*) "Error in matching ik_high variables!"
          call writeLog()
          call exit()
        end if
        
        if (iq_max .ne. currcnt%iq_max) then
          write(logInput,*) "Error in matching iq_max variables!"
          call writeLog()
          call exit()
        end if

        ! you don't need to check dk since there are round off errors in real numbers when writing and reading from files
        currcnt%nX = nX
    
        10 FORMAT (I4.4) 
        11 FORMAT (E16.8)
        
        !change the directory back
        istat=chdir('..')
        if (istat .ne. 0) then
            write(logInput,*) "Directory did not changed!!!"
            call writeLog()
            call exit()
        end if
        
      end subroutine loadMiscData
    
      
      !**************************************************************************************************************************
      ! save CNT dispersions and the crossing points and the same energy points
      !**************************************************************************************************************************
      subroutine saveTransitionPoints(cnt1,cnt2)
        use prepareForster_module
				use inputParameters
        use smallFunctions
        type(cnt), intent(in) :: cnt1,cnt2
        character(len=200) :: command
        integer(4) :: istat
        integer :: iKcm, iX, i
        
        !create and change the directory to that of the CNT
        write(command,'("mkdir ",A100)') outputDirectory
        call execute_command_line(command,wait=.true.,exitstat=i)

        istat=chdir(outputDirectory)
        if (istat .ne. 0) then
            write(logInput,*) "Directory did not changed!!!"
            call writeLog()
            call exit()
        end if
        
        !write carbon nanotube 1 k_vector
        open(unit=100,file='cnt1_kvec.dat',status="unknown")
        do iKcm=cnt1%iKcm_min,cnt1%iKcm_max
          write(100,10, advance='no') dble(iKcm)*cnt1%dk
        enddo
        close(100)
        
        !write carbon nanotube 2 k_vector
        open(unit=100,file='cnt2_kvec.dat',status="unknown")
        do iKcm=cnt2%iKcm_min,cnt2%iKcm_max
          write(100,10, advance='no') dble(iKcm)*cnt2%dk
        enddo
        close(100)
        
        !write carbon nanotube 1 Ex0_A2 dispersion
        open(unit=100,file='cnt1_Ex0_A2.dat',status="unknown")
        do iKcm=cnt1%iKcm_min,cnt1%iKcm_max
          do iX=1,cnt1%nX
            write(100,10, advance='no') cnt1%Ex0_A2(iX,iKcm)
          enddo
          write(100,10)
        enddo
        close(100)
    
        !write carbon nanotube 2 Ex0_A2 dispersion
        open(unit=100,file='cnt2_Ex0_A2.dat',status="unknown")
        do iKcm=cnt2%iKcm_min,cnt2%iKcm_max
          do iX=1,cnt2%nX
            write(100,10, advance='no') cnt2%Ex0_A2(iX,iKcm)
          enddo
          write(100,10)
        enddo
        close(100) 
        
        !write crossing points indexes
        open(unit=100,file='crossingPoints.dat',status="unknown")
        do i=1,ubound(crossingPoints,1)
          write(100,11) crossingPoints(i,1), crossingPoints(i,2), crossingPoints(i,3)
        enddo
        close(100) 
        
				!write crossing points indexes
        open(unit=100,file='sameEnergy.dat',status="unknown")
        do i=1,ubound(sameEnergy,1)
          write(100,12) sameEnergy(i,1), sameEnergy(i,2), sameEnergy(i,3), sameEnergy(i,4)
        enddo
        close(100) 
				
10			FORMAT (E16.8)
11			FORMAT (4I8, 4I8, 4I8)
12			FORMAT (4I8, 4I8, 4I8, 4I8)
					 
  
       
        
        !change the directory back
        istat=chdir('..')
        if (istat .ne. 0) then
            write(logInput,*) "Directory did not changed!!!"
            call writeLog()
            call exit()
        end if
        
			end subroutine saveTransitionPoints
			
			!**************************************************************************************************************************
      ! save total exciton density of states for a given cnt
      !**************************************************************************************************************************
      subroutine saveDOS(cnt1, cnt2)
        use prepareForster_module
				use inputParameters
        use smallFunctions
        type(cnt), intent(in) :: cnt1, cnt2
        character(len=200) :: command
        integer(4) :: istat
				real*8 :: dos
        
        integer :: iKcm, iX, i
        
        !create and change the directory to that of the CNT
        write(command,'("mkdir ",A100)') outputDirectory
        call execute_command_line(command,wait=.true.,exitstat=i)

        istat=chdir(outputDirectory)
        if (istat .ne. 0) then
            write(logInput,*) "Directory did not changed!!!"
            call writeLog()
            call exit()
        end if
        
        !write cnt1 Ex0_A2 dispersion
        open(unit=100,file='cnt1_DOS.dat',status="unknown")
        do iKcm=cnt1%iKcm_min,cnt1%iKcm_max
          do iX=1,cnt1%nX
						call calculateDOS(cnt1,iKcm,iX,DOS)
            write(100,10, advance='no') dos
          enddo
          write(100,10)
        enddo
        close(100)
        
				!write cnt2 Ex0_A2 dispersion
        open(unit=100,file='cnt2_DOS.dat',status="unknown")
        do iKcm=cnt2%iKcm_min,cnt2%iKcm_max
          do iX=1,cnt2%nX
						call calculateDOS(cnt2,iKcm,iX,DOS)
            write(100,10, advance='no') dos
          enddo
          write(100,10)
        enddo
        close(100)
				
        10 FORMAT (E16.8)
        
        !change the directory back
        istat=chdir('..')
        if (istat .ne. 0) then
            write(logInput,*) "Directory did not changed!!!"
            call writeLog()
            call exit()
        end if
        
			end subroutine saveDOS
end module dataClass