!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Declaration of input parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
module initial_final_mod
	implicit none
	private
	public :: initialize, finalize

contains	
	!**************************************************************************************************************************
	! read the CNT information from the log file of exciton calculation which is done previously
	!**************************************************************************************************************************

	subroutine initialize(cnt1,cnt2)
		use comparams, only: starttime
		use cnt_class, only: cnt
		use input_class, only: inputCNT
		use output_module, only: open_output_directory, writeLog, logInput

		type (cnt), intent(inout) :: cnt1, cnt2
		character(len=100) :: buffer
		
		call CPU_time(starttime)

		if (command_argument_count() .ne. 4) then
			write(*,*) "ERROR! Input format: main.exe cnt1_folder cnt2_folder"
			call exit()
		end if 

		! read the information of the first CNT
		call get_command_argument(1,buffer)
		buffer = trim(buffer)
		call inputCNT(cnt1,buffer)
		! set the target exciton type
		call get_command_argument(2,buffer)
		select case (trim(buffer))
		case ('Ex_A1')
			cnt1%Ex_t = cnt1%Ex_A1
			cnt1%Psi_t = cnt1%Psi_A1
			write(logInput,*) 'cnt1 target exciton = Ex_A1'
			call writeLog()
		case ('Ex0_A2')
			cnt1%Ex_t = cnt1%Ex0_A2
			cnt1%Psi_t = cnt1%Psi0_A2
			write(logInput,*) 'cnt1 target exciton = Ex0_A2'
			call writeLog()
		case ('Ex1_A2')
			cnt1%Ex_t = cnt1%Ex1_A2
			cnt1%Psi_t = cnt1%Psi1_A2
			write(logInput,*) 'cnt1 target exciton = Ex1_A2'
			call writeLog()
		case default
			write(*,*) "Could not recognize target exciton type!!!!"
			call exit()
		end select

		! read the information of the second CNT
		call get_command_argument(3,buffer)
		buffer = trim(buffer)
		call inputCNT(cnt2,buffer)
		! set the target exciton type
		call get_command_argument(4,buffer)
		select case (trim(buffer))
		case ('Ex_A1')
			cnt2%Ex_t = cnt2%Ex_A1
			cnt2%Psi_t = cnt2%Psi_A1
			write(logInput,*) 'cnt2 target exciton = Ex_A1'
			call writeLog()
		case ('Ex0_A2')
			cnt2%Ex_t = cnt2%Ex0_A2
			cnt2%Psi_t = cnt2%Psi0_A2
			write(logInput,*) 'cnt2 target exciton = Ex0_A2'
			call writeLog()
		case ('Ex1_A2')
			cnt2%Ex_t = cnt2%Ex1_A2
			cnt2%Psi_t = cnt2%Psi1_A2
			write(logInput,*) 'cnt2 target exciton = Ex1_A2'
			call writeLog()
		case default
			write(*,*) "Could not recognize target exciton type!!!!"
			call exit()
		end select

		!if Ckappa is not equal for both cnt1 and cnt2 print a warning and exit
		if (cnt1%Ckappa .ne. cnt2%Ckappa) then
			write(*,*) "cnt1 and cnt2 have different Ckappa"
			call exit()
		end if

		call open_output_directory(cnt1, cnt2)

		return
	end subroutine initialize

	subroutine finalize()
		use comparams, only: endtime
		use output_module, only: close_output_directory

		call CPU_time(endtime)
		call close_output_directory()

		return
	end subroutine finalize

	
end module initial_final_mod