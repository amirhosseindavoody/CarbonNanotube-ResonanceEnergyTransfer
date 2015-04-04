!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Declaration of input parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
module output_module
	implicit none
	private
	public :: open_output_directory, close_output_directory, writeLog

	character(len=200), public :: logInput

contains
	!*******************************************************************************
	! This subroutines opens the log file add new log and closes the file
	!*******************************************************************************
	
	subroutine writeLog()
		logical :: flgexist
		integer :: logFile = 10

		inquire(file="log.dat",exist=flgexist)
		if (flgexist) then
			open(logFile, file="log.dat", status="old", position="append", action="write")
		else
			open(logFile, file="log.dat", status="new", action="write")
		end if
		write(logFile, *) logInput
		close(logFile)

		return
	end subroutine writeLog

	!*******************************************************************************
	! This subroutines changes the working directory to the output directory for saving files
	!*******************************************************************************
	
	subroutine open_output_directory(cnt1, cnt2)
		use cnt_class, only: cnt

		type(cnt), intent(in) :: cnt1, cnt2
		integer :: istat=0
		integer :: i
		character(len=200) :: command
		integer, dimension(3) :: date, time
		character(len=100) :: outdir !this is the output directory for storing output results

		write(outdir,"('Transfer-',A,'(',I2.2,',',I2.2,')-to-',A,'(',I2.2,',',I2.2,')-Ckappa(',F3.1,')')") trim(cnt1%targetExcitonType), cnt1%n_ch, cnt1%m_ch, trim(cnt2%targetExcitonType), cnt2%n_ch, cnt2%m_ch, cnt1%Ckappa

		! specifiy the output directory
		write(command,'("rm -rf ''",A,"''")') trim(outdir) !remove the directory if it already exists
		call execute_command_line(command,wait=.true.,exitstat=i)
    	write(command,'("mkdir ''",A,"''")') trim(outdir) !create the directory again
    	call execute_command_line(command,wait=.true.,exitstat=i)
!     	write(command,'("mv log.dat ''",A,"''")') trim(outdir) !move the log file to the new output directory
!     	call execute_command_line(command,wait=.true.,exitstat=i)
    	
    	istat=chdir(trim(outdir))
		if (istat .ne. 0) then
			write(*,*) "Directory did not changed!!!"
			write(*,*) "Simulation stopped!!!"
			call exit()
		end if

		! get time and date of start of simulation
		call idate(date)
		call itime(time)

		! write simulation inputs to the log file
		write(logInput,'("Simulation started at--> DATE=",I2.2,"/",I2.2,"/",I4.4,"  TIME=",I2.2,":",I2.2,":",I2.2)') date, time
		call writeLog()

		return	
	end subroutine open_output_directory

	!*******************************************************************************
	! This subroutines changes the working directory back to root directory
	!*******************************************************************************
	subroutine close_output_directory()
		use comparams, only: starttime, endtime
		
		integer :: istat
		write(logInput,'("Run time = ",f10.3," seconds.")'),endtime-starttime
		call writeLog()

		!change the directory back
		istat=chdir('..')
		if (istat .ne. 0) then
			write(logInput,*) "Directory did not changed!!!"
			call writeLog()
			call exit()
		end if

		return	
	end subroutine close_output_directory
	
end module output_module