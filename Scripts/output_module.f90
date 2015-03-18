!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Declaration of input parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
module output_module
	use comparams
	use smallFunctions
	implicit none

contains
	subroutine open_output_directory()
		integer :: istat=0
		integer :: i
		character(len=200) :: command
		integer, dimension(3) :: date, time
		character(len=100) :: directory

		! specifiy the output directory
		write(command,'("rm -rf ''",A100,"''")') outdir !remove the directory if it already exists
		call execute_command_line(command,wait=.true.,exitstat=i)
    	write(command,'("mkdir ''",A100,"''")') outdir !create the directory again
    	call execute_command_line(command,wait=.true.,exitstat=i)
    	print *, outdir
    	write(directory,'(A)') outdir
    	istat=chdir(directory)
		if (istat .ne. 0) then
			write(logInput,*) "Directory did not changed!!!"
			call writeLog()
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

	subroutine close_output_directory()
		integer :: istat
		
		call CPU_time(endtime)
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