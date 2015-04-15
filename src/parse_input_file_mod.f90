!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Declaration of input parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
module parse_input_file_mod
	implicit none
	private
	public :: parse_input_file

contains	
	!**************************************************************************************************************************
	! parse the input file
	!**************************************************************************************************************************

	subroutine parse_input_file()
		use comparams, only: cnt1, cnt2
		use physicalConstants, only: eV
		use output_module, only: writeLog

		character(len=100) :: filename
		character(len=200) :: buffer, cnt_name, label, value
		character(len=200) :: outdir	
		character(len=200) :: indir
		integer :: istat=0
		integer :: ios=0
		integer :: pos_comma=0, pos_equal=0

		if (command_argument_count() .ne. 1) then
			write(*,*) "Input format ERROR!"
			write(*,*) "Correct input format is: main.exe inputfile.in"
			call exit()
		end if


		call get_command_argument(1,filename)
		open(unit=100,file=filename,status="old", action="read", iostat=istat)
		if (istat .ne. 0) then
			write(*,*) ""
			write(*,*) "Unable to read input file:", filename
			call exit()
		end if

		read (100,'(A)',iostat=ios) buffer
		indir = trim(adjustl(buffer))

		do while (ios == 0)
			read (100,'(A)',iostat=ios) buffer
			if (ios == 0) then
				if (buffer(1:1) .ne. '#') then
					pos_comma = scan(buffer,',')
					pos_equal = scan(buffer,'=')
					cnt_name = adjustl(buffer(1:pos_comma-1))
					label = adjustl(buffer(pos_comma+1:pos_equal-1))
					value = adjustl(buffer(pos_equal+1:))

					! set the target cnt
					select case (trim(cnt_name))
					case ('cnt1')
						select case (trim(label))
						case ('n_ch')
							read(value, *) cnt1%n_ch
						case ('m_ch')
							read(value, *) cnt1%m_ch
						case ('nkg')
							read(value, *) cnt1%nkg
						case ('nr')
							read(value, *) cnt1%nr
						case ('E_th[eV]')
							read(value, *) cnt1%E_th
							cnt1%E_th = cnt1%E_th * eV
						case ('Kcm_max[1/nm]')
							read(value, *) cnt1%Kcm_max
							cnt1%Kcm_max = cnt1%Kcm_max * 1.d9
						case ('i_sub')
							read(value, *) cnt1%i_sub
						case ('Ckappa')
							read(value, *) cnt1%Ckappa
						case ('target_exciton_type')
							read(value, *) cnt1%targetExcitonType
						case ('length[1/nm]')
							read(value, *) cnt1%length
							cnt1%length = cnt1%length * 1.d-9
						case ('center_position[nm]')
							read(value, *) cnt1%center_position
							cnt1%center_position = cnt1%center_position * 1.d-9
						end select
					case ('cnt2')
						select case (trim(label))
						case ('n_ch')
							read(value, *) cnt2%n_ch
						case ('m_ch')
							read(value, *) cnt2%m_ch
						case ('nkg')
							read(value, *) cnt2%nkg
						case ('nr')
							read(value, *) cnt2%nr
						case ('E_th[eV]')
							read(value, *) cnt2%E_th
							cnt2%E_th = cnt2%E_th * eV
						case ('Kcm_max[1/nm]')
							read(value, *) cnt2%Kcm_max
							cnt2%Kcm_max = cnt2%Kcm_max * 1.d9
						case ('i_sub')
							read(value, *) cnt2%i_sub
						case ('Ckappa')
							read(value, *) cnt2%Ckappa
						case ('target_exciton_type')
							read(value, *) cnt2%targetExcitonType
						case ('length[1/nm]')
							read(value, *) cnt2%length
							cnt2%length = cnt2%length * 1.d-9
						case ('center_position[nm]')
							read(value, *) cnt2%center_position
							cnt2%center_position = cnt2%center_position * 1.d-9
						end select
					end select
				end if
			else if (ios .gt. 0) then
				write (*,*) "Error in reading input file!"
				call exit()
			end if
		end do
		close(100)

		write(cnt1%directory,"( A, 'CNT(', I2.2, ',', I2.2, ')-nkg(', I4.4, ')-nr(', I4.4, ')-E_th(', F3.1, ')-Kcm_max(', F3.1, ')-i_sub(', I1.1, ')-Ckappa(', F3.1, ')/' )") trim(indir), cnt1%n_ch, cnt1%m_ch, cnt1%nkg, cnt1%nr, cnt1%E_th/eV, cnt1%Kcm_max*1.d-9, cnt1%i_sub, cnt1%Ckappa
		write(cnt2%directory,"( A, 'CNT(', I2.2, ',', I2.2, ')-nkg(', I4.4, ')-nr(', I4.4, ')-E_th(', F3.1, ')-Kcm_max(', F3.1, ')-i_sub(', I1.1, ')-Ckappa(', F3.1, ')/' )") trim(indir), cnt2%n_ch, cnt2%m_ch, cnt2%nkg, cnt2%nr, cnt2%E_th/eV, cnt2%Kcm_max*1.d-9, cnt2%i_sub, cnt2%Ckappa
		write(outdir,"('Transfer-',A,'(',I2.2,',',I2.2,')-to-',A,'(',I2.2,',',I2.2,')-Ckappa(',F3.1,')')") trim(cnt1%targetExcitonType), cnt1%n_ch, cnt1%m_ch, trim(cnt2%targetExcitonType), cnt2%n_ch, cnt2%m_ch, cnt1%Ckappa

		call create_outdir(outdir)

		call writeLog(cnt1%directory)
		call writeLog(cnt2%directory)
		call writeLog(outdir)

		return
	end subroutine parse_input_file

	!*******************************************************************************
	! This subroutines changes the working directory to the output directory for saving files
	!*******************************************************************************
	
	subroutine create_outdir(outdir)
		use output_module, only: writeLog

		character(len=200), intent(in) :: outdir
		integer :: istat=0
		character(len=300) :: command
		integer, dimension(3) :: date, time
		character(len=100) :: logInput

		! specifiy the output directory
		write(command,'("rm -rf ''",A,"''")') trim(outdir) !remove the directory if it already exists
		call system(command)
		write(command,'("mkdir ''",A,"''")') trim(outdir) !create the directory again
		call system(trim(command))

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
		call writeLog(logInput)

		return	
	end subroutine create_outdir

	
end module parse_input_file_mod