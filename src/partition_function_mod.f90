module partition_function_mod
	implicit none
	private
    public  :: calculate_partition_function, calculate_partition_function_table
    
    real*8, dimension(:,:), allocatable, public :: partition_function_table	

contains
	!**************************************************************************************************************************
	! calculate the partition function for target or all exciton types of a given carbon nanotube
	!**************************************************************************************************************************
	
	subroutine calculate_partition_function_table()
! 		use cnt_class, only: cnt
		use comparams, only: max_temperature, min_temperature, temperature_steps, cnt1, cnt2
		use physicalConstants, only : kb

		integer :: iTemperature
		integer :: iKcm, ix
		real*8 :: temperature
		real*8 :: partition_function

		!allocate the partition_function_table
		if (.not. allocated(partition_function_table)) then
			allocate(partition_function_table(2,temperature_steps))
			partition_function_table = 0.d0
		endif

		! read the information of A1 exciton
		if (allocated(cnt1%Ex_t)) deallocate(cnt1%Ex_t)
		allocate(cnt1%Ex_t(1:cnt1%nX_a,cnt1%iKcm_min_fine:cnt1%iKcm_max_fine))

		cnt1%nX_t=cnt1%nX_a
		
		open(unit=100,file=trim(cnt1%directory)//'Ex_A1.dat',status="old")
		do iKcm=cnt1%iKcm_min_fine,cnt1%iKcm_max_fine
			do iX=1,cnt1%nX_a
				read(100,'(E16.8)', advance='no') cnt1%Ex_t(iX,iKcm)
			enddo
			read(100,'(E16.8)')
		enddo
		close(100)

		do iTemperature = 1,temperature_steps
			temperature = min_temperature + dble(iTemperature-1) * (max_temperature-min_temperature) / dble(temperature_steps-1)
			call calculate_partition_function(cnt1, temperature, partition_function)
			! A1 excitons have both singlet and triplet states which have the same energy, thus the factor of 2.0 in the partition function.
			partition_function_table(1,iTemperature) = partition_function_table(1,iTemperature) + 2.d0 * partition_function
		enddo

		deallocate(cnt1%Ex_t)

		! read the information of singlet A2 exciton
		if (allocated(cnt1%Ex_t)) deallocate(cnt1%Ex_t)
		allocate(cnt1%Ex_t(1:cnt1%nX_a,cnt1%iKcm_min_fine:cnt1%iKcm_max_fine))
		
		cnt1%nX_t=cnt1%nX_a

		open(unit=100,file=trim(cnt1%directory)//'Ex0_A2.dat',status="old")
		do iKcm=cnt1%iKcm_min_fine,cnt1%iKcm_max_fine
			do iX=1,cnt1%nX_a
				read(100,'(E16.8)', advance='no') cnt1%Ex_t(iX,iKcm)
			enddo
			read(100,'(E16.8)')
		enddo
		close(100)

		do iTemperature = 1,temperature_steps
			temperature = min_temperature + dble(iTemperature-1) * (max_temperature-min_temperature) / dble(temperature_steps-1)
			call calculate_partition_function(cnt1, temperature, partition_function)
			partition_function_table(1,iTemperature) = partition_function_table(1,iTemperature) + partition_function
		enddo

		deallocate(cnt1%Ex_t)

		! read the information of triplet A2 exciton
		if (allocated(cnt1%Ex_t)) deallocate(cnt1%Ex_t)
		allocate(cnt1%Ex_t(1:cnt1%nX_a,cnt1%iKcm_min_fine:cnt1%iKcm_max_fine))
		
		cnt1%nX_t=cnt1%nX_a

		open(unit=100,file=trim(cnt1%directory)//'Ex1_A2.dat',status="old")
		do iKcm=cnt1%iKcm_min_fine,cnt1%iKcm_max_fine
			do iX=1,cnt1%nX_a
				read(100,'(E16.8)', advance='no') cnt1%Ex_t(iX,iKcm)
			enddo
			read(100,'(E16.8)')
		enddo
		close(100)

		do iTemperature = 1,temperature_steps
			temperature = min_temperature + dble(iTemperature-1) * (max_temperature-min_temperature) / dble(temperature_steps-1)
			call calculate_partition_function(cnt1, temperature, partition_function)
			partition_function_table(1,iTemperature) = partition_function_table(1,iTemperature) + partition_function
		enddo

		deallocate(cnt1%Ex_t)

		! read the information of singlet Ep exciton
		if (allocated(cnt1%Ex_t)) deallocate(cnt1%Ex_t)
		allocate(cnt1%Ex_t(1:cnt1%nX_a,cnt1%iKcm_min_fine:cnt1%iKcm_max_fine))
		
		cnt1%nX_t=cnt1%nX_e

		open(unit=100,file=trim(cnt1%directory)//'Ex0_Ep.dat',status="old")
		do iKcm=cnt1%iKcm_min_fine,cnt1%iKcm_max_fine
			do iX=1,cnt1%nX_e
				read(100,'(E16.8)', advance='no') cnt1%Ex_t(iX,iKcm)
			enddo
			read(100,'(E16.8)')
		enddo
		close(100)

		do iTemperature = 1,temperature_steps
			temperature = min_temperature + dble(iTemperature-1) * (max_temperature-min_temperature) / dble(temperature_steps-1)
			call calculate_partition_function(cnt1, temperature, partition_function)
			partition_function_table(1,iTemperature) = partition_function_table(1,iTemperature) + partition_function
		enddo

		deallocate(cnt1%Ex_t)

		! read the information of singlet Em exciton
		if (allocated(cnt1%Ex_t)) deallocate(cnt1%Ex_t)
		allocate(cnt1%Ex_t(1:cnt1%nX_a,cnt1%iKcm_min_fine:cnt1%iKcm_max_fine))
		
		cnt1%nX_t=cnt1%nX_e

		open(unit=100,file=trim(cnt1%directory)//'Ex0_Em.dat',status="old")
		do iKcm=cnt1%iKcm_min_fine,cnt1%iKcm_max_fine
			do iX=1,cnt1%nX_e
				read(100,'(E16.8)', advance='no') cnt1%Ex_t(iX,iKcm)
			enddo
			read(100,'(E16.8)')
		enddo
		close(100)

		do iTemperature = 1,temperature_steps
			temperature = min_temperature + dble(iTemperature-1) * (max_temperature-min_temperature) / dble(temperature_steps-1)
			call calculate_partition_function(cnt1, temperature, partition_function)
			partition_function_table(1,iTemperature) = partition_function_table(1,iTemperature) + partition_function
		enddo

		deallocate(cnt1%Ex_t)

		! read the information of triplet Ep exciton
		if (allocated(cnt1%Ex_t)) deallocate(cnt1%Ex_t)
		allocate(cnt1%Ex_t(1:cnt1%nX_a,cnt1%iKcm_min_fine:cnt1%iKcm_max_fine))
		
		cnt1%nX_t=cnt1%nX_e

		open(unit=100,file=trim(cnt1%directory)//'Ex1_Ep.dat',status="old")
		do iKcm=cnt1%iKcm_min_fine,cnt1%iKcm_max_fine
			do iX=1,cnt1%nX_e
				read(100,'(E16.8)', advance='no') cnt1%Ex_t(iX,iKcm)
			enddo
			read(100,'(E16.8)')
		enddo
		close(100)

		do iTemperature = 1,temperature_steps
			temperature = min_temperature + dble(iTemperature-1) * (max_temperature-min_temperature) / dble(temperature_steps-1)
			call calculate_partition_function(cnt1, temperature, partition_function)
			partition_function_table(1,iTemperature) = partition_function_table(1,iTemperature) + partition_function
		enddo

		deallocate(cnt1%Ex_t)
        
		! read the information of triplet Em exciton
		if (allocated(cnt1%Ex_t)) deallocate(cnt1%Ex_t)
		allocate(cnt1%Ex_t(1:cnt1%nX_a,cnt1%iKcm_min_fine:cnt1%iKcm_max_fine))
		
		cnt1%nX_t=cnt1%nX_e

		open(unit=100,file=trim(cnt1%directory)//'Ex1_Em.dat',status="old")
		do iKcm=cnt1%iKcm_min_fine,cnt1%iKcm_max_fine
			do iX=1,cnt1%nX_e
				read(100,'(E16.8)', advance='no') cnt1%Ex_t(iX,iKcm)
			enddo
			read(100,'(E16.8)')
		enddo
		close(100)

		do iTemperature = 1,temperature_steps
			temperature = min_temperature + dble(iTemperature-1) * (max_temperature-min_temperature) / dble(temperature_steps-1)
			call calculate_partition_function(cnt1, temperature, partition_function)
			partition_function_table(1,iTemperature) = partition_function_table(1,iTemperature) + partition_function
		enddo

		deallocate(cnt1%Ex_t)

		! read the information of A1 exciton
		if (allocated(cnt2%Ex_t)) deallocate(cnt2%Ex_t)
		allocate(cnt2%Ex_t(1:cnt2%nX_a,cnt2%iKcm_min_fine:cnt2%iKcm_max_fine))

		cnt2%nX_t=cnt2%nX_a
		
		open(unit=100,file=trim(cnt2%directory)//'Ex_A1.dat',status="old")
		do iKcm=cnt2%iKcm_min_fine,cnt2%iKcm_max_fine
			do iX=1,cnt2%nX_a
				read(100,'(E16.8)', advance='no') cnt2%Ex_t(iX,iKcm)
			enddo
			read(100,'(E16.8)')
		enddo
		close(100)

		do iTemperature = 1,temperature_steps
			temperature = min_temperature + dble(iTemperature-1) * (max_temperature-min_temperature) / dble(temperature_steps-1)
			call calculate_partition_function(cnt2, temperature, partition_function)
			! A1 excitons have both singlet and triplet states which have the same energy, thus the factor of 2.0 in the partition function.
			partition_function_table(2,iTemperature) = partition_function_table(2,iTemperature) + 2.d0 * partition_function
		enddo

		deallocate(cnt2%Ex_t)

		! read the information of singlet A2 exciton
		if (allocated(cnt2%Ex_t)) deallocate(cnt2%Ex_t)
		allocate(cnt2%Ex_t(1:cnt2%nX_a,cnt2%iKcm_min_fine:cnt2%iKcm_max_fine))
		
		cnt2%nX_t=cnt2%nX_a

		open(unit=100,file=trim(cnt2%directory)//'Ex0_A2.dat',status="old")
		do iKcm=cnt2%iKcm_min_fine,cnt2%iKcm_max_fine
			do iX=1,cnt2%nX_a
				read(100,'(E16.8)', advance='no') cnt2%Ex_t(iX,iKcm)
			enddo
			read(100,'(E16.8)')
		enddo
		close(100)

		do iTemperature = 1,temperature_steps
			temperature = min_temperature + dble(iTemperature-1) * (max_temperature-min_temperature) / dble(temperature_steps-1)
			call calculate_partition_function(cnt2, temperature, partition_function)
			partition_function_table(2,iTemperature) = partition_function_table(2,iTemperature) + partition_function
		enddo

		deallocate(cnt2%Ex_t)

		! read the information of triplet A2 exciton
		if (allocated(cnt2%Ex_t)) deallocate(cnt2%Ex_t)
		allocate(cnt2%Ex_t(1:cnt2%nX_a,cnt2%iKcm_min_fine:cnt2%iKcm_max_fine))
		
		cnt2%nX_t=cnt2%nX_a

		open(unit=100,file=trim(cnt2%directory)//'Ex1_A2.dat',status="old")
		do iKcm=cnt2%iKcm_min_fine,cnt2%iKcm_max_fine
			do iX=1,cnt2%nX_a
				read(100,'(E16.8)', advance='no') cnt2%Ex_t(iX,iKcm)
			enddo
			read(100,'(E16.8)')
		enddo
		close(100)

		do iTemperature = 1,temperature_steps
			temperature = min_temperature + dble(iTemperature-1) * (max_temperature-min_temperature) / dble(temperature_steps-1)
			call calculate_partition_function(cnt2, temperature, partition_function)
			partition_function_table(2,iTemperature) = partition_function_table(2,iTemperature) + partition_function
		enddo

		deallocate(cnt2%Ex_t)

		! read the information of singlet Ep exciton
		if (allocated(cnt2%Ex_t)) deallocate(cnt2%Ex_t)
		allocate(cnt2%Ex_t(1:cnt2%nX_a,cnt2%iKcm_min_fine:cnt2%iKcm_max_fine))
		
		cnt2%nX_t=cnt2%nX_e

		open(unit=100,file=trim(cnt2%directory)//'Ex0_Ep.dat',status="old")
		do iKcm=cnt2%iKcm_min_fine,cnt2%iKcm_max_fine
			do iX=1,cnt2%nX_e
				read(100,'(E16.8)', advance='no') cnt2%Ex_t(iX,iKcm)
			enddo
			read(100,'(E16.8)')
		enddo
		close(100)

		do iTemperature = 1,temperature_steps
			temperature = min_temperature + dble(iTemperature-1) * (max_temperature-min_temperature) / dble(temperature_steps-1)
			call calculate_partition_function(cnt2, temperature, partition_function)
			partition_function_table(2,iTemperature) = partition_function_table(2,iTemperature) + partition_function
		enddo

		deallocate(cnt2%Ex_t)

		! read the information of singlet Em exciton
		if (allocated(cnt2%Ex_t)) deallocate(cnt2%Ex_t)
		allocate(cnt2%Ex_t(1:cnt2%nX_a,cnt2%iKcm_min_fine:cnt2%iKcm_max_fine))
		
		cnt2%nX_t=cnt2%nX_e

		open(unit=100,file=trim(cnt2%directory)//'Ex0_Em.dat',status="old")
		do iKcm=cnt2%iKcm_min_fine,cnt2%iKcm_max_fine
			do iX=1,cnt2%nX_e
				read(100,'(E16.8)', advance='no') cnt2%Ex_t(iX,iKcm)
			enddo
			read(100,'(E16.8)')
		enddo
		close(100)

		do iTemperature = 1,temperature_steps
			temperature = min_temperature + dble(iTemperature-1) * (max_temperature-min_temperature) / dble(temperature_steps-1)
			call calculate_partition_function(cnt2, temperature, partition_function)
			partition_function_table(2,iTemperature) = partition_function_table(2,iTemperature) + partition_function
		enddo

		deallocate(cnt2%Ex_t)

		! read the information of triplet Ep exciton
		if (allocated(cnt2%Ex_t)) deallocate(cnt2%Ex_t)
		allocate(cnt2%Ex_t(1:cnt2%nX_a,cnt2%iKcm_min_fine:cnt2%iKcm_max_fine))
		
		cnt2%nX_t=cnt2%nX_e

		open(unit=100,file=trim(cnt2%directory)//'Ex1_Ep.dat',status="old")
		do iKcm=cnt2%iKcm_min_fine,cnt2%iKcm_max_fine
			do iX=1,cnt2%nX_e
				read(100,'(E16.8)', advance='no') cnt2%Ex_t(iX,iKcm)
			enddo
			read(100,'(E16.8)')
		enddo
		close(100)

		do iTemperature = 1,temperature_steps
			temperature = min_temperature + dble(iTemperature-1) * (max_temperature-min_temperature) / dble(temperature_steps-1)
			call calculate_partition_function(cnt2, temperature, partition_function)
			partition_function_table(2,iTemperature) = partition_function_table(2,iTemperature) + partition_function
		enddo

		deallocate(cnt2%Ex_t)
        
		! read the information of triplet Em exciton
		if (allocated(cnt2%Ex_t)) deallocate(cnt2%Ex_t)
		allocate(cnt2%Ex_t(1:cnt2%nX_a,cnt2%iKcm_min_fine:cnt2%iKcm_max_fine))
		
		cnt2%nX_t=cnt2%nX_e
		
		open(unit=100,file=trim(cnt2%directory)//'Ex1_Em.dat',status="old")
		do iKcm=cnt2%iKcm_min_fine,cnt2%iKcm_max_fine
			do iX=1,cnt2%nX_e
				read(100,'(E16.8)', advance='no') cnt2%Ex_t(iX,iKcm)
			enddo
			read(100,'(E16.8)')
		enddo
		close(100)

		do iTemperature = 1,temperature_steps
			temperature = min_temperature + dble(iTemperature-1) * (max_temperature-min_temperature) / dble(temperature_steps-1)
			call calculate_partition_function(cnt2, temperature, partition_function)
			partition_function_table(2,iTemperature) = partition_function_table(2,iTemperature) + partition_function
		enddo

		deallocate(cnt2%Ex_t)

		return        
	end subroutine calculate_partition_function_table

	!**************************************************************************************************************************
	! calculate the partition function for target or all exciton types of a given carbon nanotube
	!**************************************************************************************************************************
	
	subroutine calculate_partition_function(currcnt, temperature, partition_function)
		use cnt_class, only: cnt
		use comparams, only: max_temperature
		use physicalConstants, only : kb

		type(cnt), intent(in) :: currcnt
		real*8, intent(in) :: temperature
		real*8, intent(out) :: partition_function
		integer :: ix, iKcm
		real*8 :: min_energy, deltaE

		min_energy = minval(currcnt%Ex_t)
		deltaE = (-1.d0) * log(1.d-3) * kb * max_temperature
        
		partition_function = 0.d0

		! partition function due to target exciton type
		do ix = 1,currcnt%nX_t
			do iKcm = currcnt%iKcm_min_fine,currcnt%iKcm_max_fine
				if((currcnt%Ex_A1(ix,iKcm)-min_energy) .le. deltaE) then
					partition_function = partition_function + exp(-currcnt%Ex_t(ix,iKcm)/kb/temperature)
				endif
			end do
		end do
        
		return        
	end subroutine calculate_partition_function

				
end module partition_function_mod