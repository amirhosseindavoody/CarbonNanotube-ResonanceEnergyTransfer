!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate transition table
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
module transition_table_mod
	implicit none
	private
	public :: calculate_transition_table, save_transition_rates

	real*8, dimension(:,:,:), allocatable :: transitionRate	

	real*8, public :: c2cDistance
	real*8 :: theta

	complex*16, dimension(:), allocatable :: kSpaceMatrixElement_sameEnergy, kSpaceMatrixElement_crossoingPoints

	abstract interface 
		subroutine calculate_kSpaceMatrixElement(nTransitionPoints, transitionPoints, kSpaceMatrixElement)
			integer, intent(in) :: nTransitionPoints
			integer, dimension(nTransitionPoints,4), intent(in) :: transitionPoints
			complex*16, dimension(:), allocatable, intent(inout) :: kSpaceMatrixElement
		end subroutine calculate_kSpaceMatrixElement
	end interface
	
contains
	!**************************************************************************************************************************
	! calculate transition table
	!**************************************************************************************************************************
	
	subroutine calculate_transition_table (cnt1, cnt2)
		use a2a_kspace_matrix_element_mod, only: calculate_a2a_kSpaceMatrixElement
		use a2ep_kspace_matrix_element_mod, only: calculate_a2ep_kSpaceMatrixElement
		use a2em_kspace_matrix_element_mod, only: calculate_a2em_kSpaceMatrixElement
		use cnt_class, only: cnt
		use comparams, only: ppLen, min_temperature, max_temperature, temperature_steps
		use em2a_kspace_matrix_element_mod, only: calculate_em2a_kSpaceMatrixElement
		use em2ep_kspace_matrix_element_mod, only: calculate_em2ep_kSpaceMatrixElement
		use em2em_kspace_matrix_element_mod, only: calculate_em2em_kSpaceMatrixElement
		use ep2a_kspace_matrix_element_mod, only: calculate_ep2a_kSpaceMatrixElement
		use ep2ep_kspace_matrix_element_mod, only: calculate_ep2ep_kSpaceMatrixElement
		use ep2em_kspace_matrix_element_mod, only: calculate_ep2em_kSpaceMatrixElement
		use geometric_matrix_element_mod, only: calculate_finite_geometric_matrix_element, calculate_infinite_geometric_matrix_element, calculate_infinite_parallel_geometric_matrix_element
		use physicalConstants, only: pi, kb, hb, A_u
		use prepareForster_module, only: calculatePartitionFunction, calculateDOS
		use rotate_shift_mod, only: rotate_shift_cnt
		use transition_points_mod, only: findSameEnergy, sameEnergy, findCrossings, crossingPoints
		use write_log_mod, only: writeLog

		type(cnt), intent(inout) :: cnt1,cnt2
		character(len=200) :: logInput

		integer :: nSameEnergy, iS
		integer :: nCrossing, iC
		integer :: ix1, ix2, iKcm1, iKcm2
		integer :: iTemperature
		real*8 :: temperature
		real*8 :: partitionFunction1, partitionFunction2
		real*8 :: dos1, dos2
		complex*16 :: matrixElement, geometricMatrixElement

		procedure(calculate_kSpaceMatrixElement), pointer :: k_space_melement_ptr => null()

		
		select case (trim(cnt1%targetExcitonType))
		case('Ex_A1', 'Ex0_A2', 'Ex1_A2')
			select case (trim(cnt2%targetExcitonType))
			case('Ex_A1', 'Ex0_A2', 'Ex1_A2')
				k_space_melement_ptr => calculate_a2a_kSpaceMatrixElement
			case('Ex0_Ep', 'Ex1_Ep')
				k_space_melement_ptr => calculate_a2ep_kSpaceMatrixElement
			case('Ex0_Em', 'Ex1_Em')
				k_space_melement_ptr => calculate_a2em_kSpaceMatrixElement
			case default
				write(*,*) "Could not recognize target exciton type!!!"
				call exit()
			end select
		case('Ex0_Ep', 'Ex1_Ep')
			select case (trim(cnt2%targetExcitonType))
			case('Ex_A1', 'Ex0_A2', 'Ex1_A2')
				k_space_melement_ptr => calculate_ep2a_kSpaceMatrixElement
			case('Ex0_Ep', 'Ex1_Ep')
				k_space_melement_ptr => calculate_ep2ep_kSpaceMatrixElement
			case('Ex0_Em', 'Ex1_Em')
				k_space_melement_ptr => calculate_ep2em_kSpaceMatrixElement
			case default
				write(*,*) "Could not recognize target exciton type!!!"
				call exit()
			end select
		case('Ex0_Em', 'Ex1_Em')
			select case (trim(cnt2%targetExcitonType))
			case('Ex_A1', 'Ex0_A2', 'Ex1_A2')
				k_space_melement_ptr => calculate_em2a_kSpaceMatrixElement
			case('Ex0_Ep', 'Ex1_Ep')
				k_space_melement_ptr => calculate_em2ep_kSpaceMatrixElement
			case('Ex0_Em', 'Ex1_Em')
				k_space_melement_ptr => calculate_em2em_kSpaceMatrixElement
			case default
				write(*,*) "Could not recognize target exciton type!!!"
				call exit()
			end select
		case default
			write(*,*) "Could not recognize target exciton type!!!"
			call exit()
		end select

		!calculate the crossing points and points with the same energy between cnt1 and cnt2
		call findCrossings(cnt1,cnt2)
		call findSameEnergy(cnt1,cnt2)

		call k_space_melement_ptr(size(sameEnergy,1), sameEnergy, kSpaceMatrixElement_sameEnergy)
		call k_space_melement_ptr(size(crossingPoints,1), crossingPoints, kSpaceMatrixElement_crossoingPoints)

		!allocate the transition rate table for the first time
		if (.not. allocated(transitionRate)) then
			allocate(transitionRate(2,2,temperature_steps))
			transitionRate = 0.d0
		endif
		
		! calculate transfer rate for parallel orientation
		theta = 0.d0
				
		! calculate exciton transfer rate for finite length CNTs				
		if ((cnt1%length .lt. huge(1.d0)) .and. (cnt2%length .lt. huge(1.d0))) then
			write(logInput,*) 'Calculating finite transition rate: theta = 0 , c2cDistance = ', c2cDistance
			call writeLog(logInput)

			nSameEnergy = size(sameEnergy,1)

			call rotate_shift_cnt(cnt1, 0.d0, 0.d0)
			call rotate_shift_cnt(cnt2, theta, c2cDistance)
									
			do iS = 1,nSameEnergy
													
				ix1 = sameEnergy(iS,1)
				ix2 = sameEnergy(iS,2)
				iKcm1 = sameEnergy(iS,3)
				iKcm2 = sameEnergy(iS,4)

				call calculate_finite_geometric_matrix_element(iKcm1, iKcm2, geometricMatrixElement)
							
				matrixElement = geometricMatrixElement * kSpaceMatrixElement_sameEnergy(iS)
				call calculateDOS(cnt1,iKcm1,ix1,dos1)
				call calculateDOS(cnt2,iKcm2,ix2,dos2)

				do iTemperature = 1,temperature_steps
					temperature = min_temperature + real(iTemperature-1) * (max_temperature-min_temperature) / real(temperature_steps-1)

					call calculatePartitionFunction(cnt1, temperature, partitionFunction1)
					call calculatePartitionFunction(cnt2, temperature, partitionFunction2)

					transitionRate(1,1,iTemperature) = transitionRate(1,1,iTemperature) + exp(-(cnt1%Ex_t(ix1,iKcm1))/kb/temperature) * (abs(matrixElement)**2) * dos2 * (A_u**2/(4.d0*pi*pi*cnt1%radius*cnt2%radius))**2 / hb / cnt1%length / partitionFunction1 
					transitionRate(2,1,iTemperature) = transitionRate(2,1,iTemperature) + exp(-(cnt2%Ex_t(ix2,iKcm2))/kb/temperature) * (abs(matrixElement)**2) * dos1 * (A_u**2/(4.d0*pi*pi*cnt1%radius*cnt2%radius))**2 / hb / cnt2%length / partitionFunction2
				enddo

			end do

		! calculate exciton transfer rate for infinitely long CNTs
		else
			write(logInput,*) 'Calculating infinite transition rate: theta = 0 , c2cDistance = ', c2cDistance
			call writeLog(logInput)

			nCrossing = size(crossingPoints,1)

			do iC = 1,nCrossing
							
				ix1 = crossingPoints(iC,1)
				ix2 = crossingPoints(iC,2)
				iKcm1 = crossingPoints(iC,3)
				iKcm2 = crossingPoints(iC,4)

				call calculate_infinite_parallel_geometric_matrix_element(iKcm1,iKcm2,c2cDistance,geometricMatrixElement)

				matrixElement = geometricMatrixElement * kSpaceMatrixElement_crossoingPoints(iC)
				call calculateDOS(cnt1,iKcm1,ix1,dos1)
				call calculateDOS(cnt2,iKcm2,ix2,dos2)

				do iTemperature = 1,temperature_steps
					temperature = min_temperature + real(iTemperature-1) * (max_temperature-min_temperature) / real(temperature_steps-1)

					call calculatePartitionFunction(cnt1, temperature, partitionFunction1)
					call calculatePartitionFunction(cnt2, temperature, partitionFunction2)

					transitionRate(1,1,iTemperature) = transitionRate(1,1,iTemperature) + exp(-(cnt1%Ex_t(ix1,iKcm1))/kb/temperature) * (abs(matrixElement)**2) * dos1 * (A_u**2/(4.d0*pi*pi*cnt1%radius*cnt2%radius))**2 * (2*pi/cnt1%dkx) / hb / partitionFunction1
					transitionRate(2,1,iTemperature) = transitionRate(2,1,iTemperature) + exp(-(cnt2%Ex_t(ix2,iKcm2))/kb/temperature) * (abs(matrixElement)**2) * dos2 * (A_u**2/(4.d0*pi*pi*cnt1%radius*cnt2%radius))**2 * (2*pi/cnt2%dkx) / hb / partitionFunction2
				enddo

			end do
		end if
		

		theta = pi/2.d0
				
		! calculate exciton transfer rate for finite length CNTs				
		if ((cnt1%length .lt. huge(1.d0)) .and. (cnt2%length .lt. huge(1.d0))) then
			write(logInput,*) 'Calculating finite transition rate: theta = 90, c2cDistance = ', c2cDistance
			call writeLog(logInput)

			nSameEnergy = size(sameEnergy,1)

			call rotate_shift_cnt(cnt1, 0.d0, 0.d0)
			call rotate_shift_cnt(cnt2, theta, c2cDistance)
									
			do iS = 1,nSameEnergy
													
				ix1 = sameEnergy(iS,1)
				ix2 = sameEnergy(iS,2)
				iKcm1 = sameEnergy(iS,3)
				iKcm2 = sameEnergy(iS,4)

				call calculate_finite_geometric_matrix_element(iKcm1, iKcm2, geometricMatrixElement)
							
				matrixElement = geometricMatrixElement * kSpaceMatrixElement_sameEnergy(iS)
				call calculateDOS(cnt1,iKcm1,ix1,dos1)
				call calculateDOS(cnt2,iKcm2,ix2,dos2)

				do iTemperature = 1,temperature_steps
					temperature = min_temperature + real(iTemperature-1) * (max_temperature-min_temperature) / real(temperature_steps-1)

					call calculatePartitionFunction(cnt1, temperature, partitionFunction1)
					call calculatePartitionFunction(cnt2, temperature, partitionFunction2)

					transitionRate(1,2,iTemperature) = transitionRate(1,2,iTemperature) + exp(-(cnt1%Ex_t(ix1,iKcm1))/kb/temperature) * (abs(matrixElement)**2) * dos2 * (A_u**2/(4.d0*pi*pi*cnt1%radius*cnt2%radius))**2 / hb / cnt1%length / partitionFunction1 
					transitionRate(2,2,iTemperature) = transitionRate(2,2,iTemperature) + exp(-(cnt2%Ex_t(ix2,iKcm2))/kb/temperature) * (abs(matrixElement)**2) * dos1 * (A_u**2/(4.d0*pi*pi*cnt1%radius*cnt2%radius))**2 / hb / cnt2%length / partitionFunction2
				enddo		

			end do
		else
			nSameEnergy = size(sameEnergy,1)
								
			do iS = 1,nSameEnergy
						
				ix1 = sameEnergy(iS,1)
				ix2 = sameEnergy(iS,2)
				iKcm1 = sameEnergy(iS,3)
				iKcm2 = sameEnergy(iS,4)

				call calculate_infinite_geometric_matrix_element(iKcm1, iKcm2, theta, c2cDistance, geometricMatrixElement)
						
				matrixElement = geometricMatrixElement * kSpaceMatrixElement_sameEnergy(iS)
				call calculateDOS(cnt1,iKcm1,ix1,dos1)
				call calculateDOS(cnt2,iKcm2,ix2,dos2)

				do iTemperature = 1,temperature_steps
					temperature = min_temperature + real(iTemperature-1) * (max_temperature-min_temperature) / real(temperature_steps-1)

					call calculatePartitionFunction(cnt1, temperature, partitionFunction1)
					call calculatePartitionFunction(cnt2, temperature, partitionFunction2)

					transitionRate(1,2,iTemperature) = transitionRate(1,2,iTemperature) + exp(-(cnt1%Ex_t(ix1,iKcm1))/kb/temperature) * (abs(matrixElement)**2) * dos2 * (A_u**2/(4.d0*pi*pi*cnt1%radius*cnt2%radius))**2 * (sin(theta)) / hb / ppLen/ partitionFunction1
					transitionRate(2,2,iTemperature) = transitionRate(2,2,iTemperature) + exp(-(cnt2%Ex_t(ix2,iKcm2))/kb/temperature) * (abs(matrixElement)**2) * dos1 * (A_u**2/(4.d0*pi*pi*cnt1%radius*cnt2%radius))**2 * (sin(theta)) / hb / ppLen/ partitionFunction2
				enddo
			end do
		endif
		
		k_space_melement_ptr => null()

		return				
	end subroutine calculate_transition_table
	
	!**************************************************************************************************************************
	! save the calculated transition table
	!**************************************************************************************************************************
	
	subroutine save_transition_rates()

		use comparams, only: min_temperature, max_temperature, temperature_steps

		integer :: iTemperature
		real*8 ::  temperature

		open(unit=100,file='transition_rates.dat',status="unknown")
		do iTemperature = 1,temperature_steps
			temperature = min_temperature + real(iTemperature-1) * (max_temperature-min_temperature) / real(temperature_steps-1)
			write(100,'(E16.8)', advance='no') temperature
			write(100,'(E16.8)', advance='no') transitionRate(1,1,iTemperature)
			write(100,'(E16.8)', advance='no') transitionRate(2,1,iTemperature)
			write(100,'(E16.8)', advance='no') transitionRate(1,2,iTemperature)
			write(100,'(E16.8)', advance='no') transitionRate(2,2,iTemperature)
			write(100,*)
		enddo
		
		return
	end subroutine save_transition_rates

end module transition_table_mod