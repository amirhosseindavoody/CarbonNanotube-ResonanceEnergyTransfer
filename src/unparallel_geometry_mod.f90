module unparallel_geometry_mod
	implicit none
	private
	
	public  :: calculateUnparallelGeometryRate
	
contains		
	!**************************************************************************************************************************
	! calculate scattering rate from cnt1 to cnt2 when they are perpendicular
	!**************************************************************************************************************************
	
	subroutine calculateUnparallelGeometryRate(cnt1, cnt2, totalTransitionRate12, totalTransitionRate21, theta)
		use cnt_class, only: cnt
		use comparams, only : Temperature, ppLen
		use matrix_element_mod, only: geometricMatrixElement, kSpaceMatrixElement
		use physicalConstants, only : kb, pi, hb
		use prepareForster_module, only: calculatePartitionFunction, calculateDOS
		use transition_points_mod, only: sameEnergy

		type(cnt), intent(in) :: cnt1, cnt2
		real*8, intent(in) :: theta
		real*8, intent(out) :: totalTransitionRate12, totalTransitionRate21
		integer :: iC
		integer :: ix1, ix2, iKcm1, iKcm2
		integer :: nSameEnergy
		real*8 :: partitionFunction1, partitionFunction2
		real*8 :: dos1, dos2
		complex*16 :: matrixElement
				
		call calculatePartitionFunction(cnt1, partitionFunction1)
		call calculatePartitionFunction(cnt2, partitionFunction2)

		totalTransitionRate12 = 0.d0
		totalTransitionRate21 = 0.d0
		nSameEnergy = size(sameEnergy,1)
							
		do iC = 1,nSameEnergy
				
			ix1 = sameEnergy(iC,1)
			ix2 = sameEnergy(iC,2)
			iKcm1 = sameEnergy(iC,3)
			iKcm2 = sameEnergy(iC,4)
			
			matrixElement = geometricMatrixElement(iKcm1, iKcm2) * kSpaceMatrixElement(iC)
			call calculateDOS(cnt1,iKcm1,ix1,dos1)
			call calculateDOS(cnt2,iKcm2,ix2,dos2)

			totalTransitionRate12 = totalTransitionRate12 + exp(-(cnt1%Ex_t(ix1,iKcm1))/kb/Temperature) * abs(matrixElement) * dos2 * (sin(theta))**2 / hb / ppLen/ partitionFunction1 !the multiplication of cnt.dk is because the way partitionFunction is calculated it has units of 1/L while it should be unitless. * (1.d0-exp(-15.d0*theta))
			totalTransitionRate21 = totalTransitionRate21 + exp(-(cnt2%Ex_t(ix2,iKcm2))/kb/Temperature) * abs(matrixElement) * dos1 * (sin(theta))**2 / hb / ppLen/ partitionFunction2 !the multiplication of cnt.dk is because the way partitionFunction is calculated it has units of 1/L while it should be unitless. * (1.d0-exp(-15.d0*theta))

		end do
		return	
	end subroutine calculateUnparallelGeometryRate

end module unparallel_geometry_mod