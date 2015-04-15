module arbitraryAngleForster_module
	implicit none
	private
	
	public  :: calculateArbitraryForsterRate
	
contains		
	!**************************************************************************************************************************
	! calculate scattering rate from cnt1 to cnt2 when they are perpendicular
	!**************************************************************************************************************************
	
	subroutine calculateArbitraryForsterRate(cnt1, cnt2, totalTransitionRate12, totalTransitionRate21, c2cDistance, theta)
		use physicalConstants, only : kb, pi, hb
		use comparams, only : Temperature, ppLen
		use prepareForster_module
		use cnt_class, only: cnt

		type(cnt), intent(in) :: cnt1, cnt2
		real*8, intent(in) :: theta, c2cDistance
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
				
					
			if ((iKcm1 .ne. 0) .and. (iKcm2 .ne. 0)) then
				call calculateMatrixElement(cnt1,cnt2, iKcm1, iKcm2, iC, matrixElement, c2cDistance, theta)
				call calculateDOS(cnt1,iKcm1,ix1,dos1)
				call calculateDOS(cnt2,iKcm2,ix2,dos2)

				totalTransitionRate12 = totalTransitionRate12 + exp(-(cnt1%Ex_t(ix1,iKcm1))/kb/Temperature) * dble(conjg(matrixElement) * matrixElement) * dos2 * (sin(theta))**2 / hb / ppLen/ (partitionFunction1 / cnt1%dk) !the multiplication of cnt.dk is because the way partitionFunction is calculated it has units of 1/L while it should be unitless. * (1.d0-exp(-15.d0*theta))
				totalTransitionRate21 = totalTransitionRate21 + exp(-(cnt2%Ex_t(ix2,iKcm2))/kb/Temperature) * dble(conjg(matrixElement) * matrixElement) * dos1 * (sin(theta))**2 / hb / ppLen/ (partitionFunction2 / cnt2%dk) !the multiplication of cnt.dk is because the way partitionFunction is calculated it has units of 1/L while it should be unitless. * (1.d0-exp(-15.d0*theta))
			end if

		end do
		return	
	end subroutine calculateArbitraryForsterRate

	!**************************************************************************************************************************
	! calculate the matrix element for the crossing point number iC in two unparallel tube
	!**************************************************************************************************************************
	
	subroutine calculateMatrixElement(cnt1,cnt2, iKcm1, iKcm2, iC ,matrixElementFinal, c2cDistance, theta)
		use physicalConstants, only : i1, pi, eps0, q0
		use prepareForster_module, only : kSpaceMatrixElement
		use cnt_class, only: cnt

		type(cnt), intent(in) :: cnt1,cnt2
		real*8, intent(in) :: theta, c2cDistance
		complex*16, intent(out) :: matrixElementFinal
		integer, intent(in) :: iKcm1, iKcm2
		integer, intent(in) :: iC
	
		real*8 :: dk, K1, K2
		complex*16 :: Jk
		real*8 :: phi1, phi2, dPhi
		integer :: iPhi1, iPhi2, nPhi
		integer :: nx, ny
		real*8 :: dx, dy, x_max, y_max
		real*8 :: radius1, radius2
		real*8 :: arg1, arg2, arg3
				
		radius1 = cnt1%radius
		radius2	= cnt2%radius
		dk = cnt1%dk
		K1 = dble(iKcm1)*dk
		K2 = dble(iKcm2)*dk
        
		nPhi = 20
		dPhi = 2.d0*pi/dble(nPhi)
				
		x_max = 20.d-9
		y_max = x_max
			
		dx = 1.0d-10
		dy = dx
				
		nx = int(x_max/dx)
		ny = nx
				
		Jk = (0.d0,0.d0)
			
		arg1 = sqrt(K1**2+K2**2-2.d0*K1*K2*cos(theta))
		
		do iPhi1 = 1,nPhi
			phi1=dble(iPhi1)*dPhi
			do iPhi2 = 1,nPhi
				phi2=dble(iPhi2)*dPhi
				arg2 = 2.d0 * (K1*(radius2*cos(phi2)-radius1*cos(phi1)*cos(theta))+K2*(radius1*cos(phi1)-radius2*cos(phi2)*cos(theta))) / (sin(theta))
				arg3 = 2.d0 * abs((c2cDistance+radius2*sin(phi2)-radius1*sin(phi1))/(sin(theta)))
				Jk = Jk + exp(i1*dcmplx(arg2)) * dcmplx(exp(-arg3 * arg1))
			end do
		end do
				
		Jk = Jk * dPhi * dPhi * pi / arg1
				
		matrixElementFinal = kSpaceMatrixElement(iC) * Jk * dcmplx(q0*q0 / (4.d0*pi*eps0*4.d0*(pi*pi)*2.d0*pi/dk))
        return
	end subroutine calculateMatrixElement

end module arbitraryAngleForster_module