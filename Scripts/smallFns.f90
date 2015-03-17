module smallFunctions
	implicit none
	character(len=200) :: logInput
	contains
		!**********************************************************************************************************************
		! This subroutine calculates the greatest common divisor of the arguments na and nb
		!**********************************************************************************************************************
		subroutine gcd(ngcd,na,nb)
		  	integer, intent(in) :: na,nb
		   	integer, intent(out) :: ngcd
		   	integer :: ia,ib,itemp

		   	ia=na
		    ib=nb
		    do while (ib .ne. 0)
		       	itemp=ia
		       	ia=ib
		       	ib=mod(itemp,ib)
		    end do
		    ngcd=ia
		    return
		end subroutine gcd
			
		!**********************************************************************************************************************
		! This function calculates the modified bessel function of the first kind with parameter nu=0: I0(x)
		!**********************************************************************************************************************
		real*8 function bessi0(x)
			real*8 :: x
			real*8 :: ax
			real*8 :: y
			real*8, save :: p1, p2, p3, p4, p5, p6, p7, q1, q2, q3, q4, q5, q6, q7, q8, q9
				
			data p1, p2, p3, p4, p5, p6, p7 /1.d0, 3.5156229d0, 3.0899424d0, 1.2067492d0, 0.2659732d0, 0.360768d-1, 0.45813d-2/
			data q1, q2, q3, q4, q5, q6, q7, q8, q9 /0.39894228d0, 0.1328592d-1, 0.225319d-2, -0.157565d-2, 0.916281d-2, -0.2057706d-1, 0.2635537d-1, -0.1647633d-1, 0.392377d-2/
			
			if (abs(x) .lt. 3.75d0) then
				y = (x/3.75d0)**2
				bessi0 = p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))
			else
				ax = abs(x)
				y = 3.75d0/ax
				bessi0 = (exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
			end if
			return
		end function bessi0

		!**********************************************************************************************************************
		! This function calculates the modified bessel function of the second kind with parameter nu=0: K0(x)
		!**********************************************************************************************************************
		real*8 function bessk0(x)
			real*8 :: x
			real*8 :: y
			real*8, save :: p1, p2, p3, p4, p5, p6, p7, q1, q2, q3, q4, q5, q6, q7
			
			data p1, p2, p3, p4, p5, p6, p7 /-0.57721566d0, 0.42278420d0, 0.23069756d0, 0.3488590d-1, 0.262698d-2, 0.10750d-3, 0.74d-5/
			data q1, q2, q3, q4, q5, q6, q7 /1.25331414d0, -0.7832358d-1, 0.2189568d-1, -0.1062446d-1, 0.587872d-2, -0.251540d-2, 0.53208d-3/

			if (x .le. 2.d0) then
				y = x*x/4.d0
				bessk0=(-log(x/2.d0)*bessi0(x))+(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
			else
				y=2.d0/x
				bessk0=(exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*q7))))))
			end if
			return 
		end function bessk0

		!*******************************************************************************
		! This subroutines opens the log file add new log and closes the file
		!*******************************************************************************
		subroutine writeLog()
			implicit none
			logical :: exist
			integer :: logFile = 10

			inquire(file="log.dat",exist=exist)
			if (exist) then
				open(logFile, file="log.dat", status="old", position="append", action="write")
			else
				open(logFile, file="log.dat", status="new", action="write")
			end if
			write(logFile, *) logInput
			print *, logInput
			close(logFile)

			return
		end subroutine writeLog
end module smallFunctions