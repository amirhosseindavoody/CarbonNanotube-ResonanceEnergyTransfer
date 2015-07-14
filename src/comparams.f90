module comparams
	use cnt_class, only : cnt
	implicit none


	real :: starttime,endtime !duration of simulation

	real*8 :: min_temperature, max_temperature
	integer :: temperature_steps

	real*8 :: ppLen = 10.d-9 !length per perpendicular tubes

	type (cnt) :: cnt1, cnt2

end module comparams