module hydro
use ctrl
use ncutils

implicit none

contains 


subroutine timecoef(c, snow, vel_slow, vel_fast)
	type(controller), intent(in) :: c
	real*4, dimension(:,:,:), intent(in) :: snow
	real*4, dimension(:,:,:), intent(out) :: vel_slow, vel_fast
	integer k
	
	do k=1,size(snow,3)
		vel_slow(:,:,k) = &
			c%v_slow_sc_ice * snow(:,:,k) * c%ice + &
			c%v_slow_sf_ice * (1.-snow(:,:,k)) * c%ice + &
			c%v_slow_sc_land * snow(:,:,k) * (1.-c%ice) + &
			c%v_slow_sf_land * (1.-snow(:,:,k)) * (1.-c%ice)
		
		vel_fast(:,:,k) = c%alfa_f * c%slope_corr * vel_slow(:,:,k)
		vel_slow(:,:,k) = c%alfa_s * c%slope_corr * vel_slow(:,:,k)
	enddo
	
	c%vel_stats(1) = min( minval(vel_slow), c%vel_stats(1))
	c%vel_stats(2) = max( maxval(vel_slow), c%vel_stats(2))
	c%vel_stats(3) = c%vel_stats(3) + sum(vel_slow)
	c%vel_stats(4) = min( minval(vel_fast), c%vel_stats(4))
	c%vel_stats(5) = max( maxval(vel_fast), c%vel_stats(5))
	c%vel_stats(6) = c%vel_stats(6) + sum(vel_fast)
end subroutine timecoef



subroutine hydroflow(c, nx, ny, nt, nxny, ncid, varid)
	type(controller), intent(inout) :: c
	integer, intent(in) :: nx,ny,nxny,nt
	integer, intent(in), optional :: ncid, varid
	integer, pointer :: nwatersheds, index_eqn(:,:)
	real*4, dimension(:), pointer :: dis_init_slow, dis_init_fast
	real*4, dimension(nxny,nt) :: roff_input, time_slow, time_fast, r_slow, r_fast
	real*4, dimension(nx,ny,nt) :: runoff, snow, tcoef_slow, tcoef_fast, discharge_slow, discharge_fast
	integer nwshed,inflow_id,i,j,k,n,a,b

	real*4 dt,roff_input_tmp
	
	real*4, parameter :: sec_in_day = 86400.0
	
	! Runoff output time step. This has to be in the same units as
	!   the time coefficent data (like 1 day).
	dt = 1.0

	nwatersheds      => c%nwatershed_count
	index_eqn        => c%index_eqn
	dis_init_slow    => c%dis_init_slow
	dis_init_fast    => c%dis_init_fast
	
	snow = 0.
	do k=1,nt
		do j=1,ny
			do i=1,nx
				if (c%swe_depth(i,j,k).gt.0.) snow(i,j,k) = 1.
			enddo
		enddo
! Convert the runoff values from m/day to m^3/day + plus mask out ocean
		runoff(:,:,k) = c%runoff(:,:,k) * c%dx * c%dy * c%land
	enddo
	
	call timecoef(c, snow, tcoef_slow, tcoef_fast)
	tcoef_slow = c%dist / (tcoef_slow * sec_in_day)
	tcoef_fast = c%dist / (tcoef_fast * sec_in_day)
	

! Solve for the daily discharge one year at a time.

! Since initial conditions are required to calculate the discharge,
!   and these do not exist, the model first assumes the initial
!   discharge conditions in each grid cell is zero.  It then begins
!   by looping once through all of Year-1 with this zero-discharge
!   initial conditdion.  The model then takes the discharge that
!   was simulated on the last day of Year-1 and uses it for the
!   initial condition for a second simulation of Year-1.  The model
!   continues on to Year-2, and all subsequent years, using the
!   discharge from the last day of the previous year as the initial
!   condition for the current simulation year.

! Extract the runoff and time coefficient information for this
!   watershed.
            
	do n=1,nxny
		i = index_eqn(n,1)
		j = index_eqn(n,2)
		roff_input(n,:) = runoff(i,j,:)
		time_slow(n,:) = tcoef_slow(i,j,:)
		time_fast(n,:) = tcoef_fast(i,j,:)
	enddo

! Initialize the discharge arrays to put something in the grid
!   cells where nothing is ever calculated (like ocean grid cells).
	discharge_fast = 0.0
	discharge_slow = 0.0

! SOLVE THE SLOW-FLOW EQUATIONS.
! Note that dt and the time coefficients need to have the same units.
!   The slow equations are not coupled to any other grid cells, so
!   the solution order does not matter.

! Note that r_slow has dimensions of nt+1, and roff_input and
!   time_slow have dimensions of nt, so k-1 is present time for
!   roff_input and time_slow, while k-1 is the previous time step for
!   r_slow.
	
	do k=1,nt
! 		r_slow(:,k+1) = r_slow(:,k) * exp(-dt/time_slow(:,k)) + roff_input(:,k) * (1.-exp(-dt/time_slow(:,k)))
! 		r_slow(:,k+1) = (dis_init_slow-roff_input(:,k)) * (1.-exp(-dt/time_slow(:,k))) * time_slow(:,k)/dt + roff_input(:,k)
		dis_init_slow = dis_init_slow * exp(-dt/time_slow(:,k)) + roff_input(:,k) * (1.-exp(-dt/time_slow(:,k)))
		r_slow(:,k) = dis_init_slow
	enddo
		
! Extract the watershed of interest, solve for that watershed,
!   save the outputs, and move on to the next watershed.

! The solution order for each watershed is given in column 5 of
!   index_eqn, and it is the same as the read-in order for each
!   watershed.
	
	a = 1
	b = 0
	do nwshed=1,nwatersheds

! Extract the indexing data that correspond to the watershed
!  being processed.
		do while (b.lt.nxny .and. index_eqn(b+1,3).eq.nwshed)
			b = b + 1
		enddo


! SOLVE THE FAST-FLOW EQUATIONS.

! The runoff input is now the combination of the slow runoff from the
!   grid cell in question, and the inputs from any upstream grid
!   cells.  Because of this, the equations must be solved in the
!   solution order given by column 5 of the index array.  The indexing
!   array has already been sorted for this order, so it is okay to solve
!   the equations in the order extracted from the index array, 1 --> neq.
!   Note that k is the current time step for r_slow and r_fast.
		do k=1,nt
			do j=a,b
				roff_input_tmp = r_slow(j,k)
				do i=1,index_eqn(j,6)
					! Find the equation solution sequence line that corresponds to this eqn.
					do n=a,b
						if (index_eqn(n,4).eq.index_eqn(j,i+6)) inflow_id = n
					enddo
					roff_input_tmp = roff_input_tmp + r_fast(inflow_id,k)
				enddo
! 				r_fast(j,k) = (dis_init_fast(j)-roff_input_tmp) * (1. - exp(-dt/time_fast(j,k))) * time_fast(j,k)/dt + roff_input_tmp
				dis_init_fast(j) = dis_init_fast(j) * exp(-dt/time_fast(j,k)) + roff_input_tmp * (1. - exp(-dt/time_fast(j,k)))
				r_fast(j,k) = dis_init_fast(j)
			enddo
		enddo

! Build the equation-solution matrix, discharge(i,j,k).
		do n=a,b
			i = index_eqn(n,1)
			j = index_eqn(n,2)
			discharge_slow(i,j,:) = r_slow(n,:)
			discharge_fast(i,j,:) = r_fast(n,:)
		enddo
	
	a = b + 1
	enddo ! nwshed


	if (present(ncid)) then
		call check( nf90_put_var(ncid, varid,   discharge_fast, start=(/1,1,c%ts/)) )
		call check( nf90_put_var(ncid, varid+1, discharge_slow, start=(/1,1,c%ts/)) )
		call check( nf90_put_var(ncid, varid+2, tcoef_fast, start=(/1,1,c%ts/)) )
		call check( nf90_put_var(ncid, varid+3, tcoef_slow, start=(/1,1,c%ts/)) )
		call check( nf90_put_var(ncid, varid+4, runoff, start=(/1,1,c%ts/)) )
	endif

end subroutine hydroflow

end module hydro