module ctrl
implicit none

type controller
	integer, dimension(:,:), pointer :: index_eqn, nwatershed, nsoln_order
	real*4, dimension(:), pointer :: dis_init_slow, dis_init_fast
	real*4, dimension(:,:), pointer :: z, land, ice, dir, slope_corr
	real*4, dimension(:,:,:), pointer :: tcoef_slow, tcoef_fast, swe_depth, runoff, discharge_slow, discharge_fast
	real*4 :: dx, dy, dist, a
	integer, pointer :: nwatershed_count
	double precision, pointer :: vel_stats(:)
	integer :: ts = 1
	real*4 :: ice_msk = 0.5 ! ice is cells > 50%
	real*4 :: land_srf = 4.
	
! 	Constants used in the residence time calculations.  Begin by
!   defining the observed flow velocities (m/s).  We will assume
!   these are for a slope of 15% and let the 'slope_corr' factor
!   calculated below adjust these observed values up or down for
!   steeper and lower slopes, respectively.  The numbers provided
!   here are based on Mittivakkat Glacier field observations:
!   see Liston and Mernild 2012 (the original HydroFlow paper).
	real*4 :: v_slow_sc_ice = 0.12
	real*4 :: v_slow_sf_ice = 0.20
	real*4 :: v_slow_sc_land = 0.10
	real*4 :: v_slow_sf_land = 0.08
	
! 	'alfa' is multiplied by the initial slow and fast velocities to
!   increase or decrease the fluid residence time in each grid cell.
!   You can either just leave alfa_s = 1.0 and adjust the fast alfa
!   (alfa_f) to fit an observed hydrograph (because the slow-
!   component velocity has been defined using observations of some
!   kind and the fast-component velocity is largely unknown), or
!   you can adjust both alfa_s and alfa_f; I don't have enough
!   experience or data to suggest one approach is better than the
!   other.  Also see the notes below about how the 'alfa'
!   adjustments are used in the velocity and residence time
!   calculations.
!     alfa_s = 1.0
!     alfa_f = 1.0
! The alfa's below were found to provide a best-fit to the ERA-I
!   Greenland simulations when compared to the Kangerlussuaq
!   discharge observations.
	real*4 :: alfa_s = 2.0
	real*4 :: alfa_f = 2.0
	
contains
	procedure :: init => ctrl_init
	procedure :: dealloc => ctrl_dealloc
end type controller

contains

subroutine ctrl_init(this,nxny)
	class(controller) :: this
	integer, intent(in) :: nxny
	integer :: err
	allocate( &
		this%vel_stats(6), &
		this%dis_init_fast(nxny), &
		this%dis_init_slow(nxny) &
	)
end subroutine ctrl_init

subroutine ctrl_dealloc(this)
	class(controller) :: this
	integer :: err
	deallocate( &
		this%slope_corr,  &
		this%nwatershed,  &
		this%nsoln_order, &
		this%index_eqn, &
		this%nwatershed_count, &
		this%vel_stats, &
		this%dis_init_fast, &
		this%dis_init_slow &
	)
end subroutine ctrl_dealloc


end module ctrl