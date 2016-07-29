module files
use ctrl
use ncutils
use dem


implicit none

contains


subroutine read_grid(c, grid_file)
	class(controller), intent(inout) :: c
	character(len=*), intent(in) :: grid_file
	type(ncvar) :: vars(7)
	type(ncdim), pointer :: dims(:)
	integer :: nx, ny, nxny
	
	call vars(1)%init("land")
	call vars(2)%init("ice")
	call vars(3)%init("nwatershed")
	call vars(4)%init("nsoln_order")
	call vars(5)%init("slope_corr")
	call vars(6)%init("index_eqn")
	call vars(7)%init("watershed_count")
	
	call nc_read(grid_file, vars, dims)
	
	nx = dims(1)%len
	ny = dims(2)%len
		
	c%land       (1:nx,1:ny)     => vars(1)%r
	c%ice        (1:nx,1:ny)     => vars(2)%r
	c%nwatershed (1:nx,1:ny)     => vars(3)%i
	c%nsoln_order(1:nx,1:ny)     => vars(4)%i
	c%slope_corr (1:nx,1:ny)     => vars(5)%r
	c%index_eqn  (1:nx*ny, 1:14) => vars(6)%i
	c%nwatershed_count       => vars(7)%i(1)
	
	c%dx = abs( dims(1)%d(2)-dims(1)%d(1) ) * 1000.
	c%dy = abs( dims(2)%d(2)-dims(2)%d(1) ) * 1000.
	c%dist  = (c%dx + c%dy) / 2.0
end subroutine read_grid


subroutine setup_grid(c, file, dir_file, grid_file, vars, dims)
	class(controller), intent(inout) :: c
	character(len=*), intent(in) :: file, dir_file, grid_file
	type(ncvar), target, intent(inout) :: vars(:)
	type(ncdim), pointer, intent(out) :: dims(:)
	type(ncdim) :: d(4)
	type(ncvar) :: out(7), v(6)
	integer :: i, j, ncid, nx, ny, nxny
	
	v(1:2) = vars(1:2)
	call v(3)%init("SH")
	call v(4)%init("SRF")
	call v(5)%init("MSK")
	call v(6)%init("Band1")
	
	call nc_read(file, v(1:5), dims)
	call nc_read(dir_file, v(6:6))
	vars(1:2) = v(1:2)
	
	nx = dims(1)%len
	ny = dims(2)%len
	nxny = nx*ny

	c%z   (1:nx,1:ny) => v(3)%r	
	c%land(1:nx,1:ny) => v(4)%r
	c%ice (1:nx,1:ny) => v(5)%r
	c%dir (1:nx,1:ny) => v(6)%r
	
	! vegetation and surface format translation
	do j=1,ny
		do i=1,nx
			if (c%land(i,j).eq.c%land_srf) then
				c%land(i,j) = 1.
			else
				c%land(i,j) = 0.
				c%z(i,j) = 0.
			endif
			if (c%ice(i,j).gt.50.) then
				c%ice(i,j) = 1.
			else
				c%ice(i,j) = 0.
			endif
		enddo
	enddo
	
	! 'ANSWERS' direction format translation
	c%dir = 2**modulo(9 - c%dir/45.,8.)
	
	allocate( &
		c%slope_corr(nx,ny),  &
		c%nwatershed(nx,ny),  &
		c%nsoln_order(nx,ny), &
		c%index_eqn(nxny,14), &
		c%nwatershed_count &
	)
	
	c%dx = abs( dims(1)%d(2)-dims(1)%d(1) ) * 1000.
	c%dy = abs( dims(2)%d(2)-dims(2)%d(1) ) * 1000.
	c%dist  = (c%dx + c%dy) / 2.0
	
	print *,"Computing watersheds and slope corrections."
	
	! NOTE: modifies directions (dir)
	call watersheds(c, nx, ny, nxny)
	call slope_calc(c, nx, ny)
	
	if (grid_file/="") then
		print *,"Saving watersheds and slope corrections in file '" // trim(grid_file) // "'."
		
		d(1:2) = dims(1:2)
		call d(3)%init("s_ord", len=nxny)
		call d(4)%init("s_col", len=14)
	
		call out(1)%init("land", xtype=NF90_FLOAT, dimids=(/1,2/))
		call out(2)%init("ice", xtype=NF90_FLOAT, dimids=(/1,2/))
		call out(3)%init("nwatershed", xtype=NF90_INT, dimids=(/1,2/))
		call out(4)%init("nsoln_order", xtype=NF90_INT, dimids=(/1,2/))
		call out(5)%init("slope_corr", xtype=NF90_FLOAT, dimids=(/1,2/))
		call out(6)%init("index_eqn", xtype=NF90_INT, dimids=(/3,4/))
		call out(7)%init("watershed_count", xtype=NF90_INT)
	
		call define_nc(grid_file, out, d, ncid, i)
	
		call check( nf90_put_var(ncid, i,   c%land) )
		call check( nf90_put_var(ncid, i+1, c%ice) )
		call check( nf90_put_var(ncid, i+2, c%nwatershed) )
		call check( nf90_put_var(ncid, i+3, c%nsoln_order) )
		call check( nf90_put_var(ncid, i+4, c%slope_corr) )
		call check( nf90_put_var(ncid, i+5, c%index_eqn) )
		call check( nf90_put_var(ncid, i+6, c%nwatershed_count) )
	
		call check( nf90_close(ncid) )
	endif
end subroutine setup_grid


function file_name(in) result(out)
	implicit none
	character(len=*), intent(in) :: in
	character(len=:), allocatable :: out
	character(len=255) :: str
	str = in(index(in, '/', .true.)+1:)
	out = 'DIS' // trim(str(index(str, '.'):))
end function file_name

end module files