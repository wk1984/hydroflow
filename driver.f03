program driver
use ctrl
use ncutils
use files
use hydro


implicit none

type(controller) :: c
character(len=255) :: dir_file, out_file, grid_file, test_mode
character(len=255), allocatable :: in_files(:)
character(len=10) :: n_files
integer :: i, nx, ny, nxny, nt, nf, ncid, varid, test=0
type(ncvar) :: vars(7)
type(ncdim), pointer :: dims(:)
type(ncdim) :: dim(3)
logical :: grid


call get_environment_variable("n_files", n_files)
call get_environment_variable("directions_file", dir_file)
call get_environment_variable("out_file", out_file)
call get_environment_variable("grid_file", grid_file)
call get_environment_variable("test_mode", test_mode)

read(n_files,*) nf
allocate( in_files(nf) )
read (*,*) in_files
inquire(file=trim(grid_file), exist=grid)
if (test_mode/="") read(test_mode,*) test


call vars(1)%init("SHSN3", start=(/1,1,2,1/))
call vars(2)%init("RU", start=(/1,1,2,1/))
call vars(3)%init("discharge_fast", xtype=NF90_FLOAT, dimids=(/1,2,3/))
call vars(4)%init("discharge_slow", xtype=NF90_FLOAT, dimids=(/1,2,3/))
call vars(5)%init("timecoefs_fast", xtype=NF90_FLOAT, dimids=(/1,2,3/))
call vars(6)%init("timecoefs_slow", xtype=NF90_FLOAT, dimids=(/1,2,3/))
call vars(7)%init("runoff", xtype=NF90_FLOAT, dimids=(/1,2,3/))

if (grid) then
	print *,"Starting with existing grid file '" // trim(grid_file) // "'."
	call read_grid(c, grid_file)
	call nc_read(in_files(1), vars(1:2), dims)
else
	call setup_grid(c, in_files(1), dir_file, grid_file, vars(1:2), dims)
endif

nx = dims(1)%len
ny = dims(2)%len
nt = dims(4)%len
nxny = nx*ny

c%swe_depth(1:nx,1:ny,1:nt) => vars(1)%r
c%runoff   (1:nx,1:ny,1:nt) => vars(2)%r


call c%init(nxny)
c%vel_stats = (/10000000.,0.,0.,10000000.,0.,0./)


! Initialize the discharge to zero for first year, then run first year again with 
! end-of-year discharge as initial condition.
c%dis_init_slow = 0.0
c%dis_init_fast = 0.0

print *,"Initialize first year with zero discharge initial conditions."
if (test==0) call hydroflow(c, nx, ny, nt, nxny)

dim(1:2) = dims(1:2)
call dim(3)%init("TIME", len=NF90_UNLIMITED, xtype=dims(4)%xtype)

if (out_file=="") then
	call define_nc(file_name(in_files(1)), vars(3:7), dim, ncid, varid)
else
	call define_nc(out_file, vars(3:7), dim, ncid, varid)
end if

call check( nf90_put_var(ncid, 3, dims(4)%d) )

if (test==0) print *,"Repeat first year with end-of-year discharge as initial condition."
call hydroflow(c, nx, ny, nt, nxny, ncid, varid)


do i=2,nf	
	if (out_file=="") then
		call check( nf90_close(ncid) )
	else
		c%ts = c%ts + nt
	endif
	nullify(c%swe_depth)
	nullify(c%runoff)
	call nc_read(in_files(i), vars(1:2), dims)
	nt = dims(4)%len
	c%swe_depth(1:nx,1:ny,1:nt) => vars(1)%r
	c%runoff   (1:nx,1:ny,1:nt) => vars(2)%r
	
	if (out_file=="") call define_nc(file_name(in_files(i)), vars(3:4), dim, ncid, varid)
	call check( nf90_put_var(ncid, 3, dims(4)%d, start=(/c%ts/)) )
	
	print *,"Working on year ", i
	call hydroflow(c, nx, ny, nt, nxny, ncid, varid)
enddo

if (test.gt.0) then
	c%ts = c%ts + nt
	c%swe_depth(1:nx,1:ny,1) = c%swe_depth(1:nx,1:ny,nt)
	c%runoff   (1:nx,1:ny,1) = 0.
	do i=1,test
		call check( nf90_put_var(ncid, 3, (/dims(4)%d(nt)+i*(dims(4)%d(nt)-dims(4)%d(nt-1))/), start=(/c%ts/)) )
		call hydroflow(c, nx, ny, 1, nxny, ncid, varid)
		c%ts = c%ts + 1
	enddo
endif

call check( nf90_close(ncid) )

do i=1,size(dims)
	call dims(i)%dealloc()
enddo
do i=1,size(vars)
	call vars(i)%dealloc()
enddo
call c%dealloc()

end program driver