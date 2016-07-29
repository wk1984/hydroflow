module ncutils
use netcdf

implicit none



type :: var
	character(len=:), allocatable :: name
	integer :: xtype = -1
	integer, pointer :: i(:) => null()
	real*4, pointer :: r(:) => null()
	double precision, pointer :: d(:) => null()
	contains
		procedure :: init => nc_init
		procedure :: dealloc => nc_dealloc
end type var

type, extends(var) :: ncdim
	integer :: dimid, len
	logical :: record_dim = .false.
end type ncdim

type, extends(var) :: ncvar
	integer, allocatable :: dimids(:), shape(:), start(:), count(:)
end type ncvar


contains


subroutine nc_init(this, name, dimids, start, count, len, xtype)
	class(var) :: this
	character(len=*), intent(in) :: name
	integer, intent(in), optional :: dimids(:), start(:), count(:), xtype, len
	this%name = name
	if (present(xtype)) this%xtype = xtype
	select type(this)
	type is (ncvar)
		if (present(count)) allocate( this%count(size(count)), source=count )
		if (present(dimids)) allocate( this%dimids(size(dimids)), source=dimids )
		if (present(start)) allocate( this%start(size(start)), source=start )
	type is (ncdim)
		if (present(len)) this%len = len
	end select
end subroutine nc_init


subroutine nc_dealloc(this)
	class(var) :: this
	deallocate( this%name )
	if (associated(this%i)) deallocate( this%i )
	if (associated(this%r)) deallocate( this%r )
	if (associated(this%d)) deallocate( this%d )
	select type(this)
	type is (ncvar)
		if (allocated(this%dimids)) deallocate( this%dimids )
		if (allocated(this%shape)) deallocate( this%shape )
		if (allocated(this%start)) deallocate( this%start )
		if (allocated(this%count)) deallocate( this%count )
	end select
end subroutine nc_dealloc



! don't forget to deallocate out%data and dims%data
subroutine nc_read(in_file, vars, dims, id)
	character(len=*), intent(in) :: in_file
	type(ncvar), intent(inout) :: vars(:)
	type(ncdim), intent(out), pointer, optional :: dims(:)
	integer, intent(out), optional :: id
	type(ncdim), allocatable :: d(:)
	character(len=NF90_MAX_NAME) :: dname
	integer :: ncid, did, vid, v, i, j, k, ulid, ndims, nAtts, status
	integer :: dimids(NF90_MAX_DIMS), len(NF90_MAX_DIMS), udims(NF90_MAX_DIMS)
	
	! clean these variables because of fortran's automatic save attribute
	k = 0
	len = -1
	udims = -1
	
	call check( nf90_open(trim(in_file), nf90_nowrite, ncid) )
	call check( nf90_Inquire(ncid, nDimensions=ndims, unlimitedDimId=ulid) )
	
	if (present(dims)) allocate( d(ndims) )
	
	do v=1,size(vars)
		call check( nf90_inq_varid(ncid, trim(vars(v)%name), vid) )		
		call check( nf90_Inquire_Variable(ncid, vid, xtype=vars(v)%xtype, ndims=ndims, dimids=dimids, nAtts=nAtts) )
		
		if (allocated(vars(v)%dimids)) deallocate(vars(v)%dimids)
		allocate( vars(v)%dimids(ndims), source=dimids(1:ndims) )
		if (allocated(vars(v)%shape)) deallocate(vars(v)%shape)
		allocate( vars(v)%shape(ndims) )

		do j=1,ndims
			i = dimids(j)
			if (len(i) == -1) then
				if (present(dims)) then
					call check( nf90_Inquire_Dimension(ncid, i, dname, len(i)) )
					d(i)%name = trim(dname)
					status = nf90_inq_varid(ncid, trim(dname), did)
					if (status == nf90_noerr) then 
						call check( nf90_Inquire_Variable(ncid, did, xtype=d(i)%xtype) )
						call get_values(ncid, d(i), did, 1, len(i))
					endif
				else
					call check( nf90_Inquire_Dimension(ncid, i, len=len(i)) )
				endif
				k = k + 1
				udims(k) = i
			endif
			vars(v)%shape(j) = len(i)
		enddo
		
		call get_values(ncid, vars(v), vid, ndims, product(vars(v)%shape))
	enddo
	
	if (present(dims)) then
		if (associated(dims)) deallocate( dims )
		allocate( dims(k), source=d(udims(1:k)) )
		do v=1,size(vars)
			do i=1,size(vars(v)%dimids)
				do j=1,k
					if (vars(v)%dimids(i)==udims(j)) then
						vars(v)%dimids(i) = j
						exit
					endif
				enddo
			enddo
		enddo
		deallocate( d )
	endif
	
	if (present(id)) then 
		id = ncid
	else
		call check( nf90_close(ncid) )
	endif
end subroutine nc_read


subroutine get_values(ncid, v, vid, ndims, n)
	class(var), intent(inout) :: v
	integer, intent(in) :: ncid, vid, n, ndims
	integer, dimension(ndims) :: start, count

	start = 1
	select type (v)
	type is (ncvar)
		if (allocated(v%start)) start = v%start
		if (allocated(v%count)) then
			count = v%count
		else
			count = v%shape - start + 1
		endif
		v%shape = count
	type is (ncdim)
		count = n
		v%len = n
	end select

	if (v%xtype==NF90_INT) then
		if (associated( v%i )) deallocate( v%i )
		allocate( v%i(n) )
		call check( nf90_get_var(ncid, vid, v%i, start=start, count=count ) )		
	elseif (v%xtype==NF90_FLOAT) then
		if (associated( v%r )) deallocate( v%r )
		allocate( v%r(n) )
		call check( nf90_get_var(ncid, vid, v%r, start=start, count=count ) )
	elseif (v%xtype==NF90_DOUBLE) then
		if (associated( v%d )) deallocate( v%d )
		allocate( v%d(n))
		call check( nf90_get_var(ncid, vid, v%d, start=start, count=count ) )
	endif
end subroutine get_values


subroutine define_nc(out_file, vars, dims, ncid, vid)
	character(len=*), intent(in) :: out_file
	type(ncvar), intent(in) :: vars(:)
	type(ncdim), intent(in) :: dims(:)
	integer, intent(out) :: ncid
	integer, intent(out), optional :: vid
	integer :: i, id
	
	call check( nf90_create(trim(out_file), nf90_clobber, ncid) )
	
	do i=1,size(dims)
		call check( nf90_def_dim(ncid, trim(dims(i)%name), dims(i)%len, id) )
		if (dims(i)%xtype.ge.0) &
			call check( nf90_def_var(ncid, trim(dims(i)%name), dims(i)%xtype, id, vid) )
	enddo
	
	do i=1,size(vars)
		if (allocated(vars(i)%dimids)) then
			call check( nf90_def_var(ncid, trim(vars(i)%name), vars(i)%xtype, vars(i)%dimids, id) )
		else
			call check( nf90_def_var(ncid, trim(vars(i)%name), vars(i)%xtype, varid=id) )
		endif
	enddo
	
	call check( nf90_enddef(ncid) )
	
	do i=1,size(dims)
		if (dims(i)%xtype==NF90_INT .and. associated(dims(i)%i)) then
			call check( nf90_put_var(ncid, i, dims(i)%i) )
		elseif (dims(i)%xtype==NF90_FLOAT .and. associated(dims(i)%r)) then
			call check( nf90_put_var(ncid, i, dims(i)%r) )
		elseif (dims(i)%xtype==NF90_DOUBLE .and. associated(dims(i)%d)) then
			call check( nf90_put_var(ncid, i, dims(i)%d) )
		endif
	enddo
	
	vid = vid + 1
	do i=1,size(vars)
		if (vars(i)%xtype==NF90_INT .and. associated(vars(i)%i)) then
			call check( nf90_put_var(ncid, vid, vars(i)%i, count=vars(i)%shape) )
			vid = vid + 1
		elseif (vars(i)%xtype==NF90_FLOAT .and. associated(vars(i)%r)) then
			call check( nf90_put_var(ncid, vid, vars(i)%r, count=vars(i)%shape) )
			vid = vid + 1
		elseif (vars(i)%xtype==NF90_DOUBLE .and. associated(vars(i)%d)) then
			call check( nf90_put_var(ncid, vid, vars(i)%d, count=vars(i)%shape) )
			vid = vid + 1
		endif
	enddo
end subroutine define_nc



subroutine check(status)
	integer, intent (in) :: status
	if(status /= nf90_noerr) then
		print *, trim(nf90_strerror(status))
		stop "nf90 error"
	end if
end subroutine check

end module ncutils
