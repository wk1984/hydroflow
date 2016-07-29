module dem
use ctrl

implicit none

contains


! Calculate the surface slopes.
subroutine slope_calc(c, nx, ny)
	type(controller), intent(inout) :: c
	integer, intent(in) :: nx, ny
	real*4, dimension(nx,ny) :: dzdx, dzdy
	real*4 :: slope
	integer :: i, j
	
	real*4, parameter :: rad2deg = 90.0 / acos(0.0)
	
	! Find dzdx.
	do j=1,ny
		dzdx(1,j) = (c%z(2,j) - c%z(1,j)) / c%dx
		do i=2,nx-1
			dzdx(i,j) = (c%z(i+1,j) - c%z(i-1,j)) / (2.0 * c%dx)
		enddo
		dzdx(nx,j) = (c%z(nx,j) - c%z(nx-1,j)) / c%dx
	enddo

	! Find dzdy.
	do i=1,nx
		dzdy(i,1) = (c%z(i,2) - c%z(i,1)) / c%dy
		do j=2,ny-1
			dzdy(i,j) = (c%z(i,j+1) - c%z(i,j-1)) / (2.0 * c%dy)
		enddo
		dzdy(i,ny) = (c%z(i,ny) - c%z(i,ny-1)) / c%dy
	enddo

	! Calculate the terrain slope in % (actually fraction).
	do i=1,nx
		do j=1,ny
			! Some compilers will not allow dzdx and dzdy to both be 0.0 in
			!   the atan2 computation.
			!         if (abs(dzdx(i,j)).lt.1e-10) dzdx(i,j) = 1e-10
			if (abs(dzdy(i,j)).lt.1e-10) dzdy(i,j) = 1e-10

			! Calculate the slope in degrees.
			slope = rad2deg * atan(sqrt(dzdx(i,j)*dzdx(i,j) + dzdy(i,j)*dzdy(i,j)))
			slope = min(slope,45.0)
			slope = max(slope,5.0)
			slope = slope**0.8

			! Scale these values to be above and below 1.0 for slopes above
			!   and below 15% (see the velocity assumption note below).  The
			!   adjustment below gives a correction of 0.4 for a slope of
			!   5 degrees, a correction of 1.0 for a slope of 15 degrees,
			!   and a correction of 2.4 for a slope of 45 degrees.
			c%slope_corr(i,j) = slope / 15.0**0.8
		enddo
	enddo
end subroutine slope_calc



subroutine watersheds(control, nx, ny, nxny)
	type(controller), intent(inout) :: control
	integer, intent(in) :: nxny,nx,ny
	integer, pointer :: nwatershed_count
	integer, dimension(:,:), pointer :: index_final, nwatershed, nsoln_order
	real, dimension(:,:), pointer :: land, dir
	integer, dimension(nx,ny) :: index2D, nwatershed_bndry
	integer, dimension(nxny,12) :: index_rel
	integer, dimension(nxny,14) :: index_ws_sort, index_ws_sort_2, index_eqn
	integer :: nconnect(nx,ny,9)
	integer, dimension(8) :: i_shift, j_shift
	integer :: i,j,ii,jj,k,kk,maxiter,iter,ninputs,nzeros,nzeros_old,index,ieqn_count,nwshed, &
	i_count,kkk,k_start,k_end,ii_count,i_good,input_eqn,inputs,input_k_pos,n_correct, err
	real d1,d2,d3,d4,d5,d6,d7,d8,undef
	
	index_final      => control%index_eqn
	nwatershed       => control%nwatershed
	nsoln_order      => control%nsoln_order
	land             => control%land
	dir              => control%dir
	nwatershed_count => control%nwatershed_count
	
	
! Define an array which indicates the i and j shift of neighboring
!   connected ode's.
       data i_shift / 1,  1,  1,  0, -1, -1, -1,  0/
       data j_shift / 1,  0, -1, -1, -1,  0,  1,  1/

       undef = -9999.0

! Initialize the connectivity arrays. The first 'k' of this array
!   defines the number of grid cells that feed into this i,j. The
!   rest of the 'ks' are the locations of these feeder cells.
      do j=1,ny
        do i=1,nx
          do k=1,9
            nconnect(i,j,k) = 0
          enddo
        enddo
      enddo

! Zero out the directions along all border cells. This makes it
!   so there will be no inflow from any order cells.
      j = 1
      do i=1,nx
        dir(i,j) = 0.0
      enddo
      j = ny
      do i=1,nx
        dir(i,j) = 0.0
      enddo
      i = 1
      do j=1,ny
        dir(i,j) = 0.0
      enddo
      i = nx
      do j=1,ny
        dir(i,j) = 0.0
      enddo

! Sweep through each grid cell, looking for surrounding grid
!   cells that pour into it. Work on the grid cells one in from
!   the domain boundaries.
      do j=2,ny-1
        do i=2,nx-1
          d1 = dir(i+1,j+1)
          d2 = dir(i+1,j+0)
          d3 = dir(i+1,j-1)
          d4 = dir(i+0,j-1)
          d5 = dir(i-1,j-1)
          d6 = dir(i-1,j+0)
          d7 = dir(i-1,j+1)
          d8 = dir(i+0,j+1)
          if (d1.eq.16.0) then
            nconnect(i,j,1) = 1 + nconnect(i,j,1)
            index = 1 + nconnect(i,j,1)
            nconnect(i,j,index) = 1
          endif
          if (d2.eq.32.0) then
            nconnect(i,j,1) = 1 + nconnect(i,j,1)
            index = 1 + nconnect(i,j,1)
            nconnect(i,j,index) = 2
          endif
          if (d3.eq.64.0) then
            nconnect(i,j,1) = 1 + nconnect(i,j,1)
            index = 1 + nconnect(i,j,1)
            nconnect(i,j,index) = 3
          endif
          if (d4.eq.128.0) then
            nconnect(i,j,1) = 1 + nconnect(i,j,1)
            index = 1 + nconnect(i,j,1)
            nconnect(i,j,index) = 4
          endif
          if (d5.eq.1.0) then
            nconnect(i,j,1) = 1 + nconnect(i,j,1)
            index = 1 + nconnect(i,j,1)
            nconnect(i,j,index) = 5
          endif
          if (d6.eq.2.0) then
            nconnect(i,j,1) = 1 + nconnect(i,j,1)
            index = 1 + nconnect(i,j,1)
            nconnect(i,j,index) = 6
          endif
          if (d7.eq.4.0) then
            nconnect(i,j,1) = 1 + nconnect(i,j,1)
            index = 1 + nconnect(i,j,1)
            nconnect(i,j,index) = 7
          endif
          if (d8.eq.8.0) then
            nconnect(i,j,1) = 1 + nconnect(i,j,1)
            index = 1 + nconnect(i,j,1)
            nconnect(i,j,index) = 8
          endif
        enddo
      enddo

!     do k=1,9
!       print *,'k = ',k
!       do j=ny,1,-1
!         write (*,12) (nconnect(i,j,k),i=1,nx)
!       enddo
!       print *
!     enddo
!     print *

! Populate the watershed numbering array with zeros. Also build
!   an improved direction boundary array.
      do j=1,ny
        do i=1,nx
          nwatershed(i,j) = 0
        enddo
      enddo

! Identify the watershed outlets locations and watershed numbers
!   along the simulation domain border. Again, just work within
!   the border grid cells.
! Modify the outflow grid cells so that the flow is perpendicular to
!   the domain boundaries.
      nwatershed_count = 0

! Deal with the sides, top, and bottom.
      j = 2
      do i=3,nx-2
        if (dir(i,j).eq.4.0 .or. dir(i,j).eq.8.0 .or. &
          dir(i,j).eq.16.0) then
          nwatershed_count = 1 + nwatershed_count
          nwatershed(i,j) = nwatershed_count
          dir(i,j) = 8.0
        endif
      enddo

      j = ny - 1
      do i=3,nx-2
        if (dir(i,j).eq.64.0 .or. dir(i,j).eq.128.0 .or. &
          dir(i,j).eq.1.0) then
          nwatershed_count = 1 + nwatershed_count
          nwatershed(i,j) = nwatershed_count
          dir(i,j) = 128.0
        endif
      enddo

      i = 2
      do j=3,ny-2
        if (dir(i,j).eq.16.0 .or. dir(i,j).eq.32.0 .or. &
          dir(i,j).eq.64.0) then
          nwatershed_count = 1 + nwatershed_count
          nwatershed(i,j) = nwatershed_count
          dir(i,j) = 32.0
        endif
      enddo

      i = nx - 1
      do j=3,ny-2
        if (dir(i,j).eq.1.0 .or. dir(i,j).eq.2.0 .or. &
          dir(i,j).eq.4.0) then
          nwatershed_count = 1 + nwatershed_count
          nwatershed(i,j) = nwatershed_count
          dir(i,j) = 2.0
        endif
      enddo

! Deal with the corners.
      i = 2
      j = 2
      if (dir(i,j).eq.4.0 .or. dir(i,j).eq.8.0 .or. &
        dir(i,j).eq.16.0 .or. dir(i,j).eq.32.0 .or. &
        dir(i,j).eq.64.0) then
        nwatershed_count = 1 + nwatershed_count
        nwatershed(i,j) = nwatershed_count
        dir(i,j) = 16.0
      endif

      i = nx - 1
      j = 2
      if (dir(i,j).eq.1.0 .or. dir(i,j).eq.2.0 .or. &
        dir(i,j).eq.4.0 .or. dir(i,j).eq.8.0 .or. &
        dir(i,j).eq.16.0) then
        nwatershed_count = 1 + nwatershed_count
        nwatershed(i,j) = nwatershed_count
        dir(i,j) = 4.0
      endif

      i = 2
      j = ny - 1
      if (dir(i,j).eq.16.0 .or. dir(i,j).eq.32.0 .or. &
        dir(i,j).eq.64.0 .or. dir(i,j).eq.128.0 .or. &
        dir(i,j).eq.1.0) then
        nwatershed_count = 1 + nwatershed_count
        nwatershed(i,j) = nwatershed_count
        dir(i,j) = 64.0
      endif

      i = nx - 1
      j = ny - 1
      if (dir(i,j).eq.64.0 .or. dir(i,j).eq.128.0 .or. &
        dir(i,j).eq.1.0 .or. dir(i,j).eq.2.0 .or. &
        dir(i,j).eq.4.0) then
        nwatershed_count = 1 + nwatershed_count
        nwatershed(i,j) = nwatershed_count
        dir(i,j) = 1.0
      endif

!     do j=ny,1,-1
!       write (*,13) (dir(i,j),i=1,nx)
!     enddo
!     print *

! Save these as the outer boundary values, to be added back in later.
      do j=1,ny
        do i=1,nx
          nwatershed_bndry(i,j) = 0
        enddo
      enddo

      j = 2
      do i=3,nx-2
        nwatershed_bndry(i,1) = nwatershed(i,j)
      enddo

      j = ny - 1
      do i=3,nx-2
        nwatershed_bndry(i,ny) = nwatershed(i,j)
      enddo

      i = 2
      do j=3,ny-2
        nwatershed_bndry(1,j) = nwatershed(i,j)
      enddo

      i = nx - 1
      do j=3,ny-2
        nwatershed_bndry(nx,j) = nwatershed(i,j)
      enddo

      i = 2
      j = 2
      nwatershed_bndry(1,1) = nwatershed(i,j)

      i = nx - 1
      j = 2
      nwatershed_bndry(nx,1) = nwatershed(i,j)

      i = 2
      j = ny - 1
      nwatershed_bndry(1,ny) = nwatershed(i,j)

      i = nx - 1
      j = ny - 1
      nwatershed_bndry(nx,ny) = nwatershed(i,j)

! Deal with portions of the domain that dump into the ocean.
!   Also set any ocean points (type 24) to have zero direction
!   and set any negative dirs to zero.
      do j=2,ny-1
        do i=2,nx-1
          if (land(i,j)==0) then
            dir(i,j) = 0.0
          endif
          dir(i,j) = max(0.0,dir(i,j))
        enddo
      enddo

      do j=2,ny-1
        do i=2,nx-1
          if (dir(i,j).eq.0.0) then
            d1 = dir(i+1,j+1)
            d2 = dir(i+1,j+0)
            d3 = dir(i+1,j-1)
            d4 = dir(i+0,j-1)
            d5 = dir(i-1,j-1)
            d6 = dir(i-1,j+0)
            d7 = dir(i-1,j+1)
            d8 = dir(i+0,j+1)
            if (d1.eq. 16.0 .or. d2.eq.32.0 .or. d3.eq.64.0 .or. &
                d4.eq.128.0 .or. d5.eq. 1.0 .or. d6.eq. 2.0 .or. &
                d7.eq.  4.0 .or. d8.eq. 8.0) then
              nwatershed_count = 1 + nwatershed_count
              nwatershed(i,j) = nwatershed_count
            endif
          endif
        enddo
      enddo

!     print *
!     print *, nwatershed_count
!     print *

! Distribute the watershed numbers throughout the watersheds.
!   Sweep through looking for values that are greater than 0.
      maxiter = nxny
      nzeros_old = 0
      nzeros = int(1e8)
      do iter=1,maxiter
        nzeros_old = nzeros
        nzeros = 0
        do j=2,ny-1
          do i=2,nx-1
            if (nwatershed(i,j).gt.0) then
              ninputs = nconnect(i,j,1)
              do kk=1,ninputs
                ii = i + i_shift(nconnect(i,j,kk+1))
                jj = j + j_shift(nconnect(i,j,kk+1))
                nwatershed(ii,jj) = nwatershed(i,j)
              enddo
            else
              nzeros = nzeros + 1
            endif
          enddo
        enddo
!       do j=ny,1,-1
!         write (*,12) (nwatershed(i,j),i=1,nx)
!       enddo

! If things are not changing anymore, exit out of this loop.
        if (nzeros.eq.nzeros_old) goto 100
      enddo
  100 continue

! Add the boundary values back in.
      j = 1
      do i=1,nx
        nwatershed(i,j) = nwatershed_bndry(i,j)
      enddo
      j = ny
      do i=1,nx
        nwatershed(i,j) = nwatershed_bndry(i,j)
      enddo
      i = 1
      do j=2,ny-1
        nwatershed(i,j) = nwatershed_bndry(i,j)
      enddo
      i = nx
      do j=2,ny-1
        nwatershed(i,j) = nwatershed_bndry(i,j)
      enddo

!     do j=ny,1,-1
!       write (*,12) (nwatershed(i,j),i=1,nx)
!     enddo
!     print *

! Connect with those boundary points.
      do i=3,nx-2
        if (nwatershed_bndry(i,1).ne.0) then
          nconnect(i,1,1) = 1
          nconnect(i,1,2) = 8
        endif
      enddo
      do i=3,nx-2
        if (nwatershed_bndry(i,ny).ne.0) then
          nconnect(i,ny,1) = 1
          nconnect(i,ny,2) = 4
        endif
      enddo
      do j=3,ny-2
        if (nwatershed_bndry(1,j).ne.0) then
          nconnect(1,j,1) = 1
          nconnect(1,j,2) = 2
        endif
      enddo
      do j=3,ny-2
        if (nwatershed_bndry(nx,j).ne.0) then
          nconnect(nx,j,1) = 1
          nconnect(nx,j,2) = 6
        endif
      enddo

      if (nwatershed_bndry(1,1).ne.0) then
        nconnect(1,1,1) = 1
        nconnect(1,1,2) = 1
      endif
      if (nwatershed_bndry(nx,1).ne.0) then
        nconnect(nx,1,1) = 1
        nconnect(nx,1,2) = 7
      endif
      if (nwatershed_bndry(1,ny).ne.0) then
        nconnect(1,ny,1) = 1
        nconnect(1,ny,2) = 3
      endif
      if (nwatershed_bndry(nx,ny).ne.0) then
        nconnect(nx,ny,1) = 1
        nconnect(nx,ny,2) = 5
      endif

! Merge nconnect(i,j,k) with nwatershed(i,j) and i,j information
!   to build a general indexing array that contains columns of:
!   'i, j, watershed #, # of rivers in, inflow position(1-8)'

! In the first column of the data file is the i location. The
!   second column holds the j location.  The third coulumn has the
!   index 1 if the box is within the watershed and 0 if it is
!   outside the watershed.  The fourth column is a count of the
!   number of grid boxes which drain into the box in question.
!   Columns 5 through 12 (or the next 8) contain an index of the
!   inflowing neighbor cells, where the index of the box
!   immediately above is number 1, and the index number increases
!   in a clockwise manner to number 8, the block above and to the
!   left of the block of interest.  Zeros are placed in unused
!   portions of columns 5 through 12.
      do j=1,ny
        do i=1,nx
          k = i + (j - 1) * nx
          index_rel(k,1) = i
          index_rel(k,2) = j
          index_rel(k,3) = nwatershed(i,j)
          do kk=4,12
            index_rel(k,kk) = nconnect(i,j,kk-3)
          enddo
        enddo
      enddo

!     print *
!     print *, 'i, j, watershed #, # of rivers in, inflow position'
!     do k=1,nxny
!       print 20,(index_rel(k,kk),kk=1,12)
!     enddo
!     print *

! Build a 2-D array that contains an equation number describing
!   each grid cell within the watershed, and therefore its
!   position relative to the neighboring cells.  The equation
!   number must be reset for each watershed.  If the block of
!   interest is not within a numbered watershed, place a zero
!   in its position.
      do nwshed=1,nwatershed_count
        ieqn_count = 0
        do j=1,ny
          do i=1,nx
            k = i + (j - 1) * nx
            if (index_rel(k,3).ne.0) then
              if (index_rel(k,3).eq.nwshed) then
                ieqn_count = ieqn_count + 1
                index2D(i,j) = ieqn_count
              endif
            else
              index2D(i,j) = 0
            endif
          enddo
        enddo
      enddo

!     print *, '2-D matrix indicating ODE number and (i,j) position'
!     do j=ny,1,-1
!       print 12, (index2D(i,j),i=1,nx)
!     enddo
!     print *

! Create a connectivity index in terms of the 2-D ODE index array.
!   Here an array is created which describes the connectivity
!   of the grid cells via the river network, in terms of the
!   ODE equation numbers.  This is done by converting the
!   connectivity information contained in the array 'index', to
!   a form which describes the connectivity in terms of the ODE
!   equation number.  This is accomplished by using the 2-D
!   indexing matrix, index2D.

! Build a index_eqn array, similar to the index_rel array, that
!   can be modified to indicate the ode number of connected ode's.
!   Include space for the ODE equation number and the solution order
!   number.

! The final format for the watershed_index.txt output file will be:
!   i, j, watershed #,equation #,solution order,# of flow inputs,
!   inflow equation # times 8.
      do k=1,nxny
        do kk=1,3
          index_eqn(k,kk) = index_rel(k,kk)
        enddo

        i = index_rel(k,1)
        j = index_rel(k,2)
        index_eqn(k,4) = index2D(i,j)

        index_eqn(k,5) = int(undef)

        do kk=6,14
          index_eqn(k,kk) = index_rel(k,kk-2)
        enddo
      enddo

! Perform the modifications.
      do k=1,nxny
        i = index_eqn(k,1)
        j = index_eqn(k,2)

! Extract the ode index indicating a river connection.
        do kk=1,index_eqn(k,6)
          ii = i + i_shift(index_eqn(k,kk+6))
          jj = j + j_shift(index_eqn(k,kk+6))
          index_eqn(k,kk+6) = index2D(ii,jj)
        enddo
      enddo

! Print the connectinity index array.
!     print *,'connectivity index in terms of ODE number'
!     do k=1,nxny
!       print 20,(index_eqn(k,kk),kk=1,14)
!     enddo
!     print *

! Reorder the index array so that the watersheds are listed
!   in sequential order.
      i_count = 0
      do nwshed=1,nwatershed_count
        do k=1,nxny
          if (index_eqn(k,3).eq.nwshed) then
            i_count = i_count + 1
            do kk=1,14
              index_ws_sort(i_count,kk) = index_eqn(k,kk)
            enddo
          endif
        enddo
      enddo

! If there are any grid cells not in a watershed, place those
!   lines at the end of the sorted file.
      do k=1,nxny
        if (index_eqn(k,3).eq.0) then
          i_count = i_count + 1
          do kk=1,14
            index_ws_sort(i_count,kk) = index_eqn(k,kk)
          enddo
        endif
      enddo

! Determine the order the equations for each watershed set must
!   be solved, from the headwaters to the mouth.  This makes
!   sure all grid cells upstream of the cell in question has
!   been solved before the one in question is solved.
      do nwshed=1,nwatershed_count

! Identify the start and ending positions of the grid cells in
!   this watershed.
        k_start = 0
        k_end = 0
        i_count = 0
        do k=1,nxny
          if (index_ws_sort(k,3).eq.nwshed .and. k_start.eq.0) then
            k_start = k
            i_count = i_count + 1
          elseif (index_ws_sort(k,3).eq.nwshed) then
            k_end = k
            i_count = i_count + 1
          endif
        enddo

!       print *,nwshed,k_start,k_end,i_count

! The first grid cells to be solved are the ones that have no
!   cells flowing into them.
        ii_count = 0
        if (k_start.gt.0) then
			do kkk=k_start,k_end
				if (index_ws_sort(kkk,6).eq.0) then
					ii_count = ii_count + 1
					index_ws_sort(kkk,5) = ii_count
				endif
			enddo
        endif

! Now iteratively loop through the rest of the cells looking for
!   cells that have unresolved inputs.
        do iter=1,i_count

! Search for cells that do not yet have an order assigned to them.
!   For these cells, look to see whether all of the inputs to that
!   cell has been defined. If so, increment the solution order, if
!   not, leave it for the next iteration.
          do kkk=k_start,k_end
            if (index_ws_sort(kkk,5).eq.undef) then
              i_good = 0
              do inputs=1,index_ws_sort(kkk,6)
                input_eqn = index_ws_sort(kkk,6+inputs)
                input_k_pos = k_start + input_eqn - 1
                if (index_ws_sort(input_k_pos,5).ne.undef) then
                  i_good = i_good + 1
                endif
              enddo
              if (i_good.eq.index_ws_sort(kkk,6)) then
                ii_count = ii_count + 1
                index_ws_sort(kkk,5) = ii_count
              endif
            endif
          enddo
        enddo

      enddo

! Sweep through looking for watersheds with only a single grid
!   cell.  I think these will all have 'undef' in column 5.  If
!   found, strip out this watershed and place it at the end of
!   the indexing array.  It's not clear to me why some single-cell
!   watersheds are creeping up, but they cause probems because
!   there are not really a system of equations to solve.
        n_correct = 0
        do k=1,nxny
          if (index_ws_sort(k,5).eq.undef &
            .and. index_ws_sort(k,3).ne.0) then
            n_correct = n_correct + 1
            print *,'found ws with undef solution order; i,j = ', &
              index_ws_sort(k,1),index_ws_sort(k,2)
            do kk=1,2
              index_ws_sort_2(nxny-n_correct+1,kk) = index_ws_sort(k,kk)
            enddo
            do kk=3,4
              index_ws_sort_2(nxny-n_correct+1,kk) = 0
            enddo
            do kk=5,14
              index_ws_sort_2(nxny-n_correct+1,kk) = index_ws_sort(k,kk)
            enddo
          else
            do kk=1,14
              index_ws_sort_2(k-n_correct,kk) = index_ws_sort(k,kk)
            enddo
          endif
        enddo

! Place this new indexing in the working indexing array.
        do k=1,nxny
          do kk=1,14
            index_ws_sort(k,kk) = index_ws_sort_2(k,kk)
          enddo
        enddo

! Before final output, search through the watershed numbers,
!   looking watershed numbers that have no grid cells.  Delete
!   these numbers and renumber the watersheds getting a new total
!   watershed count.
      n_correct = 0
      do nwshed=1,nwatershed_count
        k_start = 0
        k_end = 0
        i_count = 0
        do k=1,nxny
          if (index_ws_sort(k,3).eq.nwshed .and. k_start.eq.0) then
            k_start = k
            i_count = i_count + 1
          elseif (index_ws_sort(k,3).eq.nwshed) then
            k_end = k
            i_count = i_count + 1
          endif
        enddo

        if (i_count.eq.0) then
!         print *,'found watershed with 0 grid cells',i_count,nwshed
          n_correct = n_correct + 1
        endif

        if (k_start.ne.0) then
          do kkk=k_start,k_end
            index_ws_sort(kkk,3) = index_ws_sort(kkk,3) - n_correct
            nwatershed_count = index_ws_sort(kkk,3)
!           print *,kkk,nwatershed_count
          enddo
        endif
      enddo

      print *
      print *, 'final watershed count = ', nwatershed_count
      print *

! Also update the nwatershed array to reflect these new watershed
!   numbers.
      do j=1,ny
        do i=1,nx
          nwatershed(i,j) = 0
          nsoln_order(i,j) = 0
        enddo
      enddo
      do k=1,nxny
        i = index_ws_sort(k,1)
        j = index_ws_sort(k,2)
        nwatershed(i,j) = index_ws_sort(k,3)
        nsoln_order(i,j) = index_ws_sort(k,5)
      enddo

! Now sort the indexing information so that they are listed in
!   order of the solution.  First make a copy of the indexing
!   array that can be modified.
      do k=1,nxny
        do kk=1,14
          index_final(k,kk) = index_ws_sort(k,kk)
        enddo
      enddo

      do nwshed=1,nwatershed_count
        k_start = 0
        k_end = 0
        i_count = 0
        do k=1,nxny
          if (index_ws_sort(k,3).eq.nwshed .and. k_start.eq.0) then
            k_start = k
            i_count = i_count + 1
          elseif (index_ws_sort(k,3).eq.nwshed) then
            k_end = k
            i_count = i_count + 1
          endif
        enddo

        do k=1,i_count
          do kkk=k_start,k_end
            if (index_ws_sort(kkk,5).eq.k) then
              do kk=1,14
                index_final(k_start+k-1,kk) = index_ws_sort(kkk,kk)
              enddo
            endif
          enddo
        enddo
      enddo

      return
      end

end module dem

