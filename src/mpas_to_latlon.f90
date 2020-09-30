! To compile 

! ifort -check bounds -warn all -o mpas_to_latlon src/mpas_to_latlon.f90
! src/mpas_filter_cells_by_area.f src/flip_to_cf.f90 src/mpas_vort_cell.f -I
! src/datetime-fortran/build/include/ -L src/datetime-fortran/build/lib/ -ldatetime

! Original code from Michael Duda duda@ucar.edu
! adapted by Dave Ahijevych ahijevyc@ucar.edu

! 
! This code interpolates fields on an MPAS mesh to a lat-lon grid.
! For fields defined on cells we find use the nearest cell and smooth by X km.
!
! Reads list of fields from standard input.
!  

! Note from Michael:
! For fields defined on vertexes (e.g. vorticity), using the lat-lon edge point, 
! we seek the nearest edge point and use that edge and the 4 edges 
! on that edge to interpolate the velocities.  
! The interpolation uses both the nearest edge and neighboring edges, and uses both 
! the normal and reconstructed-tangential velocities


! Note from ahijevych:
! For a while I had problems with ifort. nlat would change
! from 180 to 14 when I inquired the dimension id of a variable.  Strange.
! Solved when I decreased the size of the dimension for var_output_id and
! var_input_id from 100 to 40, or lower. Started changing allocatable arrays to
! static arrays because I thought that was the problem.  Recompiled in pgf90 and avoided the
! problem. 

program mpas_to_latlon
    use netcdf
    implicit none

    integer, parameter :: RKIND = selected_real_kind(12)
    integer, parameter :: max_lon_dimsize = 7200
    integer, parameter :: max_lat_dimsize = 3600
    integer, parameter :: max_nCells_dimsize = 36864002

    real (kind=RKIND) :: pii

    include 'netcdf.inc'

    integer, external :: iargc

    real (kind=RKIND), dimension(max_nCells_dimsize) :: lat_cell, lon_cell
    real (kind=RKIND), allocatable, dimension(:) :: lat_vertex, lon_vertex, areaCell
    real (kind=RKIND), allocatable, dimension(:,:) :: input_var, input_var_vertices
    real (kind=RKIND), allocatable, dimension(:,:,:) :: output_var, output_var_cf

    integer, allocatable, dimension(:,:) :: verticesOnCell, cellsOnVertex, cellsOnCell
    integer, allocatable, dimension(:,:) :: nearest_vertex_ll
    integer, dimension(max_nCells_dimsize) :: nedgesOnCell
    ! used to allocate(nedgesoncell(nCells))  ! BUG FIX (used to be nEdges)

    character (len=8)  :: now_date
    character (len=10) :: now_time
    character (len=5)  :: now_zone
    character (len=30) :: junkc
    character (len=30), dimension(100) :: interp_var_input, interp_var_output
    integer, dimension(100) :: interp_var_nz, var_output_id, var_input_id
    real (kind=RKIND) :: grid_spacing ! lat-lon mesh spacing in degrees
    real :: filter_radius_km  ! smoothing radius km (zero for no smoothing)
    integer :: nvars, ivar, ierr, ncid, ncells, ncells_id
    integer :: latcell_id, loncell_id, areaCell_id
    integer :: nvertices, nvertices_id, levels_id, nvertlevels
    integer :: verticesOnCell_id, latvertex_id, lonvertex_id
    integer :: time_id, time_axis_id, ntimes
    integer :: i, j, cellsOnVertex_id, maxedges, maxedges_id
    integer :: cellsOnCell_id, nedgesOnCell_id, ndims, natts, xtype
    ! dimids has to be allocated with enough space for at least *ndimsp integers
    ! to be returned - Dave Feb 17, 2013
    integer, dimension(NF_MAX_VAR_DIMS) :: dimids
    integer, dimension(4) :: ncount 
    integer :: ncid_ll, lat_ll_id, lon_ll_id, grid_spacing_id, filter_radius_km_id
    integer :: lat_ll_v_id, lon_ll_v_id

    real (kind=RKIND), dimension(max_lat_dimsize,max_lon_dimsize) :: lat_ll, lon_ll
    real (kind=RKIND), dimension(3,max_lat_dimsize,max_lon_dimsize) :: interp_weights
    integer, dimension(3,max_lat_dimsize,max_lon_dimsize) :: interp_cells
    integer :: nlat, nlon, start_vertex, idimn, iatt
    real (kind=RKIND) :: latitude, longitude, startLon, degreesToRadians, RadiansToDegrees, lat0, lat1
    character (len=256) :: fname, arg, savfile, outfile, tmpdir
    character (len=256) :: meshid, attname
    logical :: file_exists 


    pii = 2.0_RKIND * asin(1.0_RKIND)
    degreesToRadians = pii / 180.0_RKIND
    RadiansToDegrees = 180.0_RKIND / pii

    if (iargc() < 3) then
       write(0,*) 
       write(0,*) 'Usage:'
       write(0,*) 'mpas_to_latlon fname outfile grid_spacing filter_radius_km meshid [lat0 [lat1 [startLon]]]'
       write(0,*) 
       write(0,*) ' fname            : input MPAS filename'
       write(0,*) ' outfile          : output filename'
       write(0,*) ' grid_spacing     : output grid spacing in degrees'
       write(0,*) ' filter_radius_km : radius of circular smoothing filter in km'
       write(0,*) ' meshid           : string to identify MPAS mesh (uni, wp, us, etc.)'
       write(0,*) ' lat0             : minimum latitude in deg (default -90)'
       write(0,*) ' lat1             : maximum latitude in deg (default +90)'
       write(0,*) ' startLon         : starting longitude for output grid (default -180)'
       write(0,*) 
       write(0,*) 'Example:'
       write(0,*) 'mpas_to_latlon diagnostics.2013-09-01_00:00:00.nc latlon_0.500deg_025km/diagnostics.2013-09-01_00:00:00.nc 0.5 25 mpas3 -5. 50. -180.'
       write(0,*) 
       call exit(1) 
    end if
    call getarg(1,fname)
    call getarg(2,outfile)

    if (fname == outfile) then
        write(0,*)' mpas_to_latlon: infile equals outfile. Exiting'
        call exit(1)
    end if
    call getarg(3,arg)
    read(arg,*) grid_spacing
    call getarg(4,arg)
    read(arg,*) filter_radius_km
    call getarg(5,arg)
    read(arg,*) meshid
    lat0 = -90.0
    lat1 =  90.0
    startLon = -180.
    call getarg(6,arg)
    if (len_trim(arg).ne.0) read(arg,*) lat0 
    call getarg(7,arg)
    if (len_trim(arg).ne.0) read(arg,*) lat1
    call getarg(8,arg)
    if (len_trim(arg).ne.0) read(arg,*) startLon


    !
    !  create lat lon grid dimensions
    !

    nlat = ceiling((lat1-lat0) / grid_spacing)
    nlon = ceiling(360. / grid_spacing)
    if(nlat>max_lat_dimsize.or.nlon>max_lon_dimsize)then
       write(0,*)
       write(0,*)' mpas_to_latlon: either max_lat_dimsize or max_lon_dimsize is too small. &
           recompile after making it bigger.'
       write(0,*)
       call exit(1)
    endif


    !
    ! Populate list of variables to interpolate
    !

    ! read list of 2D vars from standard input
    write(0,*) ' mpas_to_latlon: reading from standard input (use ctrl-d to end manual entry) '
    read(5,*,iostat=ierr) junkc
    ivar = 1
    do while(ierr == 0)
       interp_var_input(ivar) = junkc 
       interp_var_output(ivar) = interp_var_input(ivar)
       ivar = ivar+1
       read(5,*,iostat=ierr) junkc
    end do
    nvars = ivar-1




    !
    !  open MPAS file and read description of MPAS mesh
    !

    ierr = nf_open(trim(fname), NF_NOWRITE, ncid)
    write(0,*)'opened '//trim(fname)//' ncid=',ncid
    if (ierr.ne.NF_NOERR) then
        write(0,*) ' nf_open error opening mpas file '//trim(fname)
        call handle_err(ierr)
    end if

    !
    ! Input must have dimensions 'nCells', 'nVertices', and 'Time'.
    !

    ierr = nf_inq_dimid(ncid, 'nCells', ncells_id)
    if(ierr.ne.NF_NOERR)then
        write(0,*) ' nf_inq_dimid err for ncells'
        call handle_err(ierr)
    end if

    ierr = nf_inq_dim(ncid,ncells_id,junkc,ncells)
    if(ierr.ne.NF_NOERR)then
        write(0,*) ' nf_inq_dim err for ncells'
        write(0,*) ' ncells for mesh ',ncells
        call handle_err(ierr)
    end if
    if(ncells.gt.max_nCells_dimsize)then
        write(0,*) 'ncells ',ncells,' exceeds limit ',max_nCells_dimsize
        write(0,*) 'increase max_nCells_dimsize parameter and recompile'
        stop
    end if

    ierr = nf_inq_dimid(ncid,'nVertices',nVertices_id)
    if(ierr.ne.NF_NOERR)then
        write(0,*) ' nf_inq_dimid err for nVertices id'
        call handle_err(ierr)
    end if

    ierr = nf_inq_dim(ncid,nVertices_id,junkc,nVertices)
    if(ierr.ne.NF_NOERR)then
        write(0,*) ' nf_inq_dim err for nVertices'
        write(0,*) ' nVertices for mesh ',nVertices
        call handle_err(ierr)
    end if

    ierr = nf_inq_dimid(ncid, 'Time', time_id)
    if(ierr.ne.NF_NOERR)then
        write(0,*) ' nf_inq_dimid err for Time id'
        call handle_err(ierr)
    end if

    ierr = nf_inq_dim(ncid, time_id, junkc, ntimes)
    if(ierr.ne.NF_NOERR.or.ntimes.ne.1)then
        write(0,*) ' err for times'
        write(0,*) ' number of times ', ntimes
        call handle_err(ierr)
    end if

    !
    ! create output file
    !
    ierr = nf_create(outfile,NF_64BIT_OFFSET, ncid_ll)
    if (ierr.ne.NF_NOERR) then 
        write(0,*) ' err for nf_create '//trim(outfile)
        write(0,*) ' do you need to create a subdirectory first?'
        call handle_err(ierr)
    end if

    ierr = nf_def_dim( ncid_ll, 'time', NF_UNLIMITED, time_id )
    if(ierr.ne.NF_NOERR)write(0,*) ' err on time_id dim def'
    call handle_err(ierr)

    ierr = nf_def_dim( ncid_ll, 'lat', nlat, lat_ll_id )
    if (ierr.ne.NF_NOERR) write(0,*) ' err defining lat dim'
    call handle_err(ierr)

    ierr = nf_def_dim( ncid_ll, 'lon', nlon, lon_ll_id )
    if (ierr.ne.NF_NOERR) write(0,*) ' err defining lon dim'
    call handle_err(ierr)

    ncount(1) = time_id
    ncount(2) = 0
    ncount(3) = 0
    ncount(4) = 0
    ierr = nf_def_var( ncid_ll, 'time', NF_DOUBLE, 1, ncount, time_axis_id )
    call handle_err(ierr)
    ierr = nf_put_att_text(ncid_ll, time_axis_id, 'long_name', len('time'),      'time')
    call handle_err(ierr)
    write(junkc,'(a11,a)') 'days since ', trim('1970-01-01 00:00:00')
    ierr = nf_put_att_text(ncid_ll, time_axis_id, 'units',     len_trim(junkc),  trim(junkc))
    call handle_err(ierr)
    ierr = nf_put_att_text(ncid_ll, time_axis_id, 'calendar',  len('gregorian'), 'gregorian')
    call handle_err(ierr)

    ncount(1) = lat_ll_id
    ncount(2) = 0
    ncount(3) = 0
    ncount(4) = 0
    ierr = nf_def_var( ncid_ll, 'lat', NF_DOUBLE, 1, ncount, lat_ll_v_id )
    if (ierr.ne.NF_NOERR) write(0,*) ' err on lat_ll_id dim def'
    call handle_err(ierr)
    ierr = nf_put_att_text(ncid_ll, lat_ll_v_id, 'long_name', len('latitude'),      'latitude')
    call handle_err(ierr)
    ierr = nf_put_att_text(ncid_ll, lat_ll_v_id, 'units',     len('degrees_north'), 'degrees_north')
    call handle_err(ierr)
    ierr = nf_put_att(ncid_ll, lat_ll_v_id, 'valid_min', NF_DOUBLE, 1, lat0)
    call handle_err(ierr)
    ierr = nf_put_att(ncid_ll, lat_ll_v_id, 'valid_max', NF_DOUBLE, 1, lat1)
    call handle_err(ierr)

    ncount(1) = lon_ll_id
    ncount(2) = 0
    ncount(3) = 0
    ncount(4) = 0
    ierr = nf_def_var( ncid_ll, 'lon', NF_DOUBLE, 1, ncount, lon_ll_v_id )
    if (ierr.ne.NF_NOERR) write(0,*) ' err on lon_ll_id dim def'
    call handle_err(ierr)
    ierr = nf_put_att_text(ncid_ll, lon_ll_v_id, 'long_name', len('longitude'),    'longitude')
    call handle_err(ierr)
    ierr = nf_put_att_text(ncid_ll, lon_ll_v_id, 'units',     len('degrees_east'), 'degrees_east')
    call handle_err(ierr)


    ierr = nf_def_var( ncid_ll, 'mesh_spacing', NF_DOUBLE, 0, 0, grid_spacing_id )
    call handle_err(ierr)
    ierr = nf_put_att_text(ncid_ll, grid_spacing_id, 'long_name', len('lat lon grid spacing (command line argument passed to mpas_to_latlon)'), &
                                                                      'lat lon grid spacing (command line argument passed to mpas_to_latlon)')
    call handle_err(ierr)
    ierr = nf_put_att_text(ncid_ll, grid_spacing_id, 'units',     len('degrees'),    'degrees')
    call handle_err(ierr)
    ierr = nf_def_var( ncid_ll, 'filter_radius_km', NF_DOUBLE,    0, 0, filter_radius_km_id )
    call handle_err(ierr)
    ierr = nf_put_att_text(ncid_ll, filter_radius_km_id, 'long_name', len('radius of influence for mpas_to_latlon smoothing filter'), &
                                                                     'radius of influence for mpas_to_latlon smoothing filter')
    call handle_err(ierr)
    ierr = nf_put_att_text(ncid_ll, filter_radius_km_id, 'units',     len('km'),    'km')
    call handle_err(ierr)


    do i=1,nvars
        ierr = nf_inq_varid(ncid,interp_var_input(i),var_input_id(i))
        call handle_err(ierr)
        write(0,*) 'got var_input_id for '//TRIM(interp_var_input(i)), var_input_id(i)
        if (ierr.ne.NF_NOERR) then 
            write(0,*) ' ierr ',ierr,' for input var ',interp_var_input(i),var_input_id(i)
            call handle_err(ierr)
        end if
        ierr = nf_inq_var(ncid, var_input_id(i), junkc, xtype, ndims, dimids, natts)
        if (ierr.ne.NF_NOERR) then 
            write(0,*) ' ierr ',ierr,' for nf_inq_var ',interp_var_input(i),var_input_id(i)
            call handle_err(ierr)
        end if


        ! TODO: pass through 1-D coordinate variables without interpolating
        ncount = -1
        ncount(1) = lon_ll_id
        ncount(2) = lat_ll_id
        ! You add 1 dimension by converting nCells or nVertices to (lat, lon). 
        ncount(ndims+1) = time_id
        interp_var_nz(i) = 1
        ! See if any dimensions are vertical. If so, define dimension in output file.
        
        do idimn=1,ndims 
           ierr = nf_inq_dim(ncid,dimids(idimn),junkc,nVertLevels)
           if (ierr.ne.NF_NOERR) write(0,*) ' err inq dimn ',idimn,junkc
           call handle_err(ierr)
           if (trim(junkc) .eq. 'nVertLevels' .or. junkc(1:10) .eq. 'nIsoLevels' .or. junkc(2:12) == '_iso_levels') then
              write(0,*) ' found input vertical dimension '//TRIM(junkc)
              if (ndims < 3) then 
                 write(0,*) '  but ndims=',ndims
                 stop ! sanity check
              end if
              ierr = nf_inq_dimid(ncid_ll, junkc, levels_id)
              if(ierr.ne.NF_NOERR) then
                ierr = nf_def_dim(ncid_ll, junkc, nVertLevels, levels_id)
                write(0,*) ' define output vertical dimension '//TRIM(junkc)
                if (ierr.ne.NF_NOERR) write(0,*) ' err defining vertical dimension for output',idimn,junkc
                call handle_err(ierr)
              end if 
              ncount(3) = levels_id
              interp_var_nz(i) = nVertLevels
           else
              write(0,*) ' input dimension '//TRIM(junkc)//' not vertical'
           end if
        end do
        if (ndims >= 3 .and. ncount(3) == -1) then
           write(0,'(A,I1,A)') ' no vertical dimension found in ',ndims,'-D input. stopping.'
           stop ! sanity check
        end if
        ! Definine output variable. Dimensions are longitude, latitude, and
        ! time. If vertical dimension present, dim 3 is vertical and dim 4
        ! is time.
        !write(0,*)'defining '//interp_var_output(i)
        !write(0,*)'ndims+1,ncount= ',ndims+1,ncount
        ! Used to force NF_DOUBLE - changed to xtype (the type of the input
        ! variable)
        ierr = nf_def_var( ncid_ll, interp_var_output(i), xtype, ndims+1, ncount, var_output_id(i) )
        if(ierr.ne.NF_NOERR) write(0,*) ' output var ',interp_var_output(i),' id ',var_output_id(i)
        call handle_err(ierr)

        ! Copy attributes from input variable to output variable
        do iatt=1,natts
            ierr = nf_inq_attname(ncid, var_input_id(i), iatt, attname)
            if (ierr.ne.NF_NOERR) then
                write(0,*) ' err for nf_inq_attname '//interp_var_input(i),iatt
                call handle_err(ierr)
            end if
            ierr = nf_copy_att(ncid, var_input_id(i), attname, ncid_ll, var_output_id(i))
            if (ierr.ne.NF_NOERR) then
                write(0,*) ' err for nf_copy_att ',attname
                call handle_err(ierr)
            end if
        end do
    end do


    ierr = nf_put_att_text(ncid_ll, NF_GLOBAL, 'Conventions', len('CF-1.0'), 'CF-1.0')
    call handle_err(ierr)
    ierr = nf_put_att_text(ncid_ll, NF_GLOBAL, 'model', len('mpas'), 'mpas')
    call handle_err(ierr)
    ierr = nf_put_att_text(ncid_ll, NF_GLOBAL, 'grid', len('interp_latlon'), 'interp_latlon')
    call handle_err(ierr)
    ierr = nf_put_att_text(ncid_ll, NF_GLOBAL, 'native_grid', len('voronoi'), 'voronoi')
    call handle_err(ierr)
    ierr = nf_put_att_text(ncid_ll, NF_GLOBAL, 'original_filename', len(trim(fname)), TRIM(fname))
    call handle_err(ierr)
    call date_and_time(DATE=now_date, TIME = now_time, ZONE = now_zone)
    ierr = nf_put_att_text(ncid_ll, NF_GLOBAL, 'run_date_and_time', len('mpas_to_latlon: '//now_date//' '//now_time//' '//now_zone), &
                                                                        'mpas_to_latlon: '//now_date//' '//now_time//' '//now_zone)
    call handle_err(ierr)

    ! Copy global attributes to output file
    ierr = nf90_inquire(ncid, nAttributes=natts) 
    if (ierr.ne.NF_NOERR) call handle_err(ierr)
    do iatt=1,natts
        ierr = nf_inq_attname(ncid, NF_GLOBAL, iatt, attname)
        if (ierr.ne.NF_NOERR) then
            write(0,*) ' err for nf_inq_attname global att ',iatt
            call handle_err(ierr)
        end if
        ierr = nf_copy_att(ncid, NF_GLOBAL, attname, ncid_ll, NF_GLOBAL)
        if (ierr.ne.NF_NOERR) then
            write(0,*) ' err for nf_copy_att ',attname
            call handle_err(ierr)
        end if
    end do

    ! netcdf - exiting define mode for output file
    ierr = nf_enddef(ncid_ll)
    if(ierr.ne.NF_NOERR)write(0,*) ' err for nf_enddef(ncid_ll)'
    call handle_err(ierr)

    !
    !  create lat lon mesh for interp
    !
    do j=1, nlon
        longitude = startLon + 0.5*grid_spacing + float(j-1)*grid_spacing
        do i=1, nlat
            latitude = (lat0+0.5*grid_spacing) + float(i-1)*grid_spacing
            lat_ll(i,j) = latitude*degreesToRadians
            lon_ll(i,j) = longitude*degreesToRadians
        end do
    end do


    ! get mpas->latlon geometry
    ! use a save file if it exists
    ! Assumed triangulation may be uniquely identified by nCells and lat range of output.
    ! This assumption burned me when I tried the mpas3 mesh, which is simply the
    ! mpas2 mesh shifted around so the fine part is over the W Pacific.
    ! Add a mesh ID (mpas, mpas2, or mpas3)
    ! Finally added a mesh ID and recompiled on Jun 5, 2014.
    call get_environment_variable("TMPDIR", tmpdir)
    ! Assumes $TMPDIR 
    ! available to store savefile. Perhaps replace with home directory or /tmp
    ! directory for other institutions to use.
    ! Tried to zero-pad lat and lon so we don't get spaces in the filename, but
    ! fortran only zero-pads integers.
    ! SP turns on leading + or - for the rest of the format string. 
    write(fmt='(A,"/",A,"_",I8.8,"_",F5.3,"deg",SP,F7.3,"N",F7.3,"N",F8.3,"E")', unit=savfile) &
      TRIM(tmpdir),trim(meshid), nCells, grid_spacing, lat0, lat1, startLon
    write(0,'(A)')'looking for save file "'//trim(savfile)//'"'
    inquire(file=savfile, exist=file_exists)
    if (file_exists) then 
        write(0,*) ' mpas_to_latlon: using save file '//trim(savfile)
        open(unit=2, file=savfile, form='unformatted')
        read(2) maxEdges, interp_weights, interp_cells
        !write(0,*) maxEdges, nCells, nVertices
        !write(0,*) 'interp_weights=',interp_weights, 'interp_cells=',interp_cells
        read(2) nEdgesOnCell
        allocate(verticesOnCell(maxEdges,nCells))
        read(2) verticesOnCell
        allocate(cellsOnCell(maxEdges,nCells))
        read(2) cellsOnCell
        allocate(areaCell(nCells))
        read(2) areaCell
        close(2)
    else

        write(0,'(A)')'did not find save file'
        ierr = nf_inq_dimid(ncid, 'maxEdges', maxedges_id)
        if (ierr.ne.NF_NOERR) then
            write(0,*) ' err for maxedges id ', NF_STRERROR(ierr)
            write(0,*) ' process init.nc first, so maxEdges and nEdges are defined in save file'
            call exit(2)
        end if

        ierr = nf_inq_dim(ncid, maxedges_id, junkc, maxedges)
        if (ierr.ne.NF_NOERR) then
            write(0,*) ' err for nf_inq_dim maxEdges ', NF_STRERROR(ierr)
            call handle_err(ierr)
            call exit(2)
        end if
        write(0,*) ' maxedges for mesh ', maxedges

        allocate(lat_vertex(nvertices))
        allocate(lon_vertex(nvertices))

        !allocate(nedgesoncell(nCells))  ! BUG FIX (used to be nEdges)
        allocate(verticesOnCell(maxedges,ncells))
        allocate(cellsOnCell(maxedges,ncells))
        allocate(cellsOnVertex(3,nVertices))
        allocate(areaCell(nCells))

        ierr = nf_inq_varid(ncid,'latCell',latcell_id)
        call handle_err(ierr)
        ierr = nf_inq_varid(ncid,'lonCell',loncell_id)
        call handle_err(ierr)
        ierr = nf_get_var_double(ncid,latcell_id,lat_cell)
        call handle_err(ierr)
        ierr = nf_get_var_double(ncid,loncell_id,lon_cell)
        call handle_err(ierr)

        ierr = nf_inq_varid(ncid,'latVertex',latVertex_id)
        call handle_err(ierr)
        ierr = nf_inq_varid(ncid,'lonVertex',lonVertex_id)
        call handle_err(ierr)
        ierr = nf_get_var_double(ncid,latVertex_id,lat_Vertex)
        call handle_err(ierr)
        ierr = nf_get_var_double(ncid,lonVertex_id,lon_Vertex)
        call handle_err(ierr)

        ierr = nf_inq_varid(ncid,'verticesOnCell',verticesOnCell_id)
        call handle_err(ierr)
        ierr = nf_get_var_int(ncid,verticesOnCell_id,verticesOnCell)
        call handle_err(ierr)

        ierr = nf_inq_varid(ncid,'cellsOnVertex',cellsOnVertex_id)
        call handle_err(ierr)
        ierr = nf_get_var_int(ncid,cellsOnVertex_id,cellsOnVertex)
        call handle_err(ierr)

        ierr = nf_inq_varid(ncid,'nEdgesOnCell',nedgesOnCell_id)
        call handle_err(ierr)
        ierr = nf_get_var_int(ncid,nedgesoncell_id,nEdgesOnCell)
        call handle_err(ierr)

        ierr = nf_inq_varid(ncid,'cellsOnCell',cellsOnCell_id)
        call handle_err(ierr)
        ierr = nf_get_var_int(ncid,cellsoncell_id,cellsOnCell)
        call handle_err(ierr)

        ierr = nf_inq_varid(ncid,'areaCell',areaCell_id)
        call handle_err(ierr)
        ierr = nf_get_var_double(ncid,areaCell_id,areaCell)
        call handle_err(ierr)


        !
        ! find nearest vertex for interpolation
        !
        allocate(nearest_vertex_ll(nlat,nlon))
        start_vertex = 1
        do j=1,nlon
            do i=1,nlat
                !write(0,*) 'find nearest vertex to',j,i
                nearest_vertex_ll(i,j) = nearest_vertex( lat_ll(i,j), lon_ll(i,j), start_vertex, ncells, nvertices, maxedges, &
                                                         nedgesOnCell, verticesOnCell, cellsOnVertex, &
                                                         lat_Cell, lon_Cell, lat_vertex, lon_vertex )
                !if (nearest_vertex_ll(i,j) .eq. 0) write(0,*), 'after nearest_vertex() the nearest vertex=0 when j=',j,' i=',i 
                start_vertex = nearest_vertex_ll(i,j)        
            end do
            start_vertex = nearest_vertex_ll(1,j)
        end do

        !
        !  compute interp weights for all the points on the lat-lon mesh
        !

        write(0,*) ' computing interp weights '
        call compute_interp_weights_all ( lat_ll, lon_ll, nearest_vertex_ll, nCells,           &
                                          nVertices, cellsOnVertex, lat_cell, lon_cell,        &
                                          interp_weights, interp_cells, nlat, nlon )

        write(0,*) ' interp weights complete '

        ! write expensive variables to save file, which can be used instead of
        ! rerunning this section of code.
        write(0,'(A)')'writing to save file'
        open(unit=2, file=savfile, form='unformatted')
        write(2) maxEdges, interp_weights, interp_cells
        write(2) nEdgesOnCell
        write(2) verticesOnCell
        write(2) cellsOnCell
        write(2) areaCell
        close(2)
        write(0,'(A)')'wrote expensive variables to save file "'//trim(savfile)//'"'
    end if

    ierr = nf_put_vara_double(ncid_ll, time_axis_id, 1, 1, days_since_1970(trim(fname)) )
    call handle_err(ierr)
    ierr = nf_put_var_double(ncid_ll, lat_ll_v_id, lat_ll(:,1)*RadiansToDegrees )
    call handle_err(ierr)
    ierr = nf_put_var_double(ncid_ll, lon_ll_v_id, lon_ll(1,:)*RadiansToDegrees )
    call handle_err(ierr)
    ierr = nf_put_var_double(ncid_ll, grid_spacing_id, grid_spacing)
    call handle_err(ierr)
    ierr = nf_put_var_double(ncid_ll, filter_radius_km_id, filter_radius_km)
    call handle_err(ierr)
 
    ! for each variable - read, interp, and write
    do i=1,nvars
        nVertLevels = interp_var_nz(i)
        allocate(input_var(nVertLevels,nCells))
        allocate(output_var(nVertLevels,nlat,nlon))
        allocate(output_var_cf(nlon,nlat,nVertLevels))

        ierr = nf_inq_var(ncid, var_input_id(i), junkc, xtype, ndims, dimids, natts)
        call handle_err(ierr)
  
        !write(0,*) '1 nlon=',nlon,' nlat=',nlat,' ncid=',ncid,' var_input_id(i)=',var_input_id(i)
        ierr = nf_inq_vardimid(ncid, var_input_id(i), dimids)
        call handle_err(ierr)
        !write(0,*) '2 dimids=',dimids, ' nVertices_id=',nVertices_id,' nlon=',nlon,' nlat=',nlat, ierr
        if (dimids(1) .eq. nVertices_id) then 
            write(0,*) trim(interp_var_output(i))//' is a vertex variable with ',nVertices,'vertices'
            allocate(input_var_vertices(nVertLevels,nVertices))
            ierr = nf_get_var_double( ncid,var_input_id(i),input_var_vertices)
            if(ierr.ne.NF_NOERR) then
               write(0,*)'could not read '//trim(interp_var_output(i))//' ierr=',ierr
               write(0,*)'ncid=',ncid,' var_input_id(i)=',var_input_id(i)
               write(0,*)NF_STRERROR(ierr)
               call exit(ierr)
            end if

            ! INPUT: input_var_vertices OUTPUT: input_var (on cells)
            !write(0,*)'calling mpas_vort_cell1'
            call mpas_vort_cell1(nEdgesOnCell, verticesOnCell, maxEdges, nVertLevels, nCells, nVertices, input_var_vertices, input_var)


            deallocate(input_var_vertices)
        else
            ierr = nf_get_var_double( ncid,var_input_id(i),input_var)
            call handle_err(ierr)
        end if

        ! Filter is based on mpas_filter_cells.f, but the number of
        ! smoothing passes is inversely proportional to sqrt(areaCell)
        ! filter_radius_km is in km.
        !write(0,*)'calling mpas_filter_cells_by_area. nEdgesOnCell(1)=',nEdgesOnCell(1), 'cellsOnCell(:,1)=',cellsOnCell(:,1),'maxEdges=',maxEdges,'nCells=',nCells
        !write(0,*)' areaCell(1)=',areaCell(1), ' filter_radius_km=', filter_radius_km
        call mpas_filter_cells_by_area(nEdgesOnCell, cellsOnCell, areaCell, maxEdges, nVertLevels, nCells, input_var, filter_radius_km)
        !write(0,*)'after mpas_filter_cells_by_area: minval(cellsOnCell)=',minval(cellsOnCell)
        !write(0,*)'after mpas_filter_cells_by_area: maxval(cellsOnCell)=',maxval(cellsOnCell)


        
        call interp_delaunay ( interp_weights, interp_cells,      &
                               input_var, output_var,             &
                               nVertLevels, nVertLevels, &
                               nCells, nlat, nlon                )

        ! lat_ll is still in radians here.
        !call flip_to_cf_and_circle_smooth( output_var, output_var_cf, nlon, nlat, nVertLevels, lat_ll(:,1), filter_radius_km)
        call flip_to_cf( output_var, output_var_cf, nlon, nlat, nVertLevels)

        ierr = nf_put_var_double(ncid_ll, var_output_id(i), output_var_cf )
        call handle_err(ierr)
   
        write(0,'(A,2E15.6,A,A)') '  mpas_to_latlon: max/min',maxval(output_var_cf), minval(output_var_cf), ' for ',interp_var_output(i)

        deallocate(input_var)
        deallocate(output_var)
        deallocate(output_var_cf)

    end do

    call handle_err( nf_close(ncid_ll) )
    write(0,*) ' mpas_to_latlon: created '//trim(outfile)
    call handle_err( nf_close(ncid) )



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 



    subroutine interp_delaunay( interp_weights, interp_cells,    &
                                value_cell, value_interp,        &
                                nVertLevels_cells, nVertlevels_ll, nCells, nlat, nlon )

        implicit none

        integer, intent(in) :: nCells, nVertLevels_cells, nVertlevels_ll, nlat, nlon
        real (kind=RKIND), dimension(nVertLevels_cells, nCells), intent(in) :: value_cell
        real (kind=RKIND), dimension(nVertLevels_ll, nlat, nlon), intent(out) :: value_interp
        real (kind=RKIND), dimension(3, nlat, nlon), intent(in) :: interp_weights
        integer, dimension(3, nlat, nlon), intent(in) :: interp_cells

        integer :: i, j, k

        do j=1,nlon
        do i=1,nlat
        do k=1,min(nVertlevels_ll, nVertlevels_cells)
           value_interp(k,i,j) =   interp_weights(1,i,j)*value_cell(k,interp_cells(1,i,j))  &
                                 + interp_weights(2,i,j)*value_cell(k,interp_cells(2,i,j))  &
                                 + interp_weights(3,i,j)*value_cell(k,interp_cells(3,i,j))
        end do
        end do
        end do

    end subroutine interp_delaunay

    ! moved subroutine flip_to_cf_and_circle_smooth to its own file. Mar 12, 2014


    subroutine compute_interp_weights_all( target_lat, target_lon, nearest_vtx,            &
                                        nCells, nVertices, cellsOnVertex,               &
                                        lat_cell, lon_cell,                             &
                                        interp_weights, interp_cells, nlat, nlon )

        implicit none

        integer, intent(in) :: nlat, nlon
        integer, dimension(nlat,nlon), intent(in) :: nearest_vtx
        real (kind=RKIND), dimension(nlat,nlon), intent(in) :: target_lat, target_lon
        integer, intent(in) :: nCells, nVertices
        integer, dimension(3,nVertices), intent(in) :: cellsOnVertex
        real (kind=RKIND), dimension(nCells), intent(in) :: lat_cell, lon_cell
        real (kind=RKIND), intent(out), dimension(3,nlat,nlon) :: interp_weights
        integer, intent(out), dimension(3,nlat,nlon) :: interp_cells
    
        real (kind=RKIND), dimension(3) :: weights
     
        integer :: i, j, k
        real (kind=RKIND), parameter :: eps = 1.e-010_RKIND
        real (kind=RKIND), parameter :: invpower = 1.0_RKIND
     
        real (kind=RKIND) :: lat0, lon0, lat1, lon1, lat2, lon2, lat3, lon3 
        real (kind=RKIND) :: area123, area012, area023, area013, sum_weights
     
        do j=1,nlon
        do i=1,nlat
        lat0 = target_lat(i,j)
        lon0 = target_lon(i,j)
        !if (j.eq.2) write(0,*) j,i,lon0,lat0,nearest_vtx(i,j)
        lat1 = lat_cell(cellsOnVertex(1,nearest_vtx(i,j)))
        lon1 = lon_cell(cellsOnVertex(1,nearest_vtx(i,j)))
        lat2 = lat_cell(cellsOnVertex(2,nearest_vtx(i,j)))
        lon2 = lon_cell(cellsOnVertex(2,nearest_vtx(i,j)))
        lat3 = lat_cell(cellsOnVertex(3,nearest_vtx(i,j)))
        lon3 = lon_cell(cellsOnVertex(3,nearest_vtx(i,j)))
     
        if (lon0 < 0) lon0 = lon0+2.*pii
        if (lon1 < 0) lon1 = lon1+2.*pii
        if (lon2 < 0) lon2 = lon2+2.*pii
        if (lon3 < 0) lon3 = lon3+2.*pii
     
        area123 = triangle_area(lat1,lon1,lat2,lon2,lat3,lon3,1.0_RKIND)
        area012 = triangle_area(lat0,lon0,lat1,lon1,lat2,lon2,1.0_RKIND)
        area023 = triangle_area(lat0,lon0,lat2,lon2,lat3,lon3,1.0_RKIND)
        area013 = triangle_area(lat0,lon0,lat1,lon1,lat3,lon3,1.0_RKIND)
    
        !
        !  check areas
        !
        weights(1) = area023/area123
        weights(2) = area013/area123
        weights(3) = area012/area123
    
        sum_weights = weights(1)+weights(2)+weights(3)
    
        weights(1) = weights(1)/sum_weights
        weights(2) = weights(2)/sum_weights
        weights(3) = weights(3)/sum_weights
    
        do k=1,3
            interp_cells(k,i,j) = cellsOnVertex(k,nearest_vtx(i,j))
            interp_weights(k,i,j) = weights(k)
        end do
    
        end do
        end do

    end subroutine compute_interp_weights_all


    subroutine convert_lx(x,y,z,radius,lat,lon)

        implicit none

        real (kind=RKIND) x,y,z,lat,lon,radius

        z = radius*sin(lat)
        x = radius*cos(lon)*cos(lat)
        y = radius*sin(lon)*cos(lat)

    end subroutine convert_lx


    integer function nearest_vertex( target_lat, target_lon, &
                                     start_vertex, &
                                     nCells, nVertices, maxEdges, &
                                     nEdgesOnCell, verticesOnCell, &
                                     cellsOnVertex, latCell, lonCell, &
                                     latVertex, lonVertex )

        implicit none
    
        real (kind=RKIND), intent(in) :: target_lat, target_lon
        integer, intent(in) :: start_vertex
        integer, intent(in) :: nCells, nVertices, maxEdges
        integer, dimension(nCells), intent(in) :: nEdgesOnCell
        integer, dimension(maxEdges,nCells), intent(in) :: verticesOnCell
        integer, dimension(3,nVertices), intent(in) :: cellsOnVertex
        real (kind=RKIND), dimension(nCells), intent(in) :: latCell, lonCell
        real (kind=RKIND), dimension(nVertices), intent(in) :: latVertex, lonVertex
    
    
        integer :: i, cell1, cell2, cell3, iCell
        integer :: iVtx
        integer :: current_vertex 
        real (kind=RKIND) :: cell1_dist, cell2_dist, cell3_dist
        real (kind=RKIND) :: current_distance, d
        real (kind=RKIND) :: nearest_distance
       
        nearest_vertex = start_vertex
        current_vertex = -1
        !write(0,*) 'in nearest_vertex(), target_lat=',target_lat,' target_lon=',target_lon 
        do while (nearest_vertex /= current_vertex)
            !write(0,*) 'in nearest_vertex(), nearest_vertex=',nearest_vertex,' current_vertex=',current_vertex 
            current_vertex = nearest_vertex
            current_distance = sphere_distance(latVertex(current_vertex), lonVertex(current_vertex), &
                                               target_lat,                target_lon,                1.0_RKIND)
            nearest_vertex = current_vertex
            nearest_distance = current_distance
            cell1 = cellsOnVertex(1,current_vertex) 
            cell2 = cellsOnVertex(2,current_vertex) 
            cell3 = cellsOnVertex(3,current_vertex) 
            cell1_dist = sphere_distance(latCell(cell1), lonCell(cell1), target_lat, target_lon, 1.0_RKIND)
            cell2_dist = sphere_distance(latCell(cell2), lonCell(cell2), target_lat, target_lon, 1.0_RKIND)
            cell3_dist = sphere_distance(latCell(cell3), lonCell(cell3), target_lat, target_lon, 1.0_RKIND)
            if (cell1_dist < cell2_dist) then
                if (cell1_dist < cell3_dist) then
                    iCell = cell1
                else
                    iCell = cell3
                end if
            else
                if (cell2_dist < cell3_dist) then
                    iCell = cell2
                else
                    iCell = cell3
                end if
            end if
            !write(0,*) 'in nearest_vertex(), current_distance=',current_distance,' iCell=',iCell
            do i = 1, nEdgesOnCell(iCell)
                iVtx = verticesOnCell(i,iCell) 
                d = sphere_distance(latVertex(iVtx), lonVertex(iVtx), target_lat, target_lon, 1.0_RKIND)
                if (d < nearest_distance) then
                    nearest_vertex = iVtx
                    nearest_distance = d
                end if
            end do
        end do
        !write(0,*) 'in nearest_vertex(), nearest_vertex, nearest_distance=',nearest_vertex,nearest_distance

    end function nearest_vertex
                               

    real (kind=RKIND) function sphere_distance(lat1, lon1, lat2, lon2, radius)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Compute the great-circle distance between (lat1, lon1) and (lat2, lon2) on a
    !   sphere with given radius.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
       implicit none
    
        real (kind=RKIND), intent(in) :: lat1, lon1, lat2, lon2, radius
    
        real (kind=RKIND) :: arg1
    
        arg1 = sqrt( sin(0.5_RKIND*(lat2-lat1))**2 +  &
                     cos(lat1)*cos(lat2)*sin(0.5_RKIND*(lon2-lon1))**2 )
        sphere_distance = 2.*radius*asin(arg1)
    
    end function sphere_distance


    real (kind=RKIND) function triangle_area(p1_lat,p1_lon,p2_lat,p2_lon,p3_lat,p3_lon,radius)

        implicit none

        real (kind=RKIND) p1_lat,p1_lon,p2_lat,p2_lon,p3_lat,p3_lon,radius
        real (kind=RKIND) a,b,c,s,e,tanqe

        a = sphere_distance(p1_lat,p1_lon,p2_lat,p2_lon,radius)
        b = sphere_distance(p2_lat,p2_lon,p3_lat,p3_lon,radius)
        c = sphere_distance(p3_lat,p3_lon,p1_lat,p1_lon,radius)
        s = 0.5_RKIND*(a+b+c)

        ! On rare occasions to due numerical precision limitations or compiler
        ! stuff, the product we take the square root of is slightly negative. It happens
        ! along lat=0 and lon=+/-90. The product has always been really small negative,
        ! like smaller than epsilon(x), so setting it to zero isn't changing it much (I guess). 
        ! It doesn't seem to do any harm to set the offending parts of the product to zero.  
        ! One might think that taking the square root of zero
        ! would still result in a NaN, but it doesn't seem to. That's why I use
        ! the max() function to fix it. Maybe instead of 0.0_RKIND, I should use
        ! TINY(X).   But I don't know what X should be.  I got this intrinsic
        ! function from www.nsc.liu.se/~boein/f77to90/a5.html It may depend on
        ! the compiler.  This is when I used PGI.  20131031
        ! Ahijevych
        if (tan(0.5_RKIND*s)*tan(0.5_RKIND*(s-a))*tan(0.5_RKIND*(s-b))*tan(0.5_RKIND*(s-c)) .le. 0) then 
            write(0,*) 'trying to square root a negative in function triangle_area'
            write(0,*) 'p1_lat=',p1_lat,' p1_lon=',p1_lon
            write(0,*) 'p2_lat=',p2_lat,' p2_lon=',p2_lon
            write(0,*) 'p3_lat=',p3_lat,' p3_lon=',p3_lon
            write(0,*) 'radius = ',radius
            write(0,*) ' a=', a, ' b=',b
            write(0,*) ' c=', c, ' s=', s
            write(0,*) ' s-a=', s-a, ' s-b=', s-b
            write(0,*) ' s-c=', s-c
            write(0,*) ' epsilon(s-a)=', epsilon(s-a)
            write(0,*) 'tan(0.5_RKIND*(s))   = ', tan(0.5_RKIND*(s))
            write(0,*) 'tan(0.5_RKIND*(s-a)) = ', tan(0.5_RKIND*(s-a))
            write(0,*) 'tan(0.5_RKIND*(s-b)) = ', tan(0.5_RKIND*(s-b))
            write(0,*) 'tan(0.5_RKIND*(s-c)) = ', tan(0.5_RKIND*(s-c))
            write(0,*) 'attempting to square root', tan(0.5_RKIND*s)*tan(0.5_RKIND*(s-a))*tan(0.5_RKIND*(s-b))*tan(0.5_RKIND*(s-c))
            !tanqe = sqrt(tan(0.5_RKIND*s)*tan(0.5_RKIND*(max(0.0_RKIND,s-a)))*tan(0.5_RKIND*(max(0.0_RKIND,s-b)))*tan(0.5_RKIND*(max(0.0_RKIND,s-c))))
            triangle_area = TINY(triangle_area)
        else
            tanqe = sqrt(tan(0.5_RKIND*s)*tan(0.5_RKIND*(s-a))*tan(0.5_RKIND*(s-b))*tan(0.5_RKIND*(s-c)))
            e = 4.*atan(tanqe)
            triangle_area = radius*radius*e
        endif


    end function triangle_area


    real (kind=RKIND) function triangle_area_1(p1_lat,p1_lon,p2_lat,p2_lon,p3_lat,p3_lon)

        implicit none

        real (kind=RKIND) :: p1_lat,p1_lon,p2_lat,p2_lon,p3_lat,p3_lon

        real (kind=RKIND) :: x1,y1,z1,x2,y2,z2,x3,y3,z3
        real (kind=RKIND) :: angle1, angle2, angle3

        call convert_lx(x1,y1,z1,p1_lat,p1_lon,1.0_RKIND)
        call convert_lx(x2,y2,z2,p2_lat,p2_lon,1.0_RKIND)
        call convert_lx(x3,y3,z3,p3_lat,p3_lon,1.0_RKIND)

        angle1 = sphere_angle(x1,y1,z1,x2,y2,z2,x3,y3,z3)
        angle2 = sphere_angle(x2,y2,z2,x1,y1,z1,x3,y3,z3)
        angle3 = sphere_angle(x3,y3,z3,x2,y2,z2,x1,y1,z1)

        triangle_area_1 = angle1+angle2+angle3 - pii

    end function triangle_area_1


    character(len=256) function basename(input)
       implicit none
       character (len=256), intent(in) :: input
       integer :: i
       logical :: back
       back = .true.
       i = scan(input, '/', back)
       basename = trim(input(i+1:))
      
    end function basename

    real (kind=RKIND) function days_since_1970(input)

       USE datetime_module, ONLY:datetime,date2num
       implicit none
       TYPE(datetime) :: a, b
       character (len=*), intent(in) :: input
       character (len=19) :: ccyy_mm_dd_hh_mm_ss
       integer :: iyear, imonth, iday, ihr, imin, isec, i, xtime_id

       integer :: ncid_tmp

       ! Open netCDF file
       ierr = nf_open(input, NF_NOWRITE, ncid_tmp)
       if (ierr.ne.NF_NOERR) then
          write(0,*) ' error opening mpas file '//input
          call HANDLE_ERR(ierr)
       end if
       ! Get variable id of 'xtime'
       ierr = nf_inq_varid(ncid_tmp,'xtime',xtime_id)
       if (ierr.eq.NF_NOERR) then
          ! Read 'xtime'
          ierr = nf_get_var_text(ncid_tmp,xtime_id,ccyy_mm_dd_hh_mm_ss)
       else
          ! 'xtime' not available, so try parsing filename for time.
          write(0,*)NF_STRERROR(ierr)
          write(0,*) ' no xtime var. Parsing filename '//input
          if (index(input,'diagnostics').gt.0) then
              i = index(input,'diagnostics')
              ccyy_mm_dd_hh_mm_ss = input(i+12:)
          else if (index(input,'diag').gt.0) then
              i = index(input,'diag')
              ccyy_mm_dd_hh_mm_ss = input(i+5:)
          else if (index(input,'mpas.output.').gt.0) then
              i = index(input,'mpas.output.')
              ccyy_mm_dd_hh_mm_ss = input(i+12:)
          else
              write(0,*)'problem with days_since_1970'
              write(0,*)'input='//input
              call exit(1)
          endif
       endif
       ierr = nf_close(ncid_tmp)

       write(0,*) '"'//ccyy_mm_dd_hh_mm_ss//'"'
       read(ccyy_mm_dd_hh_mm_ss,'(I4,1x,I2,1x,I2,1x, I2,1x,I2,1x,I2)') iyear, imonth, iday, ihr, imin, isec
       write(0,*) 'iyear=',iyear,'imonth=',imonth,'iday=',iday,'ihr=',ihr,'imin=',imin,'isec=',isec

       a = datetime(iyear,imonth,iday,ihr,imin)
       b = datetime(1970,1,1,0,0)
       write(0,*) 'a=',a.strftime("%Y %B %d %H %M %S")
       write(0,*) 'b=',b.strftime("%Y %B %d %H %M %S")
       write(0,*) 'date2num(a)=',date2num(a)
       write(0,*) 'date2num(b)=',date2num(b)
       days_since_1970 = date2num(a) - date2num(b)
       write(0,*) 'days_since_1970=',days_since_1970
      
    end function days_since_1970


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! FUNCTION SPHERE_ANGLE
    !
    ! Computes the angle between arcs AB and AC, given points A, B, and C
    ! Equation numbers w.r.t. http://mathworld.wolfram.com/SphericalTrigonometry.html
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real (kind=RKIND) function sphere_angle(ax, ay, az, bx, by, bz, cx, cy, cz)
   
        implicit none
   
        real (kind=RKIND), intent(in) :: ax, ay, az, bx, by, bz, cx, cy, cz
   
        real (kind=RKIND) :: a, b, c          ! Side lengths of spherical triangle ABC
   
        real (kind=RKIND) :: ABx, ABy, ABz    ! The components of the vector AB
        real (kind=RKIND) :: ACx, ACy, ACz    ! The components of the vector AC
   
        real (kind=RKIND) :: Dx               ! The i-components of the cross product AB x AC
        real (kind=RKIND) :: Dy               ! The j-components of the cross product AB x AC
        real (kind=RKIND) :: Dz               ! The k-components of the cross product AB x AC
   
        real (kind=RKIND) :: s                ! Semiperimeter of the triangle
        real (kind=RKIND) :: sin_angle
   
        a = acos(max(min(bx*cx + by*cy + bz*cz,1.0_RKIND),-1.0_RKIND))   ! Eqn. (3)
        b = acos(max(min(ax*cx + ay*cy + az*cz,1.0_RKIND),-1.0_RKIND))   ! Eqn. (2)
        c = acos(max(min(ax*bx + ay*by + az*bz,1.0_RKIND),-1.0_RKIND))   ! Eqn. (1)
   
        ABx = bx - ax
        ABy = by - ay
        ABz = bz - az
   
        ACx = cx - ax
        ACy = cy - ay
        ACz = cz - az
   
        Dx =   (ABy * ACz) - (ABz * ACy)
        Dy = -((ABx * ACz) - (ABz * ACx))
        Dz =   (ABx * ACy) - (ABy * ACx)
   
        s = 0.5*(a + b + c)
        sin_angle = sqrt(min(1.0_RKIND,max(0.0_RKIND,(sin(s-b)*sin(s-c))/(sin(b)*sin(c)))))   ! Eqn. (28)
   
        if ((Dx*ax + Dy*ay + Dz*az) >= 0.0) then
            sphere_angle =  2.0 * asin(max(min(sin_angle,1.0_RKIND),-1.0_RKIND))
        else
            sphere_angle = -2.0 * asin(max(min(sin_angle,1.0_RKIND),-1.0_RKIND))
        end if
   
    end function sphere_angle

     SUBROUTINE HANDLE_ERR(STATUS)
     INTEGER STATUS
     IF (STATUS .NE. NF_NOERR) THEN
       PRINT *, NF_STRERROR(STATUS)
       STOP 'Stopped'
     ENDIF
     END SUBROUTINE HANDLE_ERR
 
end program mpas_to_latlon
