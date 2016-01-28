    module io_netcdf
    
    !REVISION HISTORY:
    !Original author(s): Shamil Yakubov
    !to do:
    
    use netcdf
    use fabm, only:type_model,fabm_get_bulk_diagnostic_data
    use fabm_types, only: rk
    
    implicit none
    !all is private
    private
    !public functions
    public input_netcdf, init_netcdf, save_netcdf, close_netcdf
    !netCDF file id
    integer                                 :: nc_id
    integer, allocatable                    :: parameter_id(:)
    integer, allocatable                    :: parameter_id_diag(:)
    !parameter_ids
    integer                                 :: z_id, time_id, ice_id, iz_id
    integer                                 :: pH_id, T_id, S_id, Kz2_id
    integer                                 :: pCO2_id, Om_Ca_id, Om_Ar_id
    integer                                 :: ice_algae_id, ice_dom_id
    
    logical                                 :: first
   
    contains
    
    subroutine input_netcdf(file_name, z, dz, kz_bio, lev_max, t, s, AKs, hice, boundary_bbl_sediments, &
        width_bbl, resolution_bbl, width_bioturbation, resolution_bioturbation, &
        width_sediments, resolution_sediments, year, ice_area, heat_flux, snow_thick, t_ice)
    
    implicit none
         
    real(rk), intent(in) :: width_bbl, resolution_bbl, width_bioturbation, resolution_bioturbation, &
        width_sediments, resolution_sediments
    
    ! This is the name of the data file we will read.
    character (len = *), intent(in)                         :: file_name
    real(rk), dimension (:, :), pointer, intent(out)        :: t, s, AKs
    real(rk), dimension (:), pointer, intent(out)           :: z, dz, kz_bio, hice
    real(rk), dimension (:), pointer, intent(out)           :: ice_area, heat_flux, snow_thick, t_ice
    integer, intent(out)                                    :: lev_max, boundary_bbl_sediments
    integer, intent(in)                                     :: year
    
    integer                                     :: ncid, i, j, range, number_of_days
    integer                                     :: lat_rec, lon_rec, time_rec, h_rec
    integer                                     :: t_varid, s_varid, hice_varid, AKs_varid, depth_varid
    integer                                     :: ice_area_id, heat_flux_id, snow_thick_id, t_ice_id
    integer, dimension(nf90_max_var_dims)       :: dimids
    real(rk), allocatable, dimension(:, :)      :: t_temp, s_temp, AKs_temp
    real(rk), allocatable, dimension(:)         :: z_temp, hice_temp
    real(rk), allocatable, dimension(:)         :: ice_area_temp, heat_flux_temp, snow_thick_temp, t_ice_temp
    
    integer                                     :: bbl_count, bioturbation_count, sediments_count
    real(rk)                                    :: l !auxiliary parameter
    
    call check_err(nf90_open(FILE_NAME, NF90_NOWRITE, ncid))
    call check_err(nf90_inq_varid(ncid, "temp", t_varid))
    call check_err(nf90_inq_varid(ncid, "salt", s_varid))
    call check_err(nf90_inq_varid(ncid, "hice", hice_varid))
    call check_err(nf90_inq_varid(ncid, "aice", ice_area_id))
    call check_err(nf90_inq_varid(ncid, "shflux", heat_flux_id))
    call check_err(nf90_inq_varid(ncid, "snow_thick", snow_thick_id))
    call check_err(nf90_inq_varid(ncid, "tisrf", t_ice_id))
    call check_err(nf90_inq_varid(ncid, "AKs", AKs_varid))
    call check_err(nf90_inq_varid(ncid, "depth", depth_varid))
    call check_err(nf90_inquire_variable(ncid, t_varid, dimids = dimids))
    call check_err(nf90_inquire_dimension(ncid, dimids(1), len = h_rec))
    call check_err(nf90_inquire_dimension(ncid, dimids(2), len = time_rec))
    
    bbl_count = int(width_bbl/resolution_bbl) + h_rec + 1 ! +1 - for surface with depth 0 layer
    bioturbation_count = int(width_bioturbation/resolution_bioturbation) + bbl_count
    sediments_count = int(width_sediments/resolution_sediments) + bioturbation_count
    !temporary variables
    allocate(t_temp(h_rec, time_rec))
    allocate(s_temp(h_rec, time_rec))
    allocate(hice_temp(time_rec))
    allocate(ice_area_temp(time_rec))
    allocate(heat_flux_temp(time_rec))
    allocate(snow_thick_temp(time_rec))
    allocate(t_ice_temp(time_rec))
    allocate(AKs_temp(h_rec, time_rec))
    allocate(z_temp(h_rec))
    !state variables
    allocate(t(sediments_count, time_rec))
    allocate(s(sediments_count, time_rec))
    allocate(hice(365))
    allocate(ice_area(365))
    allocate(heat_flux(365))
    allocate(snow_thick(365))
    allocate(t_ice(365))
    allocate(AKs(sediments_count, time_rec))
    allocate(z(sediments_count))
    allocate(dz(sediments_count))
    allocate(kz_bio(sediments_count))
    !getting from file
    call check_err(nf90_get_var(ncid, t_varid, t_temp))
    call check_err(nf90_get_var(ncid, s_varid, s_temp))
    call check_err(nf90_get_var(ncid, hice_varid, hice_temp))
    call check_err(nf90_get_var(ncid, ice_area_id, ice_area_temp))
    call check_err(nf90_get_var(ncid, heat_flux_id, heat_flux_temp))
    call check_err(nf90_get_var(ncid, snow_thick_id, snow_thick_temp))
    call check_err(nf90_get_var(ncid, t_ice_id, t_ice_temp))
    call check_err(nf90_get_var(ncid, AKs_varid, AKs_temp))
    call check_err(nf90_get_var(ncid, depth_varid, z_temp))
    !creating grid
    
    range = 352 !adding 352 days to initial data (12:00 15 january 1980 become 12:00 1 january 1981)
    if (year < 1981 .and. year > 2012) then
        print *, "Wrong year, it should be between 1981 and 2012"
        stop
    end if
    if (year == 1981) then
        continue
    else
        do i = 1981, year - 1
            number_of_days = 365
            if (mod(i, 4) == 0) number_of_days = 366
            range = range + number_of_days
        end do
    end if
    !i - days from data file
    do i = 1, 365
        !temperature, salinity, turbulence coefficients in water column
        do j = h_rec, 1, -1  
            t(h_rec - j + 2, i) = t_temp(j, i + range)
            s(h_rec - j + 2, i) = s_temp(j, i + range)
            AKs(h_rec - j + 2, i) = AKs_temp(j, i + range)
        end do
        t(1, i) = t(2, i)
        s(1, i) = s(2, i)
        AKs(1, i) = AKs(2, i)
        AKs(h_rec + 1, i) = 0.02 * AKs(h_rec, i)
        !temperature, salinity, turbulence coefficients in bottom boundary layer(bbl)
        do j = h_rec + 2, bbl_count
            t(j, i) = t(h_rec + 1, i)
            s(j, i) = s(h_rec + 1, i)
            AKs(j, i) = 0.6E-6
        end do
        AKs(bbl_count, i) = 0.6E-8
        !temperature, salinity, turbulence coefficients in sediments with bioturbation
        do j = bbl_count + 1, bioturbation_count
            t(j, i) = t(h_rec, i)
            s(j, i) = s(h_rec, i)
            AKs(j, i) = 1.E-11
        end do
        !temperature, salinity, turbulence coefficients in sediments without bioturbation
        do j = bioturbation_count + 1, sediments_count
            t(j, i) = t(h_rec, i)
            s(j, i) = s(h_rec, i)
            AKs(j, i) = 1.E-11
        end do
        hice(i) = hice_temp(i + range)
        ice_area(i) = ice_area_temp(i + range)
        heat_flux(i) = heat_flux_temp(i + range)
        snow_thick(i) = snow_thick_temp(i + range)
        t_ice(i) = t_ice_temp(i + range)
    end do
    !writing depth and kz_bio to array in water column
    do j = h_rec, 1, -1
        z(h_rec - j + 2) = z_temp(j)
        kz_bio(j + 1) = 0.
    end do
    z(1) = 0.
    kz_bio(1) = 0.
    !writing depth and kz_bio to array in bbl 
    do j = h_rec + 2, bbl_count
        z(j) = z(j - 1) + resolution_bbl
        kz_bio(j) = 0.
    end do
    !for calculation module
    boundary_bbl_sediments = j
    !writing depth and kz_bio to array in sediments with bioturbation
    do j = bbl_count + 1, bioturbation_count
        z(j) = z(j - 1) + resolution_bioturbation
        kz_bio(j) = 1.E-11
    end do
    !writing depth and kz_bio to array in sediments with no bioturbation
    do j = bioturbation_count + 1, sediments_count
        l = (j - bioturbation_count + 1) / 2
        z(j) = z(j - 1) + resolution_sediments
        kz_bio(j) = 1.E-11*exp(-l) 
    end do
    !deepest lvl & dz array calculating
    lev_max = size(z)
    do j = 2, lev_max
        dz(j - 1) = z(j) - z(j - 1)
    end do
    dz(lev_max) = 0
    
    deallocate(t_temp)
    deallocate(s_temp)
    deallocate(hice_temp)
    deallocate(AKs_temp)
    deallocate(z_temp)
    deallocate(ice_area_temp)
    deallocate(heat_flux_temp)
    deallocate(snow_thick_temp)
    deallocate(t_ice_temp)
    
    call check_err(nf90_close(ncid))
    
    end subroutine input_netcdf
        
    subroutine init_netcdf(fn, nlev, model)

    implicit none
    
    !input:
    character(len = *), intent(in)          :: fn
    integer, intent(in)                     :: nlev
    type (type_model),intent(in)            :: model
    !dimension ids
    integer                                 :: z_dim_id, time_dim_id
    integer                                 :: ip, iret, ilast
    !dimension lengths
    integer, parameter                      :: time_len = NF90_UNLIMITED
    character(10)                           :: parameter_name
    integer                                 :: dim1d
    integer                                 :: dim_ids(2)
    
    first = .true.
    print *, 'NetCDF version: ', trim(nf90_inq_libvers())
    nc_id = -1
    call check_err(nf90_create(fn, NF90_CLOBBER, nc_id))
    !define the dimensions
    call check_err(nf90_def_dim(nc_id, "z", nlev, z_dim_id))
    call check_err(nf90_def_dim(nc_id, "time", time_len, time_dim_id))
    !define coordinates
    dim1d = z_dim_id
    call check_err(nf90_def_var(nc_id, "z", NF90_REAL, dim1d, z_id))
    dim1d = time_dim_id
    call check_err(nf90_def_var(nc_id, "time", NF90_REAL, dim1d, time_id))
    call check_err(nf90_def_var(nc_id, "ice", NF90_REAL, dim1d, ice_id))
    call check_err(nf90_def_var(nc_id, "ice algae", NF90_REAL, dim1d, ice_algae_id))
    call check_err(nf90_def_var(nc_id, "ice DOM", NF90_REAL, dim1d, ice_dom_id))
    !define variables
    dim_ids(1) = z_dim_id
    dim_ids(2) = time_dim_id
    allocate(parameter_id(size(model%state_variables)))
    do ip = 1, size(model%state_variables)
        ilast = index(model%state_variables(ip)%path,'/',.true.)
        call check_err(nf90_def_var(nc_id, model%state_variables(ip)%path(ilast+1:), NF90_REAL, dim_ids, parameter_id(ip)))
        call check_err(set_attributes(ncid=nc_id, id=parameter_id(ip), units=model%state_variables(ip)%units, &
            long_name=model%state_variables(ip)%long_name,missing_value=model%state_variables(ip)%missing_value))
    end do
    allocate(parameter_id_diag(size(model%diagnostic_variables)))
    do ip = 1, size(model%diagnostic_variables)
        if (model%diagnostic_variables(ip)%save) then
            ilast = index(model%diagnostic_variables(ip)%path,'/',.true.)
            call check_err(nf90_def_var(nc_id, model%diagnostic_variables(ip)%path(ilast+1:), NF90_REAL, dim_ids, parameter_id_diag(ip)))
            call check_err(set_attributes(ncid=nc_id, id=parameter_id_diag(ip), units=model%diagnostic_variables(ip)%units, &
                long_name=model%diagnostic_variables(ip)%long_name,missing_value=model%diagnostic_variables(ip)%missing_value))
        end if
    end do
    call check_err(nf90_def_var(nc_id, "T", NF90_REAL, dim_ids, T_id))
    call check_err(nf90_def_var(nc_id, "S", NF90_REAL, dim_ids, S_id))
    call check_err(nf90_def_var(nc_id, "Kz2", NF90_REAL, dim_ids, Kz2_id))
    call check_err(nf90_def_var(nc_id, "radiative_flux", NF90_REAL, dim_ids, iz_id))
    !end define
    call check_err(nf90_enddef(nc_id))
    
    end subroutine init_netcdf
    
    subroutine save_netcdf(LevMax, julianday, Cc, tem2, sal2, Kz2, model, z, ice, iz, ice_algae, ice_dom)

    implicit none
    
    integer, intent(in)                     :: LevMax, julianday
    real(rk), dimension(:, :), intent(in)   :: Cc
    real(rk), dimension(:, :), intent(in)   :: tem2, sal2, Kz2
    type (type_model), intent(in)           :: model
    real(rk), dimension(:), intent(in)      :: z, ice, iz
    real(rk), intent(in)                    :: ice_algae, ice_dom
    
    integer                                 :: edges(2), start(2), start_time(1), edges_time(1)
    real                                    :: temp_matrix(LevMax), dum(1), mud(1), foo(1), bar(1) ! nevermind what is it but it works
    integer                                 :: ip, i
        
    !write data
    edges(1) = LevMax
    edges(2) = 1
    start(1) = 1
    start(2) = julianday
    start_time(1) = julianday
    edges_time(1) = 1

    if ( first ) then
        !do ip = 1, LevMax
        !    temp_matrix(ip) = ip
        !end do
        !temp_matrix = (-1 * temp_matrix)
        !call check_err(nf90_put_var(nc_id, z_id, temp_matrix, start, edges))
        call check_err(nf90_put_var(nc_id, z_id,  z(1:LevMax), start, edges))
        first = .false.
    end if
    
    dum(1) = real(julianday)
    mud(1) = real(ice(julianday))
    foo(1) = real(ice_algae)
    bar(1) = real(ice_dom)
    if (nc_id .ne. -1) then
        call check_err(nf90_put_var(nc_id, time_id, dum, start_time, edges_time))
        call check_err(nf90_put_var(nc_id, ice_id,  mud, start_time, edges_time))
        call check_err(nf90_put_var(nc_id, ice_algae_id, foo,  start_time, edges_time))
        call check_err(nf90_put_var(nc_id, ice_dom_id,   bar,  start_time, edges_time))
        
        do ip = 1, size(model%state_variables)
            call check_err(nf90_put_var(nc_id, parameter_id(ip), Cc(1:LevMax, ip), start, edges))
        end do
        do ip = 1, size(model%diagnostic_variables)
            if (model%diagnostic_variables(ip)%save) then
                temp_matrix = fabm_get_bulk_diagnostic_data(model,ip)
                call check_err(nf90_put_var(nc_id, parameter_id_diag(ip), temp_matrix(1:LevMax), start, edges))
            end if
        end do
        call check_err(nf90_put_var(nc_id, T_id,   tem2(1:LevMax, julianday), start, edges))
        call check_err(nf90_put_var(nc_id, S_id,   sal2(1:LevMax, julianday), start, edges))
        call check_err(nf90_put_var(nc_id, Kz2_id, Kz2(1:LevMax, julianday),  start, edges))
        call check_err(nf90_put_var(nc_id, iz_id,  iz(:),  start, edges))
        call check_err(nf90_sync(nc_id))        
    end if
        
    end subroutine save_netcdf
    
    subroutine close_netcdf()
        
    if (nc_id .ne. -1) then
        call check_err(nf90_close(nc_id))
        deallocate(parameter_id)
        write (*,'(a)') "finished"
    end if
    nc_id = -1
        
    end subroutine close_netcdf
    
    integer function set_attributes(ncid,id,                         &
                                    units,long_name,                 &
                                    valid_min,valid_max,valid_range, &
                                    scale_factor,add_offset,         &
                                    FillValue,missing_value,         &
                                    C_format,FORTRAN_format)
    !
    ! !DESCRIPTION:
    !  This routine is used to set a number of attributes for
    !  variables. The routine makes heavy use of the {\tt optional} keyword.
    !  The list of recognized keywords is very easy to extend. We have
    !  included a sub-set of the COARDS conventions.
    !
    ! !USES:
    !  IMPLICIT NONE
    !
    ! !INPUT PARAMETERS:
    integer, intent(in)                     :: ncid,id
    character(len=*), optional              :: units,long_name
    real, optional                          :: valid_min,valid_max
    real, optional                          :: valid_range(2)
    real, optional                          :: scale_factor,add_offset
    double precision, optional              :: FillValue,missing_value
    character(len=*), optional              :: C_format,FORTRAN_format
    !
    ! !REVISION HISTORY:
    !  Original author(s): Karsten Bolding & Hans Burchard
    !
    ! !LOCAL VARIABLES:
    integer                                 :: iret
    real                                    :: vals(2)
    !
    !
    !-----------------------------------------------------------------------
    !
    if(present(units)) then
    iret = nf90_put_att(ncid,id,'units',trim(units))
    end if

    if(present(long_name)) then
    iret = nf90_put_att(ncid,id,'long_name',trim(long_name))
    end if

    if(present(C_format)) then
    iret = nf90_put_att(ncid,id,'C_format',trim(C_format))
    end if

    if(present(FORTRAN_format)) then
    iret = nf90_put_att(ncid,id,'FORTRAN_format',trim(FORTRAN_format))
    end if

    if(present(valid_min)) then
    vals(1) = valid_min
    iret = nf90_put_att(ncid,id,'valid_min',vals(1:1))
    end if

    if(present(valid_max)) then
    vals(1) = valid_max
    iret = nf90_put_att(ncid,id,'valid_max',vals(1:1))
    end if

    if(present(valid_range)) then
    vals(1) = valid_range(1)
    vals(2) = valid_range(2)
    iret = nf90_put_att(ncid,id,'valid_range',vals(1:2))
    end if

    if(present(scale_factor)) then
    vals(1) = scale_factor
    iret = nf90_put_att(ncid,id,'scale_factor',vals(1:1))
    end if

    if(present(add_offset)) then
    vals(1) = add_offset
    iret = nf90_put_att(ncid,id,'add_offset',vals(1:1))
    end if

    if(present(FillValue)) then
    vals(1) = FillValue
    iret = nf90_put_att(ncid,id,'_FillValue',vals(1:1))
    end if

    if(present(missing_value)) then
    vals(1) = missing_value
    iret = nf90_put_att(ncid,id,'missing_value',vals(1:1))
    end if

    set_attributes = 0
    
    return
    end function set_attributes
    
    subroutine check_err(status)

    implicit none
    
    integer, intent (in) :: status
    
    if (status .ne. NF90_NOERR) then
        print *, trim(nf90_strerror(status))
        stop
    endif
        
    end subroutine check_err
        
    end module io_netcdf
    
!-----------------------------------------------------------------------
! Copyright under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------