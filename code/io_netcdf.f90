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
    public input_netcdf, netcdf_o, netcdf_algae_o
    
    type netcdf_o
        private
        !netCDF file id
        integer                                 :: nc_id
        !parameter_ids
        integer, allocatable                    :: parameter_id(:)
        integer, allocatable                    :: parameter_id_diag(:)
        integer                                 :: z_id, time_id, iz_id
        integer                                 :: t_id, s_id, kz2_id
        !auxiliary
        logical                                 :: first
    contains
        private
        procedure, public:: init_netcdf
        procedure, public:: save_netcdf
        procedure, public:: close_netcdf
    end type netcdf_o
   
    type netcdf_algae_o
        private
        !netCDF file id
        integer                                 :: nc_id
        !parameter_ids
        integer                                 :: time_id, ice_id
        integer                                 :: z_viz_id, z_id
        integer                                 :: par_z_id
        integer                                 :: bulk_temperature_id, bulk_salinity_id, bulk_density_id
        integer                                 :: brine_temperature_id, brine_salinity_id, brine_density_id
        integer                                 :: a_carbon_id, a_nitrogen_id, a_phosphorus_id
        integer                                 :: nh4_id, no2_id, no3_id, po4_id
        integer                                 :: brine_relative_volume_id
        !auxiliary
        logical                                 :: first
    contains
        private
        procedure, public:: init_netcdf_algae
        procedure, public:: save_netcdf_algae
        procedure, public:: close_netcdf_algae
    end type netcdf_algae_o
    
    interface netcdf_o
        procedure constructor_netcdf_o
    end interface
    
    interface netcdf_algae_o
        procedure constructor_netcdf_algae_o
    end interface
    
    contains
    
    subroutine input_netcdf(file_name, z, dz, kz_bio, lev_max, t, s, AKs, hice, boundary_bbl_sediments, &
        boundary_water_bbl, width_bbl, resolution_bbl, width_bioturbation, resolution_bioturbation, &
        width_sediments, resolution_sediments, year, ice_area, heat_flux, snow_thick, t_ice)
    
    implicit none
    real(rk), intent(in) :: width_bbl, resolution_bbl, width_bioturbation, resolution_bioturbation, &
        width_sediments, resolution_sediments
    
    ! This is the name of the data file we will read.
    character (len = *), intent(in)                         :: file_name
    real(rk), dimension (:, :), pointer, intent(out)        :: t, s, AKs
    real(rk), dimension (:), pointer, intent(out)           :: z, dz, kz_bio, hice
    real(rk), dimension (:), pointer, intent(out)           :: ice_area, heat_flux, snow_thick, t_ice
    integer, intent(out)                                    :: lev_max, boundary_bbl_sediments, boundary_water_bbl
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
    boundary_water_bbl = h_rec + 2
    
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
    
    function constructor_netcdf_o()
    
    implicit none
    type(netcdf_o), pointer:: constructor_netcdf_o
    
    allocate(constructor_netcdf_o)
    
    end function constructor_netcdf_o
    
    function constructor_netcdf_algae_o()
    
    implicit none
    type(netcdf_algae_o), pointer:: constructor_netcdf_algae_o
    
    allocate(constructor_netcdf_algae_o)
    
    end function constructor_netcdf_algae_o
    
    subroutine init_netcdf(self, fn, first_lvl, last_lvl, model)

    implicit none
    class(netcdf_o):: self
    !input:
    character(len = *), intent(in)          :: fn
    integer, intent(in)                     :: first_lvl, last_lvl
    type (type_model),intent(in)            :: model
    !dimension ids
    integer                                 :: z_dim_id, time_dim_id
    integer                                 :: ip, ilast, nlev
    !dimension lengths
    integer, parameter                      :: time_len = NF90_UNLIMITED
    integer                                 :: dim1d
    integer                                 :: dim_ids(2)

    nlev = last_lvl - first_lvl + 1
    self%first = .true.
    self%nc_id = -1
    call check_err(nf90_create(fn, NF90_CLOBBER, self%nc_id))
    !define the dimensions
    call check_err(nf90_def_dim(self%nc_id, "z", nlev, z_dim_id))
    call check_err(nf90_def_dim(self%nc_id, "time", time_len, time_dim_id))
    !define coordinates
    dim1d = z_dim_id
    call check_err(nf90_def_var(self%nc_id, "z", NF90_REAL, dim1d, self%z_id))
    dim1d = time_dim_id
    call check_err(nf90_def_var(self%nc_id, "time", NF90_REAL, dim1d, self%time_id))
    !define variables
    dim_ids(1) = z_dim_id
    dim_ids(2) = time_dim_id
    allocate(self%parameter_id(size(model%state_variables)))
    do ip = 1, size(model%state_variables)
        ilast = index(model%state_variables(ip)%path,'/',.true.)
        call check_err(nf90_def_var(self%nc_id, model%state_variables(ip)%path(ilast+1:), NF90_REAL, dim_ids, self%parameter_id(ip)))
        call check_err(set_attributes(ncid=self%nc_id, id=self%parameter_id(ip), units=model%state_variables(ip)%units, &
            long_name=model%state_variables(ip)%long_name,missing_value=model%state_variables(ip)%missing_value))
    end do
    allocate(self%parameter_id_diag(size(model%diagnostic_variables)))
    do ip = 1, size(model%diagnostic_variables)
        if (model%diagnostic_variables(ip)%save) then
            ilast = index(model%diagnostic_variables(ip)%path,'/',.true.)
            call check_err(nf90_def_var(self%nc_id, model%diagnostic_variables(ip)%path(ilast+1:), NF90_REAL, dim_ids, self%parameter_id_diag(ip)))
            call check_err(set_attributes(ncid=self%nc_id, id=self%parameter_id_diag(ip), units=model%diagnostic_variables(ip)%units, &
                long_name=model%diagnostic_variables(ip)%long_name,missing_value=model%diagnostic_variables(ip)%missing_value))
        end if
    end do
    call check_err(nf90_def_var(self%nc_id, "t", NF90_REAL, dim_ids, self%t_id))
    call check_err(nf90_def_var(self%nc_id, "s", NF90_REAL, dim_ids, self%s_id))
    call check_err(nf90_def_var(self%nc_id, "Kz2", NF90_REAL, dim_ids, self%kz2_id))
    call check_err(nf90_def_var(self%nc_id, "radiative_flux", NF90_REAL, dim_ids, self%iz_id))
    !end define
    call check_err(nf90_enddef(self%nc_id))
    
    end subroutine init_netcdf
    
    subroutine init_netcdf_algae(self, fn, first_lvl, last_lvl)

    implicit none
    class(netcdf_algae_o):: self
    !input:
    character(len = *), intent(in)          :: fn
    integer, intent(in)                     :: first_lvl, last_lvl
    !dimension ids
    integer                                 :: z_dim_id, time_dim_id
    integer                                 :: ip, nlev
    !dimension lengths
    integer, parameter                      :: time_len = NF90_UNLIMITED
    integer                                 :: dim1d
    integer                                 :: dim_ids(2)

    nlev = last_lvl - first_lvl + 1
    self%first = .true.
    self%nc_id = -1
    call check_err(nf90_create(fn, NF90_CLOBBER, self%nc_id))
    !define the dimensions
    call check_err(nf90_def_dim(self%nc_id, "z", nlev, z_dim_id))
    call check_err(nf90_def_dim(self%nc_id, "time", time_len, time_dim_id))
    !define coordinates
    dim1d = z_dim_id
    call check_err(nf90_def_var(self%nc_id, "z", NF90_REAL, dim1d, self%z_viz_id))
    dim1d = time_dim_id
    call check_err(nf90_def_var(self%nc_id, "time", NF90_REAL, dim1d, self%time_id))
    call check_err(nf90_def_var(self%nc_id, "ice", NF90_REAL, dim1d, self%ice_id))
    !define variables
    dim_ids(1) = z_dim_id
    dim_ids(2) = time_dim_id
    call check_err(nf90_def_var(self%nc_id, "layers' depth", NF90_REAL, dim_ids, self%z_id))
    call check_err(nf90_def_var(self%nc_id, "bulk_temperature", NF90_REAL, dim_ids, self%bulk_temperature_id))
    call check_err(nf90_def_var(self%nc_id, "bulk_salinity", NF90_REAL, dim_ids, self%bulk_salinity_id))
    call check_err(nf90_def_var(self%nc_id, "bulk_density", NF90_REAL, dim_ids, self%bulk_density_id))
    call check_err(nf90_def_var(self%nc_id, "photosynthetic_active_radiation", NF90_REAL, dim_ids, self%par_z_id))
    call check_err(nf90_def_var(self%nc_id, "brine_temperature", NF90_REAL, dim_ids, self%brine_temperature_id))
    call check_err(nf90_def_var(self%nc_id, "brine_salinity", NF90_REAL, dim_ids, self%brine_salinity_id))
    call check_err(nf90_def_var(self%nc_id, "brine_density", NF90_REAL, dim_ids, self%brine_density_id))
    call check_err(nf90_def_var(self%nc_id, "algae_carbon", NF90_REAL, dim_ids, self%a_carbon_id))
    call check_err(nf90_def_var(self%nc_id, "algae_nitrogen", NF90_REAL, dim_ids, self%a_nitrogen_id))
    call check_err(nf90_def_var(self%nc_id, "algae_phosphorus", NF90_REAL, dim_ids, self%a_phosphorus_id))
    call check_err(nf90_def_var(self%nc_id, "nh4", NF90_REAL, dim_ids, self%nh4_id))
    call check_err(nf90_def_var(self%nc_id, "no2", NF90_REAL, dim_ids, self%no2_id))
    call check_err(nf90_def_var(self%nc_id, "no3", NF90_REAL, dim_ids, self%no3_id))
    call check_err(nf90_def_var(self%nc_id, "po4", NF90_REAL, dim_ids, self%po4_id))
    call check_err(nf90_def_var(self%nc_id, "brine_relative_volume", NF90_REAL, dim_ids, self%brine_relative_volume_id))
    !end define
    call check_err(nf90_enddef(self%nc_id))
    
    end subroutine init_netcdf_algae
    
    subroutine save_netcdf(self, first_lvl, last_lvl, lev_max, julianday, Cc, tem2, sal2, Kz2, model, z, iz)

    implicit none
    class(netcdf_o):: self
    integer, intent(in)                     :: first_lvl, last_lvl, julianday, lev_max
    real(rk), dimension(:, :), intent(in)   :: Cc
    real(rk), dimension(:, :), intent(in)   :: tem2, sal2, Kz2
    type (type_model), intent(in)           :: model
    real(rk), dimension(:), intent(in)      :: z, iz
    
    integer                                 :: ip, i
    integer                                 :: edges(2), start(2), start_time(1), edges_time(1)
    real(rk)                                :: temp_matrix(lev_max)
    real                                    :: foo(1) !nevermind what is it but it works
        
    !write data
    edges(1) = last_lvl - first_lvl + 1
    edges(2) = 1
    start(1) = 1
    start(2) = julianday
    start_time(1) = julianday
    edges_time(1) = 1

    if ( self%first ) then
        call check_err(nf90_put_var(self%nc_id, self%z_id,  z(first_lvl:last_lvl), start, edges))
        self%first = .false.
    end if
    
    foo(1) = real(julianday)
    if (self%nc_id .ne. -1) then
        call check_err(nf90_put_var(self%nc_id, self%time_id, foo, start_time, edges_time))
        
        do ip = 1, size(model%state_variables)
            call check_err(nf90_put_var(self%nc_id, self%parameter_id(ip), Cc(first_lvl:last_lvl, ip), start, edges))
        end do
        do ip = 1, size(model%diagnostic_variables)
            if (model%diagnostic_variables(ip)%save) then
                temp_matrix = fabm_get_bulk_diagnostic_data(model,ip)
                call check_err(nf90_put_var(self%nc_id, self%parameter_id_diag(ip), temp_matrix(first_lvl:last_lvl), start, edges))
            end if
        end do
        call check_err(nf90_put_var(self%nc_id, self%t_id,   tem2(first_lvl:last_lvl, julianday), start, edges))
        call check_err(nf90_put_var(self%nc_id, self%s_id,   sal2(first_lvl:last_lvl, julianday), start, edges))
        call check_err(nf90_put_var(self%nc_id, self%kz2_id, Kz2(first_lvl:last_lvl, julianday),  start, edges))
        call check_err(nf90_put_var(self%nc_id, self%iz_id,  iz(first_lvl:last_lvl),  start, edges))
        call check_err(nf90_sync(self%nc_id))        
    end if
    
    end subroutine save_netcdf
    
    subroutine save_netcdf_algae(self, ice_l, first_lvl, last_lvl, lev_max, julianday, ice)

    use ice_algae_lib
    
    implicit none
    class(netcdf_algae_o):: self
    type(ice_layer), pointer, dimension(:), intent(in)  :: ice_l
    integer, intent(in)                                 :: first_lvl, last_lvl, lev_max, julianday
    real(rk), intent(in)                                :: ice
    !internal parameters
    real(rk), dimension(lev_max)                  :: z, z_viz
    real(rk), dimension(lev_max)                  :: par_z
    real(rk), dimension(lev_max)                  :: bulk_temperature, bulk_salinity, bulk_density
    real(rk), dimension(lev_max)                  :: brine_temperature, brine_salinity, brine_density
    real(rk), dimension(lev_max)                  :: a_carbon, a_nitrogen, a_phosphorus
    real(rk), dimension(lev_max)                  :: nh4, no2, no3, po4
    real(rk), dimension(lev_max)                  :: brine_relative_volume

    integer                                 :: ip, i
    integer                                 :: edges(2), start(2), start_time(1), edges_time(1)
    real(rk)                                :: temp_matrix(lev_max)
    real                                    :: foo(1), bar(1) !nevermind what is it but it works
    
    !call ice_l%get_algae(z, par_z, bulk_temperature, bulk_salinity, &
    !    bulk_density, brine_temperature, brine_salinity, brine_density, &
    !    a_carbon, a_nitrogen, a_phosphorus, nh4, no2, no3, po4, brine_relative_volume)
    
    !write data
    edges(1) = last_lvl - first_lvl + 1
    edges(2) = 1
    start(1) = 1
    start(2) = julianday
    start_time(1) = julianday
    edges_time(1) = 1

    z_viz(1) = 1.
    do i = 2, lev_max
        z_viz(i) = z_viz(i - 1) + 1.
    end do
    
    if ( self%first .and. self%nc_id .ne. -1) then
        call check_err(nf90_put_var(self%nc_id, self%z_viz_id,  z_viz, start, edges))
        self%first = .false.
    end if
    
    foo(1) = real(julianday)
    bar(1) = real(ice)
    if (self%nc_id .ne. -1) then
        call check_err(nf90_put_var(self%nc_id, self%z_id,  z(first_lvl:last_lvl), start, edges))
        call check_err(nf90_put_var(self%nc_id, self%time_id, foo, start_time, edges_time))
        call check_err(nf90_put_var(self%nc_id, self%ice_id,  bar, start_time, edges_time))
        call check_err(nf90_put_var(self%nc_id, self%par_z_id, par_z(first_lvl:last_lvl), start, edges))
        call check_err(nf90_put_var(self%nc_id, self%bulk_temperature_id, bulk_temperature(first_lvl:last_lvl), start, edges))
        call check_err(nf90_put_var(self%nc_id, self%bulk_salinity_id, bulk_salinity(first_lvl:last_lvl), start, edges))
        call check_err(nf90_put_var(self%nc_id, self%bulk_density_id, bulk_density(first_lvl:last_lvl), start, edges))
        call check_err(nf90_put_var(self%nc_id, self%brine_temperature_id, brine_temperature(first_lvl:last_lvl), start, edges))
        call check_err(nf90_put_var(self%nc_id, self%brine_salinity_id, brine_salinity(first_lvl:last_lvl), start, edges))
        call check_err(nf90_put_var(self%nc_id, self%brine_density_id, brine_density(first_lvl:last_lvl), start, edges))
        call check_err(nf90_put_var(self%nc_id, self%a_carbon_id, a_carbon(first_lvl:last_lvl), start, edges))
        call check_err(nf90_put_var(self%nc_id, self%a_nitrogen_id, a_nitrogen(first_lvl:last_lvl), start, edges))
        call check_err(nf90_put_var(self%nc_id, self%a_phosphorus_id, a_phosphorus(first_lvl:last_lvl), start, edges))
        call check_err(nf90_put_var(self%nc_id, self%nh4_id, nh4(first_lvl:last_lvl), start, edges))
        call check_err(nf90_put_var(self%nc_id, self%no2_id, no2(first_lvl:last_lvl), start, edges))
        call check_err(nf90_put_var(self%nc_id, self%no3_id, no3(first_lvl:last_lvl), start, edges))
        call check_err(nf90_put_var(self%nc_id, self%po4_id, po4(first_lvl:last_lvl), start, edges))
        call check_err(nf90_put_var(self%nc_id, self%brine_relative_volume_id, brine_relative_volume(first_lvl:last_lvl), start, edges))
        call check_err(nf90_sync(self%nc_id))        
    end if
    
    end subroutine save_netcdf_algae
    
    subroutine close_netcdf(self)
        
    implicit none
    class(netcdf_o):: self
    if (self%nc_id .ne. -1) then
        call check_err(nf90_close(self%nc_id))
        deallocate(self%parameter_id)
        deallocate(self%parameter_id_diag)
    end if
    self%nc_id = -1
        
    end subroutine close_netcdf
    
    subroutine close_netcdf_algae(self)
        
    implicit none
    class(netcdf_algae_o):: self
    if (self%nc_id .ne. -1) then
        call check_err(nf90_close(self%nc_id))
    end if
    self%nc_id = -1
        
    end subroutine close_netcdf_algae
    
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
