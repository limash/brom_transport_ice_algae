    module brom_transport
    !REVISION HISTORY:
    !Original author(s): Evgeniy Yakushev, Shamil Yakubov
    !Shamil_11_2015 - yaml configurator was added
    !Shamil_12_2015 - ice dependencies were added
    !Shamil_01_2016 - ice_algae object was added
    !to do:
    ! - 
    
    use fabm
    use fabm_config
    use fabm_types, only: attribute_length, rk
    use io_netcdf
    use io_ascii
    use ice_algae_lib
    
    implicit none
    private
    public init_brom_transport, do_brom_transport, clear_brom_transport
    !FABM model with all data and procedures related to biogeochemistry
    type (type_model)                                               :: model
    integer                                                         :: lev_max, par_max, year, number_of_layers
    integer                                                         :: boundary_bbl_sediments, boundary_water_bbl
    !grids for depth(z), gradients of depth(dz), bioturbations coefficients, hice -  "time-averaged average ice thickness in cell"
    real(rk), pointer, dimension(:)                                 :: z, dz, kz_bio, hice
    real(rk), pointer, dimension(:)                                 :: ice_area, heat_flux, snow_thick, t_ice
    !grids for temperature, salinity, coefficient of turbulence(kz2)
    real(rk), pointer, dimension(:, :)                              :: tem2, sal2, kz2
    real(rk), allocatable, dimension(:, :)                          :: cc, dcc, wbio !cc - main array
    real(rk), allocatable, dimension(:)                             :: bound_up, bound_low
    logical, allocatable, dimension(:)                              :: use_bound_up, use_bound_low
    character(len = attribute_length), allocatable, dimension(:)    :: par_name
    real(rk), allocatable, dimension(:)                             :: iz, density, pressure
    real(rk), allocatable, dimension(:)                             :: kz_tot
    real(rk)                                                        :: dt       !time step in [day/iteration]
    real(rk)                                                        :: io, wind_speed, pco2_atm, &
                                                                       width_bbl, resolution_bbl, width_bioturbation, &
                                                                       resolution_bioturbation, width_sediments, resolution_sediments
    type(ice_layer), pointer, dimension(:)  :: ice_l=>null()
    type(netcdf_algae_o), pointer           :: netcdf_ice=>null()
    type(netcdf_o), pointer                 :: netcdf_pelagic=>null()
    type(netcdf_o), pointer                 :: netcdf_bottom=>null()

    contains

    subroutine init_brom_transport()
    !REVISION HISTORY: Shamil Yakubov 11_2015 
    !- get_brom_par from brom.yaml was added
    !- io with fabm linking was deleted
    
    implicit none
    integer                 :: i
    
    !reading brom.yaml: module io_ascii
    call init_common()
    !getting variable values from brom.yaml: module io_ascii
    width_bbl               = get_brom_par("width_bbl")
    resolution_bbl          = get_brom_par("resolution_bbl")
    width_bioturbation      = get_brom_par("width_bioturbation")
    resolution_bioturbation = get_brom_par("resolution_bioturbation")
    width_sediments         = get_brom_par("width_sediments")
    resolution_sediments    = get_brom_par("resolution_sediments")
    year                    = get_brom_par("year")
    dt                      = get_brom_par("dt")
    !default environment in absence of actual data (could be retrievd from GOTM output)
    density                 = get_brom_par("density")   !kg m-3
    wind_speed              = get_brom_par("wind_speed")!m s-1
    pressure                = z + 10                    !dbar, roughly equivalent to depth in m+ 1 bar atmospheric pressure
    pco2_atm                = get_brom_par("pco2_atm")  !ppm
    !number of layers in ice
    number_of_layers        = get_brom_par("number_of_layers")
    
    !input_netcdf, kz2 - AKs, hice -  "time-averaged average ice thickness in cell", also generates grid according to input data
    call input_netcdf('KaraSea.nc', z, dz, kz_bio, lev_max, tem2, sal2, kz2, hice, boundary_bbl_sediments, &
        boundary_water_bbl, width_bbl, resolution_bbl, width_bioturbation, resolution_bioturbation, &
        width_sediments, resolution_sediments, year, ice_area, heat_flux, snow_thick, t_ice)
    kz2(1:35, :) = kz2(1:35, :) * 0.1 !temporarily for debugging
    kz2 = kz2 * dt * 86400 !because kz2 in m*m/sec, but here we need in day * time step
    
    !initialize FABM model from fabm.yaml
    call fabm_create_model_from_yaml_file(model)
    par_max = size(model%state_variables)
    call fabm_set_domain(model, lev_max)
    ! Specify vertical index of surface and bottom
    call model%set_surface_index(2)
    call model%set_bottom_index(boundary_bbl_sediments)
    !main array and its icrement allocating    
    allocate(cc(lev_max, par_max))
    allocate(dcc(lev_max, par_max))
    allocate(wbio(lev_max, par_max))    !vertical velocity (m/s, negative for sinking)
    !bounds allocating
    allocate(bound_up(par_max))
    allocate(bound_low(par_max))
    allocate(use_bound_up(par_max))
    allocate(use_bound_low(par_max))
    use_bound_up = .false.
    use_bound_low = .false.
    !auxiliary variables
    allocate(iz(lev_max))
    allocate(density(lev_max))
    allocate(pressure(lev_max))  
    
    ! Send pointers to state variable data to FABM
    do i = 1, par_max
        call fabm_link_bulk_state_data(model, i, cc(:, i))
    end do
    
    !provide initial array slices with temperature and salinity
    !these will be resent every time julianday is updated, below
    call fabm_link_bulk_data(model, standard_variables%temperature, tem2(:, 1))
    call fabm_link_bulk_data(model, standard_variables%practical_salinity, sal2(:, 1))
    call fabm_link_bulk_data(model, standard_variables%downwelling_photosynthetic_radiative_flux, iz)           !W m-2
    !additional environmental variables
    call fabm_link_bulk_data(model, standard_variables%density, density)                                        !kg m-3
    call fabm_link_bulk_data(model, standard_variables%pressure, pressure)                                      !dbar
    call fabm_link_horizontal_data(model, standard_variables%wind_speed, wind_speed)                            !m s-1
    call fabm_link_horizontal_data(model, standard_variables%mole_fraction_of_carbon_dioxide_in_air, pco2_atm)  !ppm
    call fabm_check_ready(model)
    
    allocate(par_name(par_max))
    do i = 1, par_max
        par_name(i) = model%state_variables(i)%name
    end do
    if (get_brom_par("is_input_data") == 0) then
        !allow FABM models to use their default initialization (this sets cc)
        call fabm_initialize_state(model, 1, lev_max)
    else
        !read initials from file ans save it inside cc massive (reset cc from input.dat)
        call porting_initial_state_variables_data(lev_max, par_name, cc)
    end if
    
    !ice algae initialization
    ice_l  => ice_layer(number_of_layers)
    !initializing output
    netcdf_ice => netcdf_algae_o()
    netcdf_pelagic => netcdf_o()
    netcdf_bottom  => netcdf_o()
    call netcdf_ice%init_netcdf_algae("output_a.nc", 1, number_of_layers)
    call netcdf_pelagic%init_netcdf("output_p.nc", 1, boundary_water_bbl - 1, model)
    call netcdf_bottom%init_netcdf("output_b.nc", boundary_water_bbl, lev_max, model)
    
    end subroutine init_brom_transport
    
    subroutine do_brom_transport()
    
    !REVISION HISTORY: Shamil Yakubov 16_12_2015 
    !- excess cirles were removed
    !- to the day/dt multiplation were moved to here
    !- 12_2015 caco3 on low boundary was added
    
    use calculate, only: calculate_phys, calculate_sed
    
    implicit none
    
    integer     :: i, id, ip, k, julianday, idt
    real(rk)    :: lat_light
    real(rk)    :: kc                       !attenuation constant for the self shading effect
    real(rk)    :: k_erlov                  !extinction coefficient
    real(rk)    :: so4, mn4, fe3, caco3     !for values
    integer     :: freq_az                  !vert.turb. / bhc frequency
    integer     :: freq_sed                 !sinking / bhc frequency
    integer     :: last_day, model_year = 0
    
    real(rk) :: flux_sf(size(model%state_variables))
    !indexes of state variables
    integer         :: i_O2, i_Mn2, i_Mn3, i_Mn4, i_H2S, i_Fe2, i_Fe3, i_FeS,   &
                       i_MnS, i_DON, i_PON, i_NH4, i_NO2, i_NO3, i_S0, i_S2O3,  &
                       i_SO4, i_DIC, i_Alk, i_PO4, i_Si, i_Sipart, i_Phy, i_Zoo,&
                       i_Baae, i_Bhae, i_Baan, i_Bhan, i_Hplus, i_CaCO3, i_FeS2, i_MnCO3
        
    !getting parameters from brom.yaml
    lat_light = get_brom_par("lat_light")
    kc = get_brom_par("kc")
    k_erlov = get_brom_par("k_erlov")
    freq_az = get_brom_par("freq_az")
    freq_sed  = get_brom_par("freq_sed ")
    last_day = get_brom_par("last_day")
    so4 = get_brom_par("bound(i_SO4)")
    mn4 = get_brom_par("bound_up(i_Mn4)")
    fe3 = get_brom_par("bound_up(i_Fe3)")
    caco3 = get_brom_par("bound_low(i_CaCO3)")
    
    i_O2 = find_index(par_name, 'niva_brom_redox_O2')              
    i_NO3 = find_index(par_name, 'niva_brom_redox_NO3')             
    i_NO2 = find_index(par_name, 'niva_brom_redox_NO2')            
    i_NH4 = find_index(par_name, 'niva_brom_redox_NH4')         
    i_DON = find_index(par_name, 'niva_brom_redox_DON')         
    i_PON = find_index(par_name, 'niva_brom_redox_PON')            
    i_PO4 = find_index(par_name, 'niva_brom_redox_PO4')      
    i_Si = find_index(par_name, 'niva_brom_redox_Si')      
    i_Sipart = find_index(par_name, 'niva_brom_redox_Sipart')          
    i_Mn4 = find_index(par_name, 'niva_brom_redox_Mn4')           
    i_Mn2 = find_index(par_name, 'niva_brom_redox_Mn2')          
    i_Mn3 = find_index(par_name, 'niva_brom_redox_Mn3')          
    i_MnS = find_index(par_name, 'niva_brom_redox_MnS')
    i_MnCO3 = find_index(par_name, 'niva_brom_redox_MnCO3')       
    i_Fe3 = find_index(par_name, 'niva_brom_redox_Fe3')          
    i_Fe2 = find_index(par_name, 'niva_brom_redox_Fe2')          
    i_FeS = find_index(par_name, 'niva_brom_redox_FeS')          
    i_SO4 = find_index(par_name, 'niva_brom_redox_SO4')          
    i_S2O3 = find_index(par_name, 'niva_brom_redox_S2O3')          
    i_S0 = find_index(par_name, 'niva_brom_redox_S0')           
    i_H2S = find_index(par_name, 'niva_brom_redox_H2S')         
    i_Phy = find_index(par_name, 'niva_brom_bio_Phy')          
    i_Zoo = find_index(par_name, 'niva_brom_bio_Het')          
    i_Baae = find_index(par_name, 'niva_brom_redox_Baae')           
    i_Bhae = find_index(par_name, 'niva_brom_redox_Bhae')          
    i_Baan = find_index(par_name, 'niva_brom_redox_Baan')          
    i_Bhan = find_index(par_name, 'niva_brom_redox_Bhan')          
    i_DIC = find_index(par_name, 'niva_brom_redox_DIC')          
    i_Alk = find_index(par_name, 'niva_brom_redox_Alk')          
    i_Hplus = find_index(par_name, 'niva_brom_redox_Hplus')          
    i_CaCO3 = find_index(par_name, 'niva_brom_redox_CaCO3')          
    i_FeS2 = find_index(par_name, 'niva_brom_redox_FeS2')
    
    idt = int(1. / dt)                               !number of cycles per day
    
    do  i = 0, last_day - 1                          !BIG Cycle ("i"-days)
        !it calculates julian day and year
        julianday = max(1, i - int(i / 365) * 365 + 1)
        if (julianday == 1) then
            model_year = model_year + 1
        end if
        
        !compute surface irradiance
        io = max(0., 80. * cos((lat_light - (23.5 * sin(2. * 3.14 * (julianday - 81.) /365.))) * 3.14 / 180.)) !W m-2
        
        !resend data that depend on julianday to FABM
        call fabm_link_bulk_data(model, standard_variables%temperature, tem2(:, julianday))
        call fabm_link_bulk_data(model, standard_variables%practical_salinity, sal2(:, julianday))
        
        !boudary conditions
        if (i_SO4 /= -1) then
            use_bound_up(i_SO4) = .true.
            bound_up(i_SO4) = so4
            use_bound_low(i_SO4) = .true.
            bound_low(i_SO4) = so4
        end if
        if (i_Mn4/=-1 .and. ice_area(julianday) < 0.5) then
            use_bound_up(i_Mn4) = .true.
            bound_up(i_Mn4) = mn4
        end if
        if (i_Fe3/=-1 .and. ice_area(julianday) < 0.5) then
            use_bound_up(i_Fe3) = .true.
            bound_up(i_Fe3) = fe3
        end if
        if (i_NO3 /= -1 .and. ice_area(julianday) < 0.5) then
            use_bound_up(i_NO3) = .true.
            bound_up(i_NO3) = 0. + (1. + sin(2 * 3.14 * (julianday - 115.) / 365.)) * 3.8! max 7.5 microM at day 205 approx.
        end if
        if (i_PO4 /= -1 .and. ice_area(julianday) < 0.5) then
            use_bound_up(i_PO4) = .true.
            bound_up(i_PO4) = 0. + (1. + sin(2 * 3.14 * (julianday - 115.) / 365.)) * 0.45! max 2 microM at day 205 approx.
        end if
        if (i_Si /= -1 .and. ice_area(julianday) < 0.5) then
            use_bound_up(i_Si)  = .true.
            bound_up(i_Si)  = 0. + (1. + sin(2 * 3.14 * (julianday - 115.) / 365.)) * 2.0! max 4 microM at day 205 approx.
        end if
        if (i_CaCO3 /= -1) then
            use_bound_low(i_CaCO3) = .true.
            bound_low(i_CaCO3) = caco3
        end if
        
        if(julianday == 249) then
            continue
        end if
        
        !ice algae processes calculated once per day is here, also recalculates io for bottom of ice layer
        call ice_l%do_slow_ice(t_ice(julianday), tem2(1, julianday), sal2(1, julianday), &
             hice(julianday), io, snow_thick(julianday), julianday, lat_light)
        
        !compute irradiance at depth
        do k = 1, lev_max
            !we check that we are not in the sediments
            !irradiance changing with depth
            if (k < boundary_bbl_sediments) then
                iz(k) = io * exp(-k_erlov * z(k))
            else
                iz(k) = 0.
            end if
        end do
        
        !timesteps in the course of a day
        do id = 1, idt
            
            !Kinetic processes (time integration is needed)
            dcc = 0.
            call fabm_do(model, 1, lev_max, dcc)
            dcc = dcc * dt * 86400.
            
            !integration
            cc = max(0.00000000001, (cc + dcc))
            
            !ice algae
            call ice_l%do_ice(cc(1, i_NH4), cc(1, i_NO2), cc(1, i_NO3), cc(1, i_PO4), dt, hice(julianday))
            
            do  ip = 1, freq_az
                !compute surface fluxes in FABM
                flux_sf = 0.
                if (hice(julianday) < 0.1) then
                    call fabm_do_surface(model,flux_sf)
                    flux_sf = flux_sf * dt * 86400.
                end if
                call calculate_phys(cc, lev_max, par_max, use_bound_up, use_bound_low, bound_up, bound_low, &
                            flux_sf, boundary_bbl_sediments, kz2, julianday, kz_bio, i_O2, dz, freq_az, dt)
            end do
            
            !sedimentation of particulate matter
            !compute residual vertical velocity (sinking/floating) in FABM.
            wbio = 0.
            call fabm_get_vertical_movement(model, 1, lev_max, wbio)
            wbio = wbio * dt * 86400.
            
            !freq_sed is defined in brom.yaml
            do  ip = 1, freq_sed
                call calculate_sed(par_max, lev_max, kz2, julianday, wbio, cc, dz, freq_sed, boundary_bbl_sediments)
            end do

            !for debug reason           
            if (any(isnan(cc))) then
                stop
            end if
        end do
       
        !-------NETCDF-----------------------------------------------------------------------------------------    
        write (*,'(a, i4, a, i4)') " model year:", model_year, "; julianday:", julianday
        call netcdf_ice%save_netcdf_algae(ice_l, 1, number_of_layers, number_of_layers, julianday)
        call netcdf_pelagic%save_netcdf(1, boundary_water_bbl - 1,  lev_max, julianday, cc, tem2, sal2, Kz2, model, z, iz)
        call netcdf_bottom%save_netcdf(boundary_water_bbl, lev_max, lev_max, julianday, cc, tem2, sal2, Kz2, model, z, iz)
        if (i == last_day - 1) then
            call netcdf_ice%close_netcdf_algae()
            call netcdf_pelagic%close_netcdf()
            call netcdf_bottom%close_netcdf()
            call saving_state_variables_data(model_year, julianday, lev_max, par_max, par_name, z, cc)
        end if
        !---END-of-NETCDF--------------------------------------------------------------------------------------
        
    end do
        
    end subroutine do_brom_transport
    
    subroutine clear_brom_transport()
    !allocated here
    deallocate(cc)
    deallocate(dcc)
    deallocate(wbio) 
    deallocate(bound_up)
    deallocate(bound_low)
    deallocate(use_bound_up)
    deallocate(use_bound_low)
    deallocate(par_name)
    deallocate(iz)
    deallocate(density)
    deallocate(pressure)
    !allocated in io_netcdf
    deallocate(tem2)
    deallocate(sal2)
    deallocate(hice)
    deallocate(Kz2)
    deallocate(z)
    deallocate(dz)
    deallocate(kz_bio)
    
    end subroutine clear_brom_transport
    
    end module brom_transport

!-----------------------------------------------------------------------
! Copyright under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------