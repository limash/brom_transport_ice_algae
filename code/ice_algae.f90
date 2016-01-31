    module ice_algae_lib
    
    !REVISION HISTORY:
    !Original author(s): Shamil Yakubov
    !to do: does k_ice depend on chl?
    !       z_conv equals 1 for now
    !       for nutreints in water the response should be implemented
    
    use fabm_types, only: rk
    
    implicit none
    private
    
    integer :: number_of_layers
    !length of the day
    real(rk):: day_length
    real(rk):: prev_ice_thickness
    real(rk):: ice_growth
    real(rk):: z_s = 0.03 !bottom layer depth
    real(rk):: p_b = 1.2 !maximum rate of photosynthesis (mg C mg Chl-1 h-1)
    real(rk):: chl_to_carbon = 30. !conversion factor (mg C mg Chl-1)
    real(rk):: carbon_to_oxy = 0.375 !conversion factor (mg C mg o2-1)
    public  :: ice_layer
    
    type ice_layer
        private
        !is_bottom
        logical:: is_bottom
        !depth of layers(z) and dz
        real(rk):: z
        real(rk):: dz
        !photosynthetic active radiation
        real(rk):: par_z
        !ice bulk variables
        real(rk):: bulk_temperature
        real(rk):: bulk_salinity
        real(rk):: bulk_density
        !ice brine variables
        real(rk):: brine_temperature
        real(rk):: brine_salinity
        real(rk):: brine_density
        !different algae nutrient concentrations
        real(rk):: a_carbon
        real(rk):: a_nitrogen
        real(rk):: a_phosphorus
        !increments of different algae nutrient concentrations
        real(rk):: d_a_carbon
        real(rk):: d_a_nitrogen
        real(rk):: d_a_phosphorus
        !brine concentrations of nutrients
        real(rk):: nh4
        real(rk):: d_nh4
        real(rk):: no2
        real(rk):: d_no2
        real(rk):: no3
        real(rk):: d_no3
        real(rk):: po4
        real(rk):: d_po4
        !volume of brine channels
        real(rk):: brine_relative_volume
        !some variables that will be used iteratively 
        real(rk):: last_gpp
        real(rk):: last_f_t
        real(rk):: last_v_ammonium
        real(rk):: last_mort
    contains
        private
        procedure, public:: do_slow_ice
        procedure, public:: do_ice
        procedure, public:: get_algae
        procedure:: do_grid
        procedure:: do_par
        procedure:: do_bulk_temperature
        procedure:: do_bulk_salinity
        procedure:: do_brine_salinity
        procedure:: do_brine_density
        procedure:: do_bulk_density
        procedure:: do_brine_relative_volume
        procedure:: do_nitrogen
        procedure:: do_phosphorus
        procedure:: do_bottom
        procedure:: do_congelation
        procedure:: brine_flux_z
        procedure:: brine_flux_s
        procedure:: diffusion_flux_s
        procedure:: congelation_flux_s
        procedure:: second_derivative
        procedure:: brine_release
        procedure:: do_a_carbon
        procedure:: gpp
        procedure:: f_par
        procedure:: f_s
        procedure:: f_t
        procedure:: f_nut
        procedure:: n_cell
        procedure:: p_cell
        procedure:: resp
        procedure:: exud
        procedure:: mort
        procedure:: melt
        procedure:: do_a_nitrogen
        procedure:: uptake_n
        procedure:: release_n
        procedure:: v_ammonium
        procedure:: v_nina
        procedure:: do_a_phosphorus
        procedure:: uptake_p
        procedure:: release_p
    end type ice_layer
    
    interface ice_layer
        procedure constructor_ice_layer
    end interface
    
    contains
    
    function constructor_ice_layer(number_of_layers_in)
    
    implicit none
    class(ice_layer), dimension(:), pointer:: constructor_ice_layer
    integer, intent(in):: number_of_layers_in
    
    number_of_layers = number_of_layers_in
    allocate(constructor_ice_layer(number_of_layers))

    constructor_ice_layer%is_bottom = .false.
    constructor_ice_layer(number_of_layers)%is_bottom = .true.
    constructor_ice_layer%a_carbon = 0.
    constructor_ice_layer%a_nitrogen = 0.
    constructor_ice_layer%a_phosphorus = 0.
    
    constructor_ice_layer%d_a_carbon = 0.
    constructor_ice_layer%d_a_nitrogen = 0.
    constructor_ice_layer%d_a_phosphorus = 0.
    
    constructor_ice_layer%nh4 = 0.
    constructor_ice_layer%no2 = 0.
    constructor_ice_layer%no3 = 0.
    constructor_ice_layer%po4 = 0.
    
    constructor_ice_layer%d_nh4 = 0.
    constructor_ice_layer%d_no2 = 0.
    constructor_ice_layer%d_no3 = 0.
    constructor_ice_layer%d_po4 = 0.
    
    prev_ice_thickness = 0.5 !only for first circle, 31 dec - 0.5m
    
    end function constructor_ice_layer
    
    subroutine do_slow_ice(self, air_temp, water_temp, water_sal, ice_thickness, &
                           io, snow_thick, julian_day, lat)
    !subroutine for variables should be calculated once per day
    
    implicit none
    class(ice_layer), dimension(:):: self
    real(rk), intent(in):: air_temp, water_temp, water_sal, ice_thickness
    real(rk), intent(in):: snow_thick, lat
    integer, intent(in) :: julian_day
    real(rk), intent(inout):: io
    real(rk)               :: foo
    
    if (ice_thickness < 0.2) then
        self%z = 0.
        self%dz = 0.
        self%par_z = 0.
        self%bulk_temperature = 0.
        self%bulk_salinity = 0.
        self%bulk_density = 0.
        self%brine_temperature = 0.
        self%brine_salinity = 0.
        self%brine_density = 0.
        self%a_carbon = 0.
        self%a_nitrogen = 0.
        self%a_phosphorus = 0.
        self%d_a_carbon = 0.
        self%d_a_nitrogen = 0.
        self%d_a_phosphorus = 0.
        self%nh4 = 0.
        self%d_nh4 = 0.
        self%no2 = 0.
        self%d_no2 = 0.
        self%no3 = 0.
        self%d_no3 = 0.
        self%po4 = 0.
        self%d_po4 = 0.
        self%brine_relative_volume = 0.
        return
    end if
    
    call self%do_grid(ice_thickness)
    
    call self%do_par(io, snow_thick)
    
    call self%do_bulk_temperature(air_temp, water_temp, ice_thickness)
    self%brine_temperature = self%bulk_temperature
    
    call self%do_bulk_salinity(ice_thickness, water_sal)
    call self%do_brine_salinity()
    
    call self%do_brine_density()
    call self%do_bulk_density()
    
    call self%do_brine_relative_volume()
    
    !here we will calculate the length of the day 
    !depends on latitude and julianday
    foo = cos((julian_day + 10.) * 2. * 3.14159 / 365.25)
    day_length = real(12. - 24./3.14159 * asin(tan(lat * 3.14159 / 180.) *&
                 tan(23.5 * 3.14159 / 180.) * foo))
    !to fix complex values
    if (isnan(day_length) .and. foo > 0) day_length = 0.
    if (isnan(day_length) .and. foo < 0) day_length = 24.
    !ice_growth calculation
    ice_growth = ice_thickness - prev_ice_thickness
    prev_ice_thickness = ice_thickness
        
    end subroutine do_slow_ice
    
    subroutine do_ice(self, nh4, no2, no3, po4, dt, ice_thickness)
    
    implicit none
    class(ice_layer), dimension(:):: self
    real(rk), intent(in):: nh4 !from upper water layer
    real(rk), intent(in):: no2 !from upper water layer
    real(rk), intent(in):: no3 !from upper water layer
    real(rk), intent(in):: po4 !from upper water layer
    real(rk), intent(in):: dt
    real(rk), intent(in):: ice_thickness
    
    if (ice_thickness < 0.2) then
        return
    end if
    
    !calculate fluxes of nutrients in ice
    call self%do_nitrogen(nh4, no2, no3)
    self%nh4 = self%nh4 + self%d_nh4 * dt
    self%no2 = self%no2 + self%d_no2 * dt
    self%no3 = self%no3 + self%d_no3 * dt
    self(1)%nh4 = self(2)%nh4
    self(1)%no2 = self(2)%no2
    self(1)%no3 = self(2)%no3
    self(1)%po4 = self(2)%po4
    
    call self%do_phosphorus(po4)
    self%po4 = self%po4 + self%d_po4 * dt
    self(1)%po4 = self(2)%po4
    
    !processes rates are per hour inside procedure
    call self%do_a_carbon()
    self%a_carbon = self%a_carbon + self%d_a_carbon * dt * 24.
    call self%do_a_nitrogen()
    self%a_nitrogen = self%a_nitrogen + self%d_a_nitrogen * dt * 24.
    call self%do_a_phosphorus()
    self%a_phosphorus = self%a_phosphorus + self%d_a_phosphorus * dt * 24.
        
    end subroutine do_ice
    
    subroutine get_algae(self, z, par_z, bulk_temperature, bulk_salinity, &
        bulk_density, brine_temperature, brine_salinity, brine_density, &
        a_carbon, a_nitrogen, a_phosphorus, nh4, no2, no3, po4, brine_relative_volume)
    
    implicit none
    class(ice_layer), dimension(:):: self
    real(rk), intent(out), dimension(number_of_layers):: z
    real(rk), intent(out), dimension(number_of_layers):: par_z
    real(rk), intent(out), dimension(number_of_layers):: bulk_temperature, bulk_salinity, bulk_density
    real(rk), intent(out), dimension(number_of_layers):: brine_temperature, brine_salinity, brine_density
    real(rk), intent(out), dimension(number_of_layers):: a_carbon, a_nitrogen, a_phosphorus
    real(rk), intent(out), dimension(number_of_layers):: nh4, no2, no3, po4
    real(rk), intent(out), dimension(number_of_layers):: brine_relative_volume
    
    z = self%z
    par_z = self%par_z
    bulk_temperature = self%bulk_temperature
    bulk_salinity = self%bulk_salinity
    bulk_density = self%bulk_density
    brine_temperature = self%brine_temperature
    brine_salinity = self%brine_salinity
    brine_density = self%brine_density
    a_carbon = self%a_carbon
    a_nitrogen = self%a_nitrogen
    a_phosphorus = self%a_phosphorus
    nh4 = self%nh4
    no2 = self%no2
    no3 = self%no3
    po4 = self%po4
    brine_relative_volume = self%brine_relative_volume
    
    end subroutine get_algae
    
    subroutine do_grid(self, hice)
    !it makes grid, bottom layer is on lower edge of ice and equals 3 cm
    !other layers depth = ice_thickness  - 3 cm / number of layers - 1
    
    implicit none
    class(ice_layer), dimension(:):: self
    real(rk), intent(in):: hice
    real(rk)            :: foo, bar
    integer             :: i
    
    foo = hice - z_s !z_s is for bottom layer depth (3cm)
    bar = number_of_layers - 1
    self(number_of_layers)%z = hice
    self(bar)%z = foo
    foo = foo / bar
    
    do i = (bar - 1), 1, -1
        self(i)%z = self(i+1)%z - foo
    end do
    
    do i = number_of_layers, 2, -1
        self(i)%dz = self(i)%z - self(i-1)%z
    end do
    self(1)%dz = foo
    
    end subroutine do_grid
    
    subroutine do_par(self, io, snow_thick)
    !io in Watts, to calculate it in micromoles photons per m2*s, [w] = 4.6*[micromole photons]
    
    implicit none
    class(ice_layer), dimension(:):: self
    real(rk), intent(in)    :: snow_thick
    real(rk), intent(inout) :: io
    real(rk)                :: io_e             !io in micromoles
    real(rk)                :: par_alb, par_scat
    real(rk), parameter     :: alb_ice = 0.744  !ice albedo
    real(rk), parameter     :: alb_snow = 0.9   !snow albedo
    real(rk), parameter     :: k_ice = 0.930    !extinction coeff. for ice m-1
    real(rk), parameter     :: k_snow = 4.3     !extinction coeff. for snow m-1
    real(rk), parameter     :: io_ice = 0.97    !fraction of rad transmitted through ice
    
    io_e = io / 4.6
    if (snow_thick <= 0.005) then !albedo influence
        par_alb = io * (1 - alb_ice)
    else
        par_alb = io * (1 - alb_snow) * exp(-k_snow * snow_thick)
    end if
    
    par_scat = par_alb * io_ice !after scattered surface of ice
    self%par_z = par_scat * exp(-k_ice * self%z)
    io = 4.6 * self(number_of_layers)%par_z
    
    end subroutine do_par
    
    subroutine do_bulk_temperature(self, air_temp, water_temp, ice_thickness)
    ![C]
    
    implicit none
    class(ice_layer), dimension(:):: self
    real(rk), intent(in):: air_temp, water_temp, ice_thickness
    integer:: i
    
    self%bulk_temperature = air_temp + ((water_temp - air_temp) * self%z) / ice_thickness
    do i = 1, number_of_layers
        if (self(i)%bulk_temperature > -0.2) self(i)%bulk_temperature = -0.2
    end do
    
    end subroutine do_bulk_temperature
    
    subroutine do_bulk_salinity(self, ice_thickness, water_sal)
    !bulk salinity through ice [ppt]
    
    implicit none
    class(ice_layer), dimension(:):: self
    real(rk), intent(in):: ice_thickness
    real(rk), intent(in):: water_sal
    integer:: i
    real(rk):: foo
    real(rk):: z_p
    
    foo = 0.
    do i = 1, number_of_layers
        foo = foo + self(i)%dz
        z_p = foo / ice_thickness
        self(i)%bulk_salinity = 19.539 * (z_p**(2)) - 19.93 * z_p + 8.913
    end do
    self(number_of_layers)%bulk_salinity = water_sal
    
    end subroutine do_bulk_salinity
    
    subroutine do_brine_salinity(self)
    !brine salinity through ice in [ppt]
    
    implicit none
    class(ice_layer), dimension(:):: self
    integer:: i
    
    do i = 1, number_of_layers
        if (self(i)%bulk_temperature < 0 .and. self(i)%bulk_temperature >=  -22.9) then
            self(i)%brine_salinity = -3.9921 + (-22.700   * self(i)%bulk_temperature) +&
                                               (-1.0015   * self(i)%bulk_temperature**2) +&
                                               (-0.019956 * self(i)%bulk_temperature**3)
        else if (self(i)%bulk_temperature < -22.9 .and. self(i)%bulk_temperature >=  -44) then
            self(i)%brine_salinity = 206.24 + (-1.8907    * self(i)%bulk_temperature) +&
                                              (-0.060868  * self(i)%bulk_temperature**2) +&
                                              (-0.0010247 * self(i)%bulk_temperature**3)
        else if (self(i)%bulk_temperature < -44) then
            self(i)%brine_salinity = -4442.1 + (-277.86  * self(i)%bulk_temperature) +&
                                               (-5.501   * self(i)%bulk_temperature**2) +&
                                               (-0.03669 * self(i)%bulk_temperature**3)
        end if
    end do
    
    self(number_of_layers)%brine_salinity = self(number_of_layers)%bulk_salinity
    
    end subroutine do_brine_salinity
    
    subroutine do_brine_density(self)
    !brine density through ice [g*m-3]
    
    implicit none
    class(ice_layer), dimension(:):: self
    real(rk):: c ![g*cm-1*ppt-1]
    
    c = 8e-4
    self%brine_density = (1. + c*self%brine_salinity)*1e6
    
    end subroutine do_brine_density
    
    subroutine do_bulk_density(self)
    !bulk density through ice [g*m-3]
    
    implicit none
    class(ice_layer), dimension(:):: self
    real(rk):: dens_pure !density of pure ice
    integer:: i
    real(rk):: foo
    
    dens_pure = 912000. ![g*m-3]
    
    do i = 1, number_of_layers
        foo = self(i)%brine_salinity
        if (foo < 2.) foo = 2.! function has spikes below x=1
        self(i)%bulk_density = dens_pure * self(i)%brine_density * foo /&
        (self(i)%brine_density * foo - self(i)%bulk_salinity *&
        (self(i)%brine_density - dens_pure))
    end do
    
    end subroutine do_bulk_density
    
    subroutine do_brine_relative_volume(self)
    !at the bottom layer it assumed equal to 0.5
    
    implicit none
    class(ice_layer), dimension(:):: self
    integer:: i
    real(rk):: foo
    
    do i = 1, number_of_layers - 1
        foo = self(i)%brine_salinity
        if (foo < 10.) foo = 10.!to fix above 1 values
        self(i)%brine_relative_volume = (self(i)%bulk_density * self(i)%bulk_salinity) /&
                                     (self(i)%brine_density  * foo)
    end do
    self(number_of_layers)%brine_relative_volume = 0.5
    
    end subroutine do_brine_relative_volume
    
    subroutine do_nitrogen(self, nh4, no2, no3)
    !calculate fluxes of nitrogen in brine channels
    !per day
    
    implicit none
    class(ice_layer), dimension(:):: self
    real(rk), intent(in):: nh4 !from upper water layer
    real(rk), intent(in):: no2 !from upper water layer
    real(rk), intent(in):: no3 !from upper water layer
    integer:: i
    
    !calculates d_nh4 for bottom layer
    call self%do_bottom(nh4, self(number_of_layers)%d_nh4, self(number_of_layers)%a_nitrogen, &
                        self%v_ammonium(number_of_layers), self%release_n(number_of_layers), &
                        self%nh4)
    call self%do_bottom(no2, self(number_of_layers)%d_no2, self(number_of_layers)%a_nitrogen, &
                        (self%v_nina(number_of_layers) / 2.), self%release_n(number_of_layers), &
                        self%no2)
    call self%do_bottom(no3, self(number_of_layers)%d_no3, self(number_of_layers)%a_nitrogen, &
                        (self%v_nina(number_of_layers) / 2.), self%release_n(number_of_layers), &
                        self%no3)
    do i = number_of_layers - 1, 2, -1
        call self%do_congelation(self(i)%d_nh4, self(i)%a_nitrogen, &
                        self%v_ammonium(i), self%release_n(i), i, &
                        self%nh4)
        call self%do_congelation(self(i)%d_no2, self(i)%a_nitrogen, &
                        (self%v_nina(number_of_layers) / 2.), self%release_n(i), i, &
                        self%no2)
        call self%do_congelation(self(i)%d_no3, self(i)%a_nitrogen, &
                        (self%v_nina(number_of_layers) / 2.), self%release_n(i), i, &
                        self%no3)
    end do
    
    end subroutine do_nitrogen
    
    subroutine do_phosphorus(self, po4)
    !calculate fluxes of phosphorus in brine channels
    !per day
    
    implicit none
    class(ice_layer), dimension(:):: self
    real(rk), intent(in):: po4
    integer:: i
    
    !calculates d_nh4 for bottom layer
    call self%do_bottom(po4, self(number_of_layers)%d_po4, self(number_of_layers)%a_phosphorus, &
                        (self%uptake_p(number_of_layers) * 24.), self%release_p(number_of_layers), &
                        self%po4)
    do i = number_of_layers - 1, 2, -1
        call self%do_congelation(self(i)%d_po4, self(i)%a_phosphorus, &
                        (self%uptake_p(i) * 24.), self%release_p(i), i, &
                        self%po4)
    end do
    
    end subroutine do_phosphorus
    
    subroutine do_bottom(self, nut_in_water, d_nut, a_nut, uptake, release, nut_brine)
    !calculate bottom fluxes
    !per day
    
    implicit none
    class(ice_layer), dimension(:):: self
    real(rk), dimension(:)        :: nut_brine
    real(rk), intent(in) :: nut_in_water, uptake, release
    real(rk), intent(in) :: a_nut !in algae nutrients
    real(rk), intent(out):: d_nut !increment of calculated brine nutrient
    
    d_nut = self%brine_flux_z(number_of_layers, nut_brine, nut_in_water) +&
            self%brine_flux_s(number_of_layers, nut_brine, nut_in_water) +&
            self%diffusion_flux_s(number_of_layers, nut_brine, nut_in_water) +&
            self%congelation_flux_s(number_of_layers, nut_brine, nut_in_water) -&
            a_nut * (uptake - release)
    
    end subroutine do_bottom
    
    subroutine do_congelation(self, d_nut, a_nut, uptake, release, layer, nut_brine)
    !brine concentrarion changes
    !per day
    
    implicit none
    class(ice_layer), dimension(:):: self
    real(rk), dimension(:)        :: nut_brine
    integer              :: layer
    real(rk), intent(in) :: uptake, release
    real(rk), intent(in) :: a_nut !in algae nutrients
    real(rk), intent(out):: d_nut
    
    d_nut = self%brine_flux_z(layer, nut_brine) -&
            a_nut * (uptake - release)
    
    end subroutine do_congelation
    
    function brine_flux_z(self, layer, nut_brine, nut_in_water)
    !brine flux for all layers - per day
    
    implicit none
    class(ice_layer), dimension(:)  :: self
    real(rk), dimension(:)          :: nut_brine
    integer                         :: layer
    real(rk)                        :: brine_flux_z
    real(rk), intent(in), optional  :: nut_in_water !for bottom layer
    real(rk)                        :: gravity_drainage = 1.e-8! [m * s-1]
    real(rk)                        :: z_conv = 1. !vertical distance over which sea ice is influenced by gravity drainage
    
    if (present(nut_in_water)) then
        brine_flux_z = self(layer)%brine_relative_volume *&
        gravity_drainage * z_conv *&
        self%second_derivative(nut_brine(layer - 1), nut_brine(layer), nut_in_water, self(layer - 1)%dz, self(layer)%dz)
    else
        brine_flux_z = self(layer)%brine_relative_volume *&
        gravity_drainage * z_conv *&
        self%second_derivative(nut_brine(layer - 1), nut_brine(layer), nut_brine(layer + 1), self(layer - 1)%dz, self(layer)%dz)
    end if
    
    brine_flux_z = 86400. * brine_flux_z !convert per second to per day
    
    end function brine_flux_z
    
    function brine_flux_s(self, layer, nut_brine, nut_in_water)
    !brine flux for bottom - per day
    
    implicit none
    class(ice_layer), dimension(:)  :: self
    real(rk), dimension(:)          :: nut_brine
    integer                         :: layer
    real(rk)                        :: brine_flux_s
    real(rk), intent(in)            :: nut_in_water !for bottom layer
    
    brine_flux_s = self(layer)%brine_relative_volume *&
        self%brine_release() * z_s *&
        self%second_derivative(nut_brine(layer - 1), nut_brine(layer), nut_in_water, self(layer - 1)%dz, self(layer)%dz)
    brine_flux_s = 86400. * brine_flux_s !convert per second to per day
    
    end function brine_flux_s
    
    function diffusion_flux_s(self, layer, nut_brine, nut_in_water)
    !brine diffusion flux for bottom - per day
    
    implicit none
    class(ice_layer), dimension(:)  :: self
    real(rk), dimension(:)          :: nut_brine
    integer                         :: layer
    real(rk)                        :: diffusion_flux_s
    real(rk), intent(in)            :: nut_in_water !for bottom layer
    real(rk)                        :: k_wi = 1.e-5! [m2 s-1] diffusion coef. between bottom layer and sea water
    
    diffusion_flux_s = self(layer)%brine_relative_volume *&
        k_wi *&
        self%second_derivative(nut_brine(layer - 1), nut_brine(layer), nut_in_water, self(layer - 1)%dz, self(layer)%dz)
    diffusion_flux_s = 86400. * diffusion_flux_s !convert per second to per day
    
    end function diffusion_flux_s
    
    function congelation_flux_s(self, layer, nut_brine, nut_in_water)
    !brine diffusion flux for bottom - per day
    
    implicit none
    class(ice_layer), dimension(:)  :: self
    real(rk), dimension(:)          :: nut_brine
    integer                         :: layer
    real(rk)                        :: congelation_flux_s
    real(rk), intent(in)            :: nut_in_water !for bottom layer
    
    congelation_flux_s = (ice_growth * (nut_in_water - nut_brine(layer)) * nut_in_water) /&
                         ((z_s + ice_growth) * self(layer)%brine_relative_volume * 86400.)
    
    end function congelation_flux_s
    
    function second_derivative(self, pre_layer, layer, next_layer, pre_h, next_h)
    !second_derivative
    
    implicit none
    class(ice_layer), dimension(:)  :: self
    real(rk), intent(in)            :: pre_layer, layer, next_layer, pre_h, next_h
    real(rk)                        :: second_derivative
    real(rk)                        :: h
    
    h = (pre_h + next_h) / 2.
    second_derivative = (pre_layer - 2. * layer + next_layer) / h**2
    
    end function second_derivative
    
    function brine_release(self)
    !brine_release [m s-1]
    
    implicit none
    class(ice_layer), dimension(:)  :: self
    real(rk)                        :: brine_release
    
    brine_release = (9.667e-9 + 4.49e-6 * ice_growth - 1.39e-7 * ice_growth**2) / 100.
    
    end function brine_release
    
    subroutine do_a_carbon(self)
    !brine concentrations changes of ice algae carbon
    !mg C per hour
    
    implicit none
    class(ice_layer), dimension(:):: self
    integer:: i
    real(rk):: recruit = 0.01 !recruitment of algae on bottom

    do i = number_of_layers, 1, -1
        if (self(i)%is_bottom == .true.) then
            self(i)%d_a_carbon = self(i)%a_carbon * (self%gpp(i)- self%resp(i) -&
                      self%exud(i) - self%mort(i) - self%melt()) + recruit
        else
            self(i)%d_a_carbon = self(i)%a_carbon * (self%gpp(i)- self%resp(i) -&
                      self%exud(i) - self%mort(i))
        end if
    end do
    
    end subroutine do_a_carbon
    
    function gpp(self, layer)
    !gross photosynthetic production - per hour
    
    implicit none
    class(ice_layer), dimension(:):: self
    integer :: layer
    real(rk):: gpp
    
    gpp = p_b * self%f_par(layer) * self%f_s(layer) * self%f_t(layer) *&
          self%f_nut(layer) / chl_to_carbon
    self(layer)%last_gpp = gpp
    
    end function gpp
    
    function f_par(self, layer)
    !light limitation, dimensionless
    implicit none
    class(ice_layer), dimension(:):: self
    integer :: layer
    real(rk):: f_par
    !initial slope or photosynthetic efficiency
    real(rk), parameter:: alpha = 0.227 ![mg C mg Chl-1 microE m2s]
    !degree of photoinhibition
    real(rk), parameter:: betta = 0. ![mg C mg Chl-1 h-1(microM photons m-2 s-1)-1]
    
    f_par = (1 - exp(-1 * alpha * self(layer)%par_z / p_b)) &
               * exp(-1 * betta * self(layer)%par_z / p_b)
    
    end function f_par
    
    function f_s(self, layer)
    !salinity limitation, dimensionless
    implicit none
    class(ice_layer), dimension(:):: self
    integer :: layer
    real(rk):: f_s
    
    f_s = 1.1e-2 + 3.012e-2  * self(layer)%brine_salinity    +&
                   1.0342e-3 * self(layer)%brine_salinity**2 +&
                   4.6033e-5 * self(layer)%brine_salinity**3 +&
                   4.926e-7  * self(layer)%brine_salinity**4 +&
                   1.659e-9  * self(layer)%brine_salinity**5
    
    end function f_s
    
    function f_t(self, layer)
    !temperature limitation, dimensionless
    implicit none
    class(ice_layer), dimension(:):: self
    integer :: layer
    real(rk):: f_t
    !reference temperature
    real(rk):: t_0 = 0.
    !temperature augmentation rate
    real(rk):: temp_aug_rate = 0.0663
    
    f_t = exp(temp_aug_rate * (self(layer)%brine_temperature - t_0))
    self(layer)%last_f_t = f_t
    
    end function f_t
    
    function f_nut(self, layer)
    !nutrient limitation, dimensionless
    implicit none
    class(ice_layer), dimension(:):: self
    integer :: layer
    real(rk):: f_nut
    real(rk):: n_cell
    real(rk):: p_cell
    !half saturation constant for growth limited by nitrogen cell quota
    real(rk), parameter:: kn_cell = 0.028 ![mg N mg C-1]
    !half saturation constant for growth limited by phosphorus cell quota
    real(rk), parameter:: kp_cell = 0.004 ![mg P mg C-1]
    
    n_cell = self%n_cell(layer)
    p_cell = self%p_cell(layer)
    
    if (self(layer)%a_carbon == 0.) then
        f_nut = 0.
    else
        f_nut = min(n_cell / (kn_cell + n_cell), p_cell / (kp_cell + p_cell))
    end if
    
    end function f_nut
    
    function n_cell(self, layer)
    !nitrogen cell quota - [mg N mg C-1]
    implicit none
    class(ice_layer), dimension(:):: self
    integer :: layer
    real(rk):: n_cell
    
    n_cell = self(layer)%a_nitrogen / self(layer)%a_carbon
    
    end function n_cell
    
    function p_cell(self, layer)
    !phosphorus cell quota - [mg P mg C-1]
    implicit none
    class(ice_layer), dimension(:):: self
    integer :: layer
    real(rk):: p_cell
    
    p_cell = self(layer)%a_phosphorus / self(layer)%a_carbon
    
    end function p_cell
    
    function resp(self, layer)
    !respiration - per hour
    
    implicit none
    class(ice_layer), dimension(:):: self
    integer :: layer
    real(rk):: resp
    !maintenance respiration
    real(rk), parameter:: r0 = 0.01 ![mmol o2 mg chl-1 h-1]
    !linear coefficient of increase in biomass-specific dark respiration
    !with gross photosynthesis
    real(rk), parameter:: r_dark = 0.3
    !ratio between respiration in the light and
    !respiration in the dark (dimensionless)
    real(rk), parameter:: dl_ratio = 2.
    real(rk):: resp_day, resp_night
    
    resp_day = (r0 + r_dark * dl_ratio * self(layer)%last_gpp) *&
        (carbon_to_oxy / chl_to_carbon) *&
        16. * day_length * self(layer)%last_f_t
    
    resp_night = (r0 + r_dark * self(layer)%last_gpp) *&
        (carbon_to_oxy / chl_to_carbon) *&
        16. * (24. - day_length) * self(layer)%last_f_t
    
    !to make resp coeffic in per hour units
    resp = (resp_day + resp_night) / 24.
    
    end function resp
    
    function exud(self, layer)
    !ice algae exudation rate - per hour
    implicit none
    class(ice_layer), dimension(:):: self
    integer :: layer
    real(rk):: exud
    !fraction exudated
    real(rk), parameter:: exud_rate = 0.1
    
    exud = exud_rate * self(layer)%last_gpp
    
    end function exud
    
    function mort(self, layer)
    !mortality loss - per hour
    implicit none
    class(ice_layer), dimension(:):: self
    integer :: layer
    real(rk):: mort
    !algae mortality at temperature 0 C
    real(rk), parameter:: phy_mort = 9.23e-5 ![h-1]
    
    mort = phy_mort * self(layer)%last_f_t
    self(layer)%last_mort = mort
    
    end function mort
    
    function melt(self)
    !melting rate - per hour
    implicit none
    class(ice_layer), dimension(:):: self
    real(rk):: melt
    
    melt =  -1. * ice_growth / z_s / 24. !24 is to transform in per hour units
    if (melt < 0) melt = 0.
    
    end function melt
    
    subroutine do_a_nitrogen(self)
    !brine concentrations changes of ice algae nitrogen
    !mg N per hour
    
    implicit none
    class(ice_layer), dimension(:):: self
    integer:: i
    real(rk):: recruit = 0.01 !recruitment of algae on bottom

    do i = number_of_layers, 1, -1
        if (self(i)%is_bottom == .true.) then
            self(i)%d_a_nitrogen = self(i)%a_nitrogen * (self%uptake_n(i) - self%release_n(i) -&
                      self(i)%last_mort - self%melt()) + recruit
        else
            self(i)%d_a_nitrogen = self(i)%a_nitrogen * (self%uptake_n(i) - self%release_n(i) -&
                      self(i)%last_mort)
        end if
    end do
    
    end subroutine do_a_nitrogen
    
    function uptake_n(self, layer)
    !a lot of problems could be here - uptake of N by algae - per hour
    implicit none
    class(ice_layer), dimension(:):: self
    integer :: layer
    real(rk):: uptake_n
    
    uptake_n = (self%v_ammonium(layer) + self%v_nina(layer)) / 24. !to per hour transform
    
    end function uptake_n
    
    function v_ammonium(self, layer)
    !ammonium uptake rate - per day
    implicit none
    class(ice_layer), dimension(:):: self
    integer :: layer
    real(rk):: v_ammonium
    real(rk):: v_max_n = 1.08! [1/day] maximal uptake rate
    real(rk):: k_ammonium = 2.94! [microM N] half saturation constant for ammonium uptake
    real(rk):: n_max = 0.53! [mg N * mg C-1]
    real(rk):: n_min = 0.1! [mg N * mg C-1]
    real(rk):: max_n_p = 291! max N/P in algae
    real(rk):: foo
    
    foo = self(layer)%a_nitrogen / self(layer)%a_phosphorus
    if (self(layer)%a_nitrogen > n_min .and. self(layer)%a_nitrogen < n_max &
        .and. foo < max_n_p) then
        v_ammonium = (v_max_n * self(layer)%nh4) / (k_ammonium + self(layer)%nh4) *&
            (1 - (self%n_cell(layer) / n_max))
    else
        v_ammonium = 0.
    end if
    self(layer)%last_v_ammonium = v_ammonium
    
    end function v_ammonium
    
    function v_nina(self, layer)
    !nitrite and nitrate uptake rate - per day
    implicit none
    class(ice_layer), dimension(:):: self
    integer :: layer
    real(rk):: v_nina
    real(rk):: v_max_n = 1.08! [1/day] maximal uptake rate
    real(rk):: k_nina = 30! [microM N] half saturation constant for no2 and no3 uptake
    real(rk):: n_max = 0.53! [mg N * mg C-1] maximum nitrogen cell quota
    real(rk):: n_min = 0.1! [mg N * mg C-1]
    real(rk):: max_n_p = 291! max N/P in algae
    real(rk):: foo
    
    foo = self(layer)%a_nitrogen / self(layer)%a_phosphorus
    if (self(layer)%a_nitrogen > n_min .and. self(layer)%a_nitrogen < n_max &
        .and. foo < max_n_p) then
        v_nina = max(0., (v_max_n - self(layer)%last_v_ammonium)) *&
            (self(layer)%no2 + self(layer)%no3) / (k_nina + self(layer)%no2 + self(layer)%no3) *&
            (1 - (self%n_cell(layer) / n_max))
    else
        v_nina = 0.
    end if
    
    end function v_nina
    
    function release_n(self, layer)
    !release - release everything more then n_max - mg N
    
    implicit none
    class(ice_layer), dimension(:):: self
    integer :: layer
    real(rk):: release_n
    real(rk):: n_max = 0.53! [mg N * mg C-1] maximum nitrogen cell quota
    
    if (self%n_cell(layer) > n_max .and. self(layer)%a_nitrogen > 0.) then
        release_n = (self(layer)%a_nitrogen - n_max * self(layer)%a_carbon) / self(layer)%a_nitrogen
    else
        release_n = 0.
    end if
    
    end function release_n
    
    subroutine do_a_phosphorus(self)
    !brine concentrations changes of ice algae po4
    !mg N per hour
    
    implicit none
    class(ice_layer), dimension(:):: self
    integer:: i
    real(rk):: recruit = 0.01 !recruitment of algae on bottom

    do i = number_of_layers, 1, -1
        if (self(i)%is_bottom == .true.) then
            self(i)%d_a_phosphorus = self(i)%a_phosphorus * (self%uptake_p(i) - self%release_p(i) -&
                      self(i)%last_mort - self%melt()) + recruit
        else
            self(i)%d_a_phosphorus = self(i)%a_phosphorus * (self%uptake_p(i) - self%release_p(i) -&
                      self(i)%last_mort)
        end if
    end do
    
    end subroutine do_a_phosphorus
    
    function uptake_p(self, layer)
    !a lot of problems could be here - uptake of N by algae - per hour
    implicit none
    class(ice_layer), dimension(:):: self
    integer :: layer
    real(rk):: uptake_p
    real(rk):: phos_min = 0.002! [mg P mg C-1] minimum phosphorus cell quota
    real(rk):: phos_max = 0.08! [mg P mg C-1] maximum phosphorus cell quota
    real(rk):: min_n_p = 4.! min N/P in algae
    real(rk):: v_max_p = 1.08! [d-1] maximal uptake rate of po4
    real(rk):: k_p = 2.! [microM P] half saturation constant for po4 uptake
    real(rk):: foo
    
    foo = self(layer)%a_nitrogen / self(layer)%a_phosphorus
    if (self(layer)%a_phosphorus > phos_min .and. self(layer)%a_phosphorus < phos_max &
        .and. foo > min_n_p) then
        uptake_p = (v_max_p * self(layer)%po4) / (k_p + self(layer)%po4) *&
            (1 - (self%p_cell(layer) / phos_max))
    else
        uptake_p = 0.
    end if
    uptake_p = uptake_p / 24. !to per hour transform
        
    end function uptake_p
    
    function release_p(self, layer)
    !release - release everything more then n_max - mg N
    
    implicit none
    class(ice_layer), dimension(:):: self
    integer :: layer
    real(rk):: release_p
    real(rk):: phos_max = 0.08! [mg P mg C-1] maximum phosphorus cell quota
    
    if (self%p_cell(layer) > phos_max .and. self(layer)%a_phosphorus > 0.) then
        release_p = (self(layer)%a_phosphorus - phos_max * self(layer)%a_carbon) / self(layer)%a_phosphorus
    else
        release_p = 0.
    end if
    
    end function release_p
    
    end module ice_algae_lib
    
!-----------------------------------------------------------------------
! Copyright under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------