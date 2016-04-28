#define NUMBER_OF_LAYERS_IN_ICE 5

module ice_algae_lib
!to do: does k_ice depend on chl? - update at every model
!step implementation
!     f_par integration in layers with algae for more
!precise light limitation
!     z_conv equals 1 for now - correction needed
!     for nutreints in water the response should be implemented
!     recruitment of algae - melnikov
!     Si, detritus and oxygen implementation!
!     parts of nitrite and nitrates of v_nina
    
    use fabm_types, only: rk
    
    implicit none
    private
    
    integer :: number_of_layers = NUMBER_OF_LAYERS_IN_ICE 

    real(rk):: day_length
    real(rk):: prev_ice_thickness
    real(rk):: ice_growth
    real(rk):: ice_growth_temp
    real(rk):: a_b !algae position
    logical:: trigger !ice control
    logical:: trigger_melting !melting control

    real(rk):: alb_ice = 0.744  !ice albedo
    real(rk):: alb_snow = 0.9   !snow albedo
    real(rk):: k_ice = 0.930    !extinction coeff. for ice m-1
    real(rk):: k_snow = 4.3     !extinction coeff. for snow m-1
    real(rk):: io_ice_par = 0.97    !fraction of rad transmitted through ice

    real(rk):: z_s = 0.03 !bottom layer width
    !maximum rate of photosynthesis (mg C mg Chl-1 h-1)
    real(rk):: p_b = 1.2
    real(rk):: chl_to_carbon = 30. !conversion factor (mg C mg Chl-1)
    real(rk):: carbon_to_oxy = 0.375 !conversion factor (mg C mg o2-1)
    
    !vertical distance over which sea ice is influenced by gravity drainage
    real(rk):: z_conv = 1.        
    real(rk):: gravity_drainage = 1.e-8 ![m * s-1]
    
    ![m2 s-1] diffusion coef. between bottom layer and sea water
    real(rk):: k_wi = 1.e-5
    real(rk):: v_max_n = 1.08![1/day] maximal uptake rate of nitrogen
    ![microM N] half saturation constant for ammonium uptake
    real(rk):: k_ammonium = 2.94
    real(rk):: n_max = 0.53![mg N * mg C-1]
    real(rk):: n_min = 0.1![mg N * mg C-1]
    real(rk):: max_n_p = 291!max N/P in algae
    ![microM N] half saturation constant for no2 and no3 uptake
    real(rk):: k_nina = 30

    real(rk):: phos_min = 0.002![mg P mg C-1] minimum phosphorus cell quota
    real(rk):: phos_max = 0.08![mg P mg C-1] maximum phosphorus cell quota
    real(rk):: min_n_p = 4.!min N/P in algae
    real(rk):: v_max_p = 1.08![d-1] maximal uptake rate of po4
    real(rk):: k_p = 2.![microM P] half saturation constant for po4 uptake

    real(rk):: recruit = 0.01!recruitment of algae on bottom

    !initial slope or photosynthetic efficiency
    real(rk):: alpha = 0.227![mg C mg Chl-1 microE m2s]
    !degree of photoinhibition
    real(rk):: betta = 0.![mg C mg Chl-1 h-1(microM photons m-2 s-1)-1]

    !half saturation constant for growth limited by nitrogen cell quota
    real(rk):: kn_cell = 0.028![mg N mg C-1]
    !half saturation constant for growth limited by phosphorus cell quota
    real(rk):: kp_cell = 0.004![mg P mg C-1]

    !reference temperature
    real(rk):: t_0 = 0.
    !temperature augmentation rate
    real(rk):: temp_aug_rate = 0.0663
    
    !maintenance respiration
    real(rk):: r0 = 0.01 ![mmol o2 mg chl-1 h-1]
    !linear coefficient of increase in biomass-specific dark respiration
    !with gross photosynthesis
    real(rk):: r_dark = 0.3
    !ratio between respiration in the light and
    !respiration in the dark (dimensionless)
    real(rk):: dl_ratio = 2.
        
    !fraction exudated
    real(rk):: exud_rate = 0.1

    !algae mortality at temperature 0 C
    real(rk):: phy_mort = 9.23e-5 ![h-1]

    real(rk), dimension(NUMBER_OF_LAYERS_IN_ICE):: nh4_m, no2_m, no3_m
    real(rk), dimension(NUMBER_OF_LAYERS_IN_ICE):: po4_m
    real(rk), dimension(NUMBER_OF_LAYERS_IN_ICE):: dz_m, z_m
    real(rk), dimension(NUMBER_OF_LAYERS_IN_ICE):: a_carbon_m
    real(rk), dimension(NUMBER_OF_LAYERS_IN_ICE):: a_nitrogen_m
    real(rk), dimension(NUMBER_OF_LAYERS_IN_ICE):: a_phosphorus_m

    public :: ice_layer
    
    type ice_layer
        private
        !depth of layers(z)
        real(rk):: z
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
        procedure, public:: rewrite_algae
        procedure, public:: do_rec_algae
        procedure:: do_grid
        procedure:: do_par
        procedure:: do_congelation_algae
        procedure:: do_melting_algae
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
    
    function constructor_ice_layer(number_of_layers_out)
    
        type(ice_layer), dimension(:), pointer:: constructor_ice_layer
        integer, intent(out):: number_of_layers_out
        
        number_of_layers_out = number_of_layers
        allocate(constructor_ice_layer(number_of_layers))
        nh4_m = 0.
        no2_m = 0.
        no3_m = 0.
        po4_m = 0.
        dz_m = 0.
        z_m = 0.
        a_carbon_m = 0.
        a_nitrogen_m = 0.
        a_phosphorus_m = 0.

        constructor_ice_layer%z = 0.
        constructor_ice_layer%par_z = 0.
        constructor_ice_layer%bulk_temperature = 0.
        constructor_ice_layer%bulk_salinity = 0.
        constructor_ice_layer%bulk_density = 0.
        constructor_ice_layer%brine_temperature = 0.
        constructor_ice_layer%brine_salinity = 0.
        constructor_ice_layer%brine_density = 0.

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
        
        constructor_ice_layer%brine_relative_volume = 0.
        
        constructor_ice_layer%last_gpp = 0.
        constructor_ice_layer%last_f_t = 0.
        constructor_ice_layer%last_v_ammonium = 0.
        constructor_ice_layer%last_mort = 0.

        prev_ice_thickness = 0.5 !only for first circle, 31 dec - 0.5m
        a_b = 0.

        trigger = .false.
        trigger_melting = .false.

    end function constructor_ice_layer

    subroutine do_slow_ice(self, lvl, air_temp, water_temp, water_sal, &
            ice_thickness, io, io_ice, snow_thick, julian_day, lat)
    !subroutine for variables should be calculated once per day
        
        class(ice_layer):: self
        integer,  intent(in):: lvl
        real(rk), intent(in):: air_temp, water_temp
        real(rk), intent(in):: water_sal, ice_thickness
        real(rk), intent(in):: snow_thick, lat
        integer,  intent(in):: julian_day
        real(rk), intent(in):: io
        real(rk), intent(out)  :: io_ice
        real(rk)               :: foo
        
        if (ice_thickness < 0.2) then
            self%z = 0.
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

            self%last_v_ammonium = 0.
            self%last_gpp = 0.
            self%last_f_t = 0.
            self%last_mort = 0.

            nh4_m = 0.
            no2_m = 0.
            no3_m = 0.
            po4_m = 0.
            dz_m = 0.
            z_m = 0.
            a_carbon_m = 0.
            a_nitrogen_m = 0.
            a_phosphorus_m = 0.

            a_b = 0.
            ice_growth = 0.
            prev_ice_thickness = 0.

            io_ice = io

            return
        end if
    
        call self%do_grid(lvl, ice_thickness)
        
        call self%do_par(lvl, io, io_ice, snow_thick)
        
        call self%do_bulk_temperature(air_temp, water_temp, ice_thickness)
        self%brine_temperature = self%bulk_temperature
        
        call self%do_bulk_salinity(lvl, ice_thickness, water_sal)
        call self%do_brine_salinity(lvl)
        
        call self%do_brine_density()
        call self%do_bulk_density()
        
        call self%do_brine_relative_volume(lvl)
        
        !calculate the length of the day 
        !depends on latitude and julianday
        foo = cos((julian_day + 10.) * 2. * 3.14159 / 365.25)
        day_length = real(12. - 24./3.14159 * &
            asin(tan(lat * 3.14159 / 180.) * &
            tan(23.5 * 3.14159 / 180.) * foo))
        !to fix complex values
        if (isnan(day_length) .and. foo > 0) day_length = 0.
        if (isnan(day_length) .and. foo < 0) day_length = 24.
        
    end subroutine do_slow_ice
    
    subroutine do_rec_algae(self, lvl, ice_thickness, da_c, before)
    
        class(ice_layer):: self
        integer,  intent(in):: lvl
        real(rk), intent(in):: ice_thickness
        real(rk), intent(out)  :: da_c
        logical, intent(in)    :: before

        !ice_growth/melting calculation
        if ((trigger .eqv. .false.) .and. (before .eqv. .true.)) then
            ice_growth = ice_thickness - prev_ice_thickness
            prev_ice_thickness = ice_thickness
            trigger = .true.
        end if
        ice_growth_temp = ice_growth

        da_c = 0.
        if ((ice_growth_temp > 0.) .and. (before .eqv. .false.)) then
            call self%do_congelation_algae(lvl, ice_growth_temp)
            if (lvl == 1) trigger = .false.
        else if ((ice_growth_temp <= 0.) .and. (before .eqv. .true.)) then
            call self%do_melting_algae(lvl, ice_growth_temp, da_c)
            if (lvl == 1) trigger_melting = .false.
            if (lvl == 1) trigger = .false.
        end if

    end subroutine do_rec_algae
    
    subroutine do_ice(self, lvl, nh4, no2, no3, po4, dt, ice_thickness)
    
        class(ice_layer):: self
        integer,  intent(in):: lvl
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
        !process rates are per day inside procedure
        call self%do_nitrogen(lvl, nh4, no2, no3)

        self%nh4 = self%nh4 + self%d_nh4 * dt
        nh4_m(lvl) = nh4_m(lvl) + self%d_nh4 * dt

        self%no2 = self%no2 + self%d_no2 * dt
        no2_m(lvl) = no2_m(lvl) + self%d_no2 * dt

        self%no3 = self%no3 + self%d_no3 * dt
        no3_m(lvl) = no3_m(lvl) + self%d_no3 * dt

        !process rates are per day inside procedure
        call self%do_phosphorus(lvl, po4)
        self%po4 = self%po4 + self%d_po4 * dt
        po4_m(lvl) = po4_m(lvl) + self%d_po4 * dt
        
        if (lvl == 2) then
            nh4_m(lvl - 1) = nh4_m(lvl)
            no2_m(lvl - 1) = no2_m(lvl)
            no3_m(lvl - 1) = no3_m(lvl)
            po4_m(lvl - 1) = po4_m(lvl)
        end if
        
        !process rates are per hour inside procedure
        call self%do_a_carbon(lvl)
        self%a_carbon = self%a_carbon + self%d_a_carbon * dt * 24.
        call self%do_a_nitrogen(lvl)
        self%a_nitrogen = self%a_nitrogen + self%d_a_nitrogen * dt * 24.
        call self%do_a_phosphorus(lvl)
        self%a_phosphorus = self%a_phosphorus + self%d_a_phosphorus * dt * 24.
        
    end subroutine do_ice
    
    subroutine get_algae(self, z, par_z, bulk_temperature, bulk_salinity, &
        bulk_density, brine_temperature, brine_salinity, brine_density, &
        a_carbon, a_nitrogen, a_phosphorus, nh4, no2, no3, po4, brine_relative_volume)
    
        class(ice_layer):: self
        real(rk), intent(out):: z
        real(rk), intent(out):: par_z
        real(rk), intent(out):: bulk_temperature, bulk_salinity, bulk_density
        real(rk), intent(out):: brine_temperature, brine_salinity, brine_density
        real(rk), intent(out):: a_carbon, a_nitrogen, a_phosphorus
        real(rk), intent(out):: nh4, no2, no3, po4
        real(rk), intent(out):: brine_relative_volume
        
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

    subroutine rewrite_algae(self, lvl, befor_after)
    !rewrites algae layers

        class(ice_layer):: self
        integer, intent(in):: lvl
        integer, intent(in):: befor_after
        
        if (befor_after == 0) then
            a_carbon_m(lvl) = self%a_carbon
        else
            self%a_carbon = a_carbon_m(lvl)
        end if
        
    end subroutine rewrite_algae
    
    subroutine do_grid(self, lvl, hice)
    !it makes grid, bottom layer is on lower edge of ice and equals 3 cm
    !other layers depth = ice_thickness  - 3 cm / number of layers - 1
    !calculates z and dz for current layer and current day
        
        class(ice_layer):: self
        integer,  intent(in):: lvl
        real(rk), intent(in):: hice
        real(rk)            :: non_b_layers, delta
        
        non_b_layers = number_of_layers - 1.

        !z_s is for bottom layer depth (3cm)
        delta = (hice - z_s) / non_b_layers
        if (lvl == number_of_layers) then
            self%z = hice
            z_m(lvl) = self%z
            dz_m(lvl) = z_s
        else
            self%z = lvl * delta
            z_m(lvl) = self%z
            dz_m(lvl) = delta
        end if
        
    end subroutine do_grid
    
    subroutine do_par(self, lvl, io, io_ice_var, snow_thick)
    !io in Watts, to calculate it in micromoles photons per m2*s =>
    !=> [w] = 4.6*[micromole photons]
    !calculates irradiance par_z
    !Grenfell and Maykutt 1977 indicate that the magnitude and shape
    !of the albedo curves depend strongly on the amount of liquid
    !water present in the upper part of the ice, so it fluctuates
    !throught year (true also for extinction coefficient)
        
        class(ice_layer):: self
        integer,  intent(in)    :: lvl
        real(rk), intent(in)    :: snow_thick
        real(rk), intent(in)    :: io
        real(rk), intent(out)   :: io_ice_var
        real(rk)                :: io_e             !io in micromoles
        real(rk)                :: par_alb, par_scat
        
        io_e = io / 4.6
        !if (snow_thick <= 0.005) then !albedo influence
        par_alb = io_e * (1. - alb_ice)
        !else
        !    par_alb = io_e * (1. - alb_snow) * exp(-k_snow * snow_thick)
        !end if
        
        par_scat = par_alb * io_ice_par !after scattered surface of ice
        self%par_z = par_scat * exp(-k_ice * self%z)

        if (lvl == number_of_layers) then
            io_ice_var = 4.6 * self%par_z
        else
            io_ice_var = -1.
        end if
        
    end subroutine do_par

    subroutine do_congelation_algae(self, lvl, ice_growth_temp_in)
    !evaluate transport of the algae
    !caused by freezing/melting

        class(ice_layer):: self
        integer, intent(in):: lvl
        real(rk), intent(in):: ice_growth_temp_in

        real(rk):: delta1, delta2, cache
        real(rk):: ice_growth_temp
        real(rk):: motility

        ice_growth_temp = ice_growth_temp_in

        !for initializing algae position after summer for example
        if ((lvl == (number_of_layers - 1)) .and. a_b == 0.) then
            a_b = self%z
        end if

        if (lvl /= number_of_layers) then
        !ice increasing
        !recalculate layer by layer
        !uses new layers
            if (a_b >= self%z) then
                return
            else
                if ((lvl == (number_of_layers - 1)) .and. &
                    (ice_growth_temp < 0.015)) then
                    motility = (1 - (ice_growth_temp / 0.015)) / 100.
                    a_b = min(self%z, a_b + motility)
                end if
                
                delta1 = self%z - a_b
                delta2 = min(1., delta1 / dz_m(lvl+1)) !part of the lower layer
                a_carbon_m(lvl) = a_carbon_m(lvl) + delta2 * a_carbon_m(lvl+1)
                a_carbon_m(lvl+1) = a_carbon_m(lvl+1) -&
                    delta2 * a_carbon_m(lvl+1)
            end if
        end if

    end subroutine do_congelation_algae

    subroutine do_melting_algae(self, lvl, ice_growth_temp_in, da_c)
    !evaluate transport of the algae
    !caused by freezing/melting

        class(ice_layer):: self
        integer, intent(in):: lvl
        real(rk), intent(in):: ice_growth_temp_in
        real(rk), intent(out):: da_c

        real(rk):: delta0, delta1, delta2, cache
        real(rk):: ice_growth_temp
        integer:: i, m, k

        ice_growth_temp = ice_growth_temp_in
        da_c = 0.

        if(ice_growth_temp < 0. .and. (trigger_melting .eqv. .false.)) then
        !recalculate all layers once, so we need trigger
        !uses old layers
            ice_growth_temp = abs(ice_growth_temp)
            delta0 = z_m(lvl) - ice_growth_temp
            if (delta0 <= a_b) then
                do k = lvl, 1, -1
                    da_c = da_c + a_carbon_m(k)
                    a_carbon_m(k) = 0.
                end do
                a_b = delta0 - dz_m(lvl)
            else if (delta0 > a_b .and. ice_growth_temp > dz_m(lvl)) then
                i = lvl
                do while (ice_growth_temp > dz_m(i))
                    da_c = da_c + a_carbon_m(i)
                    a_carbon_m(i) = 0.
                    ice_growth_temp = ice_growth_temp - dz_m(i)
                    i = i - 1
                end do
                delta1 = ice_growth_temp/dz_m(i)
                da_c = da_c + delta1*a_carbon_m(i)
                a_carbon_m(i) = &
                        a_carbon_m(i) - a_carbon_m(i)*delta1

                delta2 = dz_m(i) - ice_growth_temp
                cache = delta2 - dz_m(lvl)
                if (dz_m(lvl) >= delta2) then
                !write to bottom
                    m = 0
                    do k = i, 1, -1
                        a_carbon_m(lvl-m) = a_carbon_m(i)
                        a_carbon_m(i) = 0.
                        m = m + 1
                    end do
                else
                    m = 0
                    do k = i, 2, -1
                        if (i == k) then
                            a_carbon_m(lvl-m) = a_carbon_m(lvl-m)+&
                                a_carbon_m(k)*dz_m(lvl-m)/delta2
                            if (lvl-m-1 == k) then
                                a_carbon_m(lvl-m-1) = a_carbon_m(k)*cache/delta2
                            else
                                a_carbon_m(lvl-m-1) = a_carbon_m(lvl-m-1)+&
                                    a_carbon_m(k)*cache/delta2
                            end if
                            m = m + 1
                        else
                            a_carbon_m(lvl-m) = a_carbon_m(lvl-m)+&
                                a_carbon_m(k)*(dz_m(k)-cache)/dz_m(k)
                            if (lvl-m-1 == k) then
                                a_carbon_m(lvl-m-1) = a_carbon_m(k)*cache/dz_m(k)
                            else
                                a_carbon_m(lvl-m-1) = a_carbon_m(lvl-m-1)+&
                                    a_carbon_m(k)*cache/dz_m(k)
                            end if
                            m = m + 1
                        end if
                    end do
                end if
            else
                delta1 = ice_growth_temp/dz_m(lvl)
                da_c = da_c + delta1*a_carbon_m(lvl)
                a_carbon_m(lvl) = &
                    a_carbon_m(lvl) - a_carbon_m(lvl)*delta1
                do i = lvl, 2, -1
                    a_carbon_m(i) = &
                        a_carbon_m(i) + a_carbon_m(i-1)*&
                        ice_growth_temp/dz_m(i-1)
                    a_carbon_m(i-1) = &
                        a_carbon_m(i-1) - a_carbon_m(i-1)*&
                        ice_growth_temp/dz_m(i-1)
                end do
            end if
            if (a_b > (delta0 - dz_m(lvl))) a_b = delta0 - dz_m(lvl)
            ice_growth_temp = 0.
            trigger_melting = .true.
        else
            return
        end if

    end subroutine do_melting_algae
    
    subroutine do_bulk_temperature(self, air_temp, water_temp, ice_thickness)
    ![C]
    
        class(ice_layer):: self
        real(rk), intent(in):: air_temp, water_temp, ice_thickness
        
        self%bulk_temperature = air_temp + ((water_temp - air_temp) * self%z) / ice_thickness
        if (self%bulk_temperature > -0.2) self%bulk_temperature = -0.2
    
    end subroutine do_bulk_temperature
    
    subroutine do_bulk_salinity(self, lvl, ice_thickness, water_sal)
    !bulk salinity through ice [ppt]
    
        class(ice_layer):: self
        integer,  intent(in):: lvl
        real(rk), intent(in):: ice_thickness
        real(rk), intent(in):: water_sal
        real(rk):: z_p !ratio btwn the distance from the ice surface and ice thickness
        
        if (lvl == number_of_layers) then
            self%bulk_salinity = water_sal
        else
            z_p = self%z / ice_thickness
            self%bulk_salinity = 19.539 * (z_p**2) - 19.93 * z_p + 8.913
        end if
    
    end subroutine do_bulk_salinity
    
    subroutine do_brine_salinity(self, lvl)
    !brine salinity through ice in [ppt]
    
        class(ice_layer):: self
        integer,  intent(in):: lvl
        
        
        if (lvl == number_of_layers) then
            self%brine_salinity = self%bulk_salinity
        else
            if (self%bulk_temperature < 0 .and. self%bulk_temperature >=  -22.9) then
                self%brine_salinity = -3.9921 + (-22.700   * self%bulk_temperature) +&
                                                (-1.0015   * self%bulk_temperature**2) +&
                                                (-0.019956 * self%bulk_temperature**3)
            else if (self%bulk_temperature < -22.9 .and. self%bulk_temperature >=  -44) then
                self%brine_salinity = 206.24 + (-1.8907    * self%bulk_temperature) +&
                                               (-0.060868  * self%bulk_temperature**2) +&
                                               (-0.0010247 * self%bulk_temperature**3)
            else if (self%bulk_temperature < -44) then
                self%brine_salinity = -4442.1 + (-277.86  * self%bulk_temperature) +&
                                                (-5.501   * self%bulk_temperature**2) +&
                                                (-0.03669 * self%bulk_temperature**3)
            end if
        end if
    
    end subroutine do_brine_salinity
    
    subroutine do_brine_density(self)
    !brine density through ice [g*m-3]
    
        class(ice_layer):: self
        real(rk):: c ![g*cm-1*ppt-1]
        
        c = 8e-4
        self%brine_density = (1. + c*self%brine_salinity)*1e6
    
    end subroutine do_brine_density
    
    subroutine do_bulk_density(self)
    !bulk density through ice [g*m-3]
    
        class(ice_layer):: self
        real(rk):: dens_pure !density of pure ice
        real(rk):: bs
        
        dens_pure = 912000. ![g*m-3]
        bs = self%brine_salinity
        if (bs < 2.) bs = 2.! function has spikes below x=1
        self%bulk_density = dens_pure * self%brine_density * bs /&
        (self%brine_density * bs - self%bulk_salinity *&
        (self%brine_density - dens_pure))
    
    end subroutine do_bulk_density
    
    subroutine do_brine_relative_volume(self, lvl)
    !at the bottom layer it assumed equal to 0.5
    
        class(ice_layer):: self
        integer,  intent(in):: lvl
        real(rk):: bs
        
        if (lvl == number_of_layers) then
            self%brine_relative_volume = 0.5
        else
            bs = self%brine_salinity
            if (bs < 10.) bs = 10.!to fix above 1 values
            self%brine_relative_volume = (self%bulk_density * self%bulk_salinity) /&
                                         (self%brine_density  * bs)
        end if
    
    end subroutine do_brine_relative_volume
    
    subroutine do_nitrogen(self, lvl, nh4, no2, no3)
    !calculate fluxes of nitrogen in brine channels
    !per day
    
        class(ice_layer):: self
        integer,  intent(in):: lvl
        real(rk), intent(in):: nh4 !from upper water layer
        real(rk), intent(in):: no2 !from upper water layer
        real(rk), intent(in):: no3 !from upper water layer
        
        real(rk):: v_ammonium, v_nina, release_n

        v_ammonium = self%v_ammonium()
        v_nina = self%v_nina()
        release_n = self%release_n()
        
        !calculates nitrogen increment for bottom layer
        if (lvl == number_of_layers) then
            call self%do_bottom(nh4, self%d_nh4, self%a_nitrogen, &
                                v_ammonium, release_n, &
                                nh4_m)
            call self%do_bottom(no2, self%d_no2, self%a_nitrogen, &
                                (v_nina / 2.), release_n, &
                                no2_m)
            call self%do_bottom(no3, self%d_no3, self%a_nitrogen, &
                                (v_nina / 2.), release_n, &
                                no3_m)
        !and for rest of layers
        else
            call self%do_congelation(self%d_nh4, self%a_nitrogen, &
                            v_ammonium, release_n, lvl, &
                            nh4_m)
            call self%do_congelation(self%d_no2, self%a_nitrogen, &
                            (v_nina / 2.), release_n, lvl, &
                            no2_m)
            call self%do_congelation(self%d_no3, self%a_nitrogen, &
                            (v_nina / 2.), release_n, lvl, &
                            no3_m)
        end if
    
    end subroutine do_nitrogen
    
    subroutine do_phosphorus(self, lvl, po4)
    !calculate fluxes of phosphorus in brine channels
    !per day
    
        class(ice_layer):: self
        integer,  intent(in):: lvl
        real(rk), intent(in):: po4

        real(rk):: uptake_p, release_p

        uptake_p = self%uptake_p()
        release_p = self%release_p()
        
        !calculates d_po4 for bottom layer
        if (lvl == number_of_layers) then
            call self%do_bottom(po4, self%d_po4, self%a_phosphorus, &
                                (uptake_p * 24.), release_p, &
                                po4_m)
        !and for rest of layers
        else
            call self%do_congelation(self%d_po4, self%a_phosphorus, &
                           (uptake_p * 24.), release_p, lvl, &
                           po4_m)
        end if
    
    end subroutine do_phosphorus
    
    subroutine do_bottom(self, nut_in_water, d_nut, a_nut, uptake, release, nut_brine)
    !calculate bottom fluxes
    !per day
    
        class(ice_layer):: self
        real(rk), intent(in) :: nut_in_water, uptake, release
        real(rk), intent(in) :: a_nut !in algae nutrients
        real(rk), intent(out):: d_nut !increment of calculated brine nutrient
        real(rk), dimension(:), intent(in) :: nut_brine
        
        d_nut = self%brine_flux_z(number_of_layers, nut_brine, nut_in_water) +&
                self%brine_flux_s(number_of_layers, nut_brine, nut_in_water) +&
                self%diffusion_flux_s(number_of_layers, nut_brine, nut_in_water) +&
                self%congelation_flux_s(number_of_layers, nut_brine, nut_in_water) -&
                a_nut * (uptake - release)
    
    end subroutine do_bottom
    
    subroutine do_congelation(self, d_nut, a_nut, uptake, release, lvl, nut_brine)
    !brine concentrarion changes
    !per day
    
        class(ice_layer):: self
        integer, intent(in)  :: lvl
        real(rk), intent(in) :: uptake, release
        real(rk), intent(in) :: a_nut !in algae nutrients
        real(rk), intent(out):: d_nut
        real(rk), dimension(:), intent(in) :: nut_brine
        
        d_nut = self%brine_flux_z(lvl, nut_brine) -&
                a_nut * (uptake - release)
    
    end subroutine do_congelation
    
    function brine_flux_z(self, lvl, nut_brine, nut_in_water)
    !brine flux for all lvls - per day
    
        class(ice_layer):: self
        real(rk)                        :: brine_flux_z
        integer, intent(in)             :: lvl
        real(rk), intent(in), optional  :: nut_in_water !for bottom lvl
        real(rk), dimension(:), intent(in) :: nut_brine

        if (present(nut_in_water)) then
            brine_flux_z = self%brine_relative_volume *&
            gravity_drainage * z_conv *&
            self%second_derivative(nut_brine(lvl - 1), nut_brine(lvl),&
            nut_in_water, dz_m(lvl - 1), dz_m(lvl))
        else
            brine_flux_z = self%brine_relative_volume *&
            gravity_drainage * z_conv *&
            self%second_derivative(nut_brine(lvl - 1), nut_brine(lvl),&
            nut_brine(lvl + 1), dz_m(lvl - 1), dz_m(lvl))
        end if
    
        brine_flux_z = 86400. * brine_flux_z !convert per second to per day
    
    end function brine_flux_z
    
    function brine_flux_s(self, lvl, nut_brine, nut_in_water)
    !brine flux for bottom - per day
    
        class(ice_layer):: self
        integer, intent(in)             :: lvl
        real(rk)                        :: brine_flux_s
        real(rk), intent(in)            :: nut_in_water !for bottom layer
        real(rk), dimension(:), intent(in) :: nut_brine
        
        brine_flux_s = self%brine_relative_volume *&
        self%brine_release() * z_s *&
        self%second_derivative(nut_brine(lvl - 1), nut_brine(lvl),&
        nut_in_water, dz_m(lvl - 1), dz_m(lvl))

        brine_flux_s = 86400. * brine_flux_s !convert per second to per day
    
    end function brine_flux_s
    
    function diffusion_flux_s(self, lvl, nut_brine, nut_in_water)
    !brine diffusion flux for bottom - per day
    
        class(ice_layer):: self
        integer, intent(in)             :: lvl
        real(rk)                        :: diffusion_flux_s
        real(rk), intent(in)            :: nut_in_water !for bottom layer
        real(rk), dimension(:), intent(in) :: nut_brine
        
        diffusion_flux_s = self%brine_relative_volume * k_wi *&
        self%second_derivative(nut_brine(lvl - 1), nut_brine(lvl),&
        nut_in_water, dz_m(lvl - 1), dz_m(lvl))

        diffusion_flux_s = 86400. * diffusion_flux_s !convert per second to per day
    
    end function diffusion_flux_s
    
    function congelation_flux_s(self, lvl, nut_brine, nut_in_water)
    !brine diffusion flux for bottom - per day
    
        class(ice_layer):: self
        integer, intent(in)             :: lvl
        real(rk)                        :: congelation_flux_s
        real(rk), intent(in)            :: nut_in_water !for bottom layer
        real(rk), dimension(:), intent(in) :: nut_brine
        
        congelation_flux_s = (ice_growth * (nut_in_water - nut_brine(lvl)) * nut_in_water) /&
                             ((z_s + ice_growth) * self%brine_relative_volume * 86400.)
    
    end function congelation_flux_s
    
    function second_derivative(self, pre_layer, layer, next_layer, pre_h, next_h)
    !second_derivative
    
        class(ice_layer):: self
        real(rk), intent(in)            :: pre_layer, layer, next_layer, pre_h, next_h
        real(rk)                        :: second_derivative
        real(rk)                        :: h
        
        h = (pre_h + next_h) / 2.
        second_derivative = (pre_layer - 2. * layer + next_layer) / h**2
    
    end function second_derivative
    
    function brine_release(self)
    !brine_release [m s-1]
    
        class(ice_layer):: self
        real(rk)         :: brine_release
        
        brine_release = (9.667e-9 + 4.49e-6 * ice_growth - 1.39e-7 * ice_growth**2) / 100.
    
    end function brine_release
    
    subroutine do_a_carbon(self, lvl)
    !brine concentrations changes of ice algae carbon
    !mg C per hour
    
        class(ice_layer):: self
        integer, intent(in):: lvl
        real(rk):: gpp = 0., resp = 0., &
            exud = 0., mort = 0., melt = 0.

        gpp = self%gpp()
        resp = self%resp()
        exud = self%exud()
        mort = self%mort()

        if (lvl == number_of_layers) then
            melt = self%melt()
            self%d_a_carbon = self%a_carbon * (gpp - resp -&
                      exud - mort - melt) + recruit
        else
            self%d_a_carbon = self%a_carbon * (gpp - resp -&
                      exud - mort)
        end if
    
    end subroutine do_a_carbon
    
    function gpp(self)
    !gross photosynthetic production - per hour
    !dimensionless
    
        class(ice_layer):: self
        real(rk):: gpp
        real(rk):: f_par = 0., f_s = 0., &
            f_t = 0., f_nut = 0.
        
        f_par = self%f_par()
        f_s = self%f_s()
        f_t = self%f_t()
        f_nut = self%f_nut()
        
        gpp = p_b * f_par * f_s * f_t *&
              f_nut / chl_to_carbon
        self%last_gpp = gpp
    
    end function gpp
    
    function f_par(self)
    !light limitation, dimensionless

        class(ice_layer):: self
        real(rk):: f_par
        
        f_par = (1. - exp(-1. * alpha * self%par_z / p_b)) &
                   * exp(-1. * betta * self%par_z / p_b)
        
    end function f_par
    
    function f_s(self)
    !salinity limitation, dimensionless

        class(ice_layer):: self
        real(rk):: f_s
        
        f_s = 1.1e-2 + 3.012e-2  * self%brine_salinity    +&
                       1.0342e-3 * self%brine_salinity**2 +&
                       4.6033e-5 * self%brine_salinity**3 +&
                       4.926e-7  * self%brine_salinity**4 +&
                       1.659e-9  * self%brine_salinity**5
    
    end function f_s
    
    function f_t(self)
    !temperature limitation, dimensionless

        class(ice_layer):: self
        real(rk):: f_t

        f_t = exp(temp_aug_rate * (self%brine_temperature - t_0))
        self%last_f_t = f_t
    
    end function f_t
    
    function f_nut(self)
    !nutrient limitation, dimensionless

        class(ice_layer):: self
        real(rk):: f_nut
        real(rk):: n_cell
        real(rk):: p_cell

        n_cell = self%n_cell()
        p_cell = self%p_cell()
        
        if (self%a_carbon == 0.) then
            f_nut = 0.
        else
            f_nut = min(n_cell / (kn_cell + n_cell), p_cell / (kp_cell + p_cell))
        end if
    
    end function f_nut
    
    function n_cell(self)
    !nitrogen cell quota - [mg N mg C-1]
    
        class(ice_layer):: self
        real(rk):: n_cell
        
        n_cell = self%a_nitrogen / self%a_carbon
        if(isnan(n_cell)) n_cell = 0.
    
    end function n_cell
    
    function p_cell(self)
    !phosphorus cell quota - [mg P mg C-1]
    
        class(ice_layer):: self
        real(rk):: p_cell
        
        p_cell = self%a_phosphorus / self%a_carbon
        if(isnan(p_cell)) p_cell = 0.
    
    end function p_cell
    
    function resp(self)
    !respiration - per hour
    
        class(ice_layer):: self
        real(rk):: resp
        real(rk):: resp_day, resp_night
        
        resp_day = (r0 + r_dark * dl_ratio * self%last_gpp) *&
            (carbon_to_oxy / chl_to_carbon) *&
            16. * day_length * self%last_f_t
        
        resp_night = (r0 + r_dark * self%last_gpp) *&
            (carbon_to_oxy / chl_to_carbon) *&
            16. * (24. - day_length) * self%last_f_t
        
        !to make resp coeffic in per hour units
        resp = (resp_day + resp_night) / 24.
    
    end function resp
    
    function exud(self)
    !ice algae exudation rate - per hour
    
        class(ice_layer):: self
        real(rk):: exud
        
        exud = exud_rate * self%last_gpp
    
    end function exud
    
    function mort(self)
    !mortality loss - per hour

        class(ice_layer):: self
        real(rk):: mort

        mort = phy_mort * self%last_f_t! or nutrient concentration should be here?

        self%last_mort = mort
    
    end function mort
    
    function melt(self)
    !melting rate - per hour

        class(ice_layer):: self
        real(rk):: melt
        
        melt =  -1. * ice_growth / z_s / 24. !24 is to transform in per hour units
        if (melt < 0) melt = 0.
    
    end function melt
    
    subroutine do_a_nitrogen(self, lvl)
    !brine concentrations changes of ice algae nitrogen
    !mg N per hour
    
        class(ice_layer):: self
        integer, intent(in):: lvl

        if (lvl == number_of_layers) then
            self%d_a_nitrogen = self%a_nitrogen * (self%uptake_n() - self%release_n() -&
                      self%last_mort - self%melt()) + recruit
        else
            self%d_a_nitrogen = self%a_nitrogen * (self%uptake_n() - self%release_n() -&
                      self%last_mort)
        end if
    
    end subroutine do_a_nitrogen
    
    function uptake_n(self)
    !a lot of problems could be here - uptake of N by algae - per hour

        class(ice_layer):: self
        real(rk):: uptake_n
        
        uptake_n = (self%v_ammonium() + self%v_nina()) / 24. !to per hour transform
    
    end function uptake_n
    
    function v_ammonium(self)
    !ammonium uptake rate - per day
    
        class(ice_layer):: self
        real(rk):: v_ammonium

        real(rk):: foo
        
        foo = self%a_nitrogen / self%a_phosphorus
        if (self%a_nitrogen > n_min .and. self%a_nitrogen < n_max &
            .and. foo < max_n_p) then
            v_ammonium = (v_max_n * self%nh4) / (k_ammonium + self%nh4) *&
                (1 - (self%n_cell() / n_max))
        else
            v_ammonium = 0.
        end if

        self%last_v_ammonium = v_ammonium
    
    end function v_ammonium
    
    function v_nina(self)
    !nitrite and nitrate uptake rate - per day

        class(ice_layer):: self
        real(rk):: v_nina

        real(rk):: foo
        
        foo = self%a_nitrogen / self%a_phosphorus
        if (self%a_nitrogen > n_min .and. self%a_nitrogen < n_max &
            .and. foo < max_n_p) then
            v_nina = max(0., (v_max_n - self%last_v_ammonium)) *&
                (self%no2 + self%no3) / (k_nina + self%no2 + self%no3) *&
                (1 - (self%n_cell() / n_max))
        else
            v_nina = 0.
        end if
    
    end function v_nina
    
    function release_n(self)
    !release - release everything more then n_max - mg N
    
        class(ice_layer):: self
        real(rk):: release_n
        
        if (self%n_cell() > n_max .and. self%a_nitrogen > 0.) then
            release_n = (self%a_nitrogen - n_max * self%a_carbon) / self%a_nitrogen
        else
            release_n = 0.
        end if
    
    end function release_n
    
    subroutine do_a_phosphorus(self, lvl)
    !brine concentrations changes of ice algae po4
    !mg N per hour
    
        class(ice_layer):: self
        integer, intent(in):: lvl

        if (lvl == number_of_layers) then
            self%d_a_phosphorus = self%a_phosphorus * (self%uptake_p() - self%release_p() -&
                      self%last_mort - self%melt()) + recruit
        else
            self%d_a_phosphorus = self%a_phosphorus * (self%uptake_p() - self%release_p() -&
                      self%last_mort)
        end if
    
    end subroutine do_a_phosphorus
    
    function uptake_p(self)
    !a lot of problems could be here - uptake of N by algae - per hour

        class(ice_layer):: self
        real(rk):: uptake_p
        real(rk):: foo
        
        foo = self%a_nitrogen / self%a_phosphorus
        if (self%a_phosphorus > phos_min .and. self%a_phosphorus < phos_max &
            .and. foo > min_n_p) then
            uptake_p = (v_max_p * self%po4) / (k_p + self%po4) *&
                (1. - (self%p_cell() / phos_max))
        else
            uptake_p = 0.
        end if
        uptake_p = uptake_p / 24. !to per hour transform
        
    end function uptake_p
    
    function release_p(self)
    !release - release everything more then n_max - mg N
    
        class(ice_layer):: self
        real(rk):: release_p
        
        if (self%p_cell() > phos_max .and. self%a_phosphorus > 0.) then
            release_p = (self%a_phosphorus - phos_max * self%a_carbon) / self%a_phosphorus
        else
            release_p = 0.
        end if
    
    end function release_p
    
end module ice_algae_lib
    
!-----------------------------------------------------------------------
! Copyright under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------