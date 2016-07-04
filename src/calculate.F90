module calculate

    use fabm_types, only: rk

    implicit none
    private

    public calculate_phys, calculate_sed

contains

    subroutine calculate_phys(cc, lev_max, par_max, use_bound_up, &
    use_bound_low, bound_up, bound_low, surf_flux, &
    boundary_bbl_sediments, kz2, julianday, kz_bio, i_O2, &
    dz, freq_az, dcc, fick)


        integer, intent(in)                                 :: lev_max
        integer, intent(in)                                 :: par_max
        real(rk), dimension(lev_max, par_max), intent(inout):: dcc, fick
        real(rk), dimension(:, :), intent(inout) :: cc
        real(rk), dimension(:), intent(in)       :: surf_flux
        real(rk), dimension(:), intent(in)       :: bound_up
        real(rk), dimension(:), intent(in)       :: bound_low
        logical, dimension(:), intent(in)        :: use_bound_low
        logical, dimension(:), intent(in)        :: use_bound_up
        real(rk), intent(in)                     :: kz2(:, :)
        real(rk), intent(in)                     :: kz_bio(:)
        real(rk), intent(in)                     :: dz(:)
        integer, intent(in)                      :: boundary_bbl_sediments
        integer, intent(in)                      :: julianday
        integer, intent(in)                      :: i_O2
        integer, intent(in)                      :: freq_az

        integer                                  :: k, ip

        do k = 2, (lev_max - 1)
        !boundary conditions
            !upper boundary
            if (k .eq. 2) then
                do ip = 1, par_max
                    if (use_bound_up(ip)) then
                        cc(k - 1, ip) = bound_up(ip)
                    else
                        cc(k - 1, ip) = cc(k, ip)
                    end if
                enddo
            end if
            !low Boundary
            if (k .eq. (lev_max - 1)) then
                do ip = 1, par_max
                    if (use_bound_low(ip)) then
                        cc(k + 1, ip) = bound_low(ip)
                    else
                        cc(k + 1, ip) = cc(k, ip)
                    end if
                end do
            end if
            !turbulence
            !fluxes of parameters, bioturbation (kz_bio) depend on
            !Zoo or O2 availability at the surface
            fick(k, :) = (kz2(k, julianday) + kz_bio(k) &
                * cc(boundary_bbl_sediments - 1, i_O2) / &
                (cc(boundary_bbl_sediments - 1, i_O2) + 5.)) &
                * (cc(k + 1, :) - cc(k, :)) / dz(k)
            fick(k - 1, :) = (kz2(k - 1, julianday) + kz_bio(k - 1) &
                * cc(boundary_bbl_sediments - 1, i_O2) / &
                (cc(boundary_bbl_sediments - 1, i_O2) + 5.)) &
                * (cc(k, :) - cc(k - 1, :)) / dz(k - 1)
            dcc(k, :) = (fick(k, :) - fick(k - 1, :)) / &
                ((dz(k) + dz(k - 1)) / 2.)
            !fluxes of O2 and DIC from the air:
            if (k .eq. 2) then
                dcc(k, :) = dcc(k, :) + surf_flux(:) / dz(k)
            end if
        end do

    end subroutine calculate_phys

    subroutine calculate_sed(par_max, lev_max, wbio, cc, dz, &
    boundary_bbl_sediments, dcc)

        real(rk), dimension(:, :), intent(inout) :: cc
        real(rk), dimension(:, :), intent(in)    :: wbio
        real(rk), intent(in)                     :: dz(:)
        integer, intent(in)                      :: par_max
        integer, intent(in)                      :: lev_max
        integer, intent(in)                      :: boundary_bbl_sediments
        real(rk), dimension(lev_max, par_max), intent(inout):: dcc

        real(rk)                                 :: w_u(par_max)
        real(rk)                                 :: w_d(par_max)
        real(rk), parameter                      :: w_buruing = 0.
        integer                                  :: k

        do  k = 2, lev_max - 1
        !boundary conditions

            !from "up": we check either we are in water or sediment
            if (k < boundary_bbl_sediments) then
                !in FABM: negative w means sinking!
                w_u(:) = -wbio(k, :)
            else
                w_u(:) = w_buruing
            end if

            !to "down"
            if ((k + 1) < boundary_bbl_sediments ) then
                w_d(:) = -wbio(k, :)
            else
                w_d(:) = w_buruing
            end if

            !upper boundary
            if (k < 3) then
                w_u = 0.0
            end if

            !multiplication by 0.5 was deleted
            dcc(k, :) = (w_u * cc(k - 1, :) / ((dz(k - 1) + dz(k)) /2.) - &
                w_d * cc(k, :) / ((dz(k - 1) + dz(k)) /2.))
        end do

    end subroutine calculate_sed

end module calculate

!-----------------------------------------------------------------------
! Copyright under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
