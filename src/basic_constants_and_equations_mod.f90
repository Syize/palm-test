!> @file basic_constants_and_equations_mod.f90
!--------------------------------------------------------------------------------------------------!
! This file is part of the PALM model system.
! PALM is free software: you can redistribute it and/or modify it under the terms of the GNU General
! Public License as published by the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! PALM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
! implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
! Public License for more details.
! You should have received a copy of the GNU General Public License along with PALM. If not, see
! <http://www.gnu.org/licenses/>.
! Copyright 1997-2021 Leibniz Universitaet Hannover
! Copyright 2022-2022 pecanode GmbH
!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This module contains all basic (physical) constants and functions for the calculation of
!> diagnostic quantities.
!--------------------------------------------------------------------------------------------------!
module basic_constants_and_equations_mod

    use kinds

    implicit none

    real(wp), parameter :: c_p = 1005.0_wp !< heat capacity of dry air (J kg-1 K-1)
    real(wp), parameter :: c_w = 4185.0_wp !< heat capacity of water at 0°C (J kg-1 K-1)
    real(wp), parameter :: degc_to_k = 273.15_wp !< temperature (in K) of 0 deg C (K)
    real(wp), parameter :: g = 9.81_wp !< gravitational acceleration (m s-2)
    real(wp), parameter :: kappa = 0.4_wp !< von Karman constant
    real(wp), parameter :: k_boltzmann = 1.380649e-23_wp !< Boltzmann's constant (J K-1)
    real(wp), parameter :: l_m = 0.33e+06_wp !< latent heat of water melting (J kg-1)
    real(wp), parameter :: l_v = 2.5e+06_wp !< latent heat of water vaporization (J kg-1)
    real(wp), parameter :: l_s = l_m + l_v !< latent heat of water sublimation (J kg-1)
    real(wp), parameter :: molecular_weight_of_nacl = 0.05844_wp !< mol. m. NaCl (kg mol-1)
    real(wp), parameter :: molecular_weight_of_c3h4o4 = 0.10406_wp !< mol. m. malonic acid (kg mol-1)
    real(wp), parameter :: molecular_weight_of_nh4no3 = 0.08004_wp !< mol. m. ammonium sulfate (kg mol-1)
    real(wp), parameter :: molecular_weight_of_water = 0.01801528_wp !< mol. m. H2O (kg mol-1)
    real(wp), parameter :: pi = 3.141592654_wp !< PI
    !$ACC DECLARE COPYIN(pi)
    real(wp), parameter :: rgas_univ = 8.31446261815324_wp !< universal gas constant (J K-1 mol-1)
    real(wp), parameter :: rho_i = 916.7_wp !> density of pure ice (kg m-3)
    real(wp), parameter :: rho_l = 1.0e3_wp !< density of water (kg m-3)
    real(wp), parameter :: rho_nacl = 2165.0_wp !< density of NaCl (kg m-3)
    real(wp), parameter :: rho_c3h4o4 = 1600.0_wp !< density of malonic acid (kg m-3)
    real(wp), parameter :: rho_nh4no3 = 1720.0_wp !< density of ammonium sulfate (kg m-3)
    real(wp), parameter :: r_d = 287.0_wp !< sp. gas const. dry air (J kg-1 K-1)
    real(wp), parameter :: r_v = 461.51_wp !< sp. gas const. water vapor (J kg-1 K-1)
    real(wp), parameter :: sigma_sb = 5.67037e-08_wp !< Stefan-Boltzmann constant
    real(wp), parameter :: solar_constant = 1368.0_wp !< solar constant at top of atmosphere
    real(wp), parameter :: vanthoff_nacl = 2.0_wp !< van't Hoff factor for NaCl
    real(wp), parameter :: vanthoff_c3h4o4 = 1.37_wp !< van't Hoff factor for malonic acid
    real(wp), parameter :: vanthoff_nh4no3 = 2.31_wp !< van't Hoff factor for ammonium sulfate

    real(wp), parameter :: p_0 = 100000.0_wp !< standard pressure reference state

    real(wp), parameter :: cp_d_rd = c_p / r_d !< precomputed c_p / r_d
    real(wp), parameter :: degrees_to_radiants = pi / 180.0_wp !< conversion factor between degrees and radiants
    real(wp), parameter :: g_d_cp = g / c_p !< precomputed g / c_p
    real(wp), parameter :: lv_d_cp = l_v / c_p !< precomputed l_v / c_p
    real(wp), parameter :: ls_d_cp = l_s / c_p !< precomputed l_s / c_p
    real(wp), parameter :: lv_d_rd = l_v / r_d !< precomputed l_v / r_d
    real(wp), parameter :: radiants_to_degrees = 180.0_wp / pi !< conversion factor between radiants and degrees
    real(wp), parameter :: rd_d_rv = r_d / r_v !< precomputed r_d / r_v
    real(wp), parameter :: rd_d_cp = r_d / c_p !< precomputed r_d / c_p

    real(wp) :: molecular_weight_of_solute = molecular_weight_of_nacl !< mol. m. NaCl (kg mol-1)
    real(wp) :: rho_s = rho_nacl !< density of NaCl (kg m-3)
    real(wp) :: vanthoff = vanthoff_nacl !< van't Hoff factor for NaCl

    save

    private magnus_0d, &
        magnus_1d, &
        magnus_tl_0d, &
        magnus_tl_1d, &
        magnus_0d_ice, &
        magnus_1d_ice, &
        ideal_gas_law_rho_0d, &
        ideal_gas_law_rho_1d, &
        ideal_gas_law_rho_pt_0d, &
        ideal_gas_law_rho_pt_1d, &
        exner_function_0d, &
        exner_function_1d, &
        exner_function_invers_0d, &
        exner_function_invers_1d, &
        barometric_formula_0d, &
        barometric_formula_1d, &
        get_relative_humidity_equilibrium_0d, &
        get_relative_humidity_equilibrium_1d, &
        get_relative_humidity_supersaturated_0d, &
        get_relative_humidity_supersaturated_1d

    interface convert_utm_to_geographic
        module procedure convert_utm_to_geographic
    end interface convert_utm_to_geographic

    interface magnus
        module procedure magnus_0d
        module procedure magnus_1d
    end interface magnus

    interface magnus_tl
        module procedure magnus_tl_0d
        module procedure magnus_tl_1d
    end interface magnus_tl

    interface magnus_ice
        module procedure magnus_0d_ice
        module procedure magnus_1d_ice
    end interface magnus_ice

    interface ideal_gas_law_rho
        module procedure ideal_gas_law_rho_0d
        module procedure ideal_gas_law_rho_1d
    end interface ideal_gas_law_rho

    interface ideal_gas_law_rho_pt
        module procedure ideal_gas_law_rho_pt_0d
        module procedure ideal_gas_law_rho_pt_1d
    end interface ideal_gas_law_rho_pt

    interface exner_function
        module procedure exner_function_0d
        module procedure exner_function_1d
    end interface exner_function

    interface exner_function_invers
        module procedure exner_function_invers_0d
        module procedure exner_function_invers_1d
    end interface exner_function_invers

    interface barometric_formula
        module procedure barometric_formula_0d
        module procedure barometric_formula_1d
    end interface barometric_formula

    interface get_relative_humidity
        module procedure get_relative_humidity_equilibrium_0d
        module procedure get_relative_humidity_equilibrium_1d
        module procedure get_relative_humidity_supersaturated_0d
        module procedure get_relative_humidity_supersaturated_1d
    end interface get_relative_humidity

    !-- Public function and routines
    public convert_utm_to_geographic

contains

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Convert UTM coordinates into geographic latitude and longitude. Conversion is based on the work
    !> of Krüger (1912) DOI: 10.2312/GFZ.b103-krueger28 and Karney (2013) DOI: 10.1007/s00190-012-0578-z
    !> Based on a JavaScript of the geodesy function library written by chrisveness
    !> https://github.com/chrisveness/geodesy
    !--------------------------------------------------------------------------------------------------!
    subroutine convert_utm_to_geographic(crs, eutm, nutm, lon, lat)

        integer(iwp) :: j !< loop index

        real(wp), intent(in) :: eutm !< easting (UTM)
        real(wp), intent(out) :: lat !< geographic latitude in degree
        real(wp), intent(out) :: lon !< geographic longitude in degree
        real(wp), intent(in) :: nutm !< northing (UTM)

        real(wp) :: a !< 2*pi*a is the circumference of a meridian
        real(wp) :: cos_eta_s !< cos(eta_s)
        real(wp) :: delta_i !<
        real(wp) :: delta_tau_i !<
        real(wp) :: e !< eccentricity
        real(wp) :: eta !<
        real(wp) :: eta_s !<
        real(wp) :: n !< 3rd flattening
        real(wp) :: n2 !< n^2
        real(wp) :: n3 !< n^3
        real(wp) :: n4 !< n^4
        real(wp) :: n5 !< n^5
        real(wp) :: n6 !< n^6
        real(wp) :: nu !<
        real(wp) :: nu_s !<
        real(wp) :: sin_eta_s !< sin(eta_s)
        real(wp) :: sinh_nu_s !< sinush(nu_s)
        real(wp) :: tau_i !<
        real(wp) :: tau_i_s !<
        real(wp) :: tau_s !<
        real(wp) :: x !< adjusted easting
        real(wp) :: y !< adjusted northing

        real(wp), dimension(6) :: beta !< 6th order Krüger expressions

        real(wp), dimension(8), intent(in) :: crs !< coordinate reference system, consists of
        !< (/semi_major_axis,
        !<   inverse_flattening,
        !<   longitude_of_prime_meridian,
        !<   longitude_of_central_meridian,
        !<   scale_factor_at_central_meridian,
        !<   latitude_of_projection_origin,
        !<   false_easting,
        !<   false_northing /)

        x = eutm - crs(7) ! remove false easting
        y = nutm - crs(8) ! remove false northing
        !-- From Karney 2011 Eq 15-22, 36:
        e = sqrt(1.0_wp / crs(2) * (2.0_wp - 1.0_wp / crs(2)))
        n = 1.0_wp / crs(2) / (2.0_wp - 1.0_wp / crs(2))
        n2 = n * n
        n3 = n * n2
        n4 = n * n3
        n5 = n * n4
        n6 = n * n5

        a = crs(1) / (1.0_wp + n) * (1.0_wp + 0.25_wp * n2 + 0.015625_wp * n4 + 3.90625e-3_wp * n6)

        nu = x / (crs(5) * a)
        eta = y / (crs(5) * a)

        !-- According to Krüger (1912), eq. 26*
        beta(1) = 0.5_wp * n &
                  - 2.0_wp / 3.0_wp * n2 &
                  + 37.0_wp / 96.0_wp * n3 &
                  - 1.0_wp / 360.0_wp * n4 &
                  - 81.0_wp / 512.0_wp * n5 &
                  + 96199.0_wp / 604800.0_wp * n6

        beta(2) = 1.0_wp / 48.0_wp * n2 &
                  + 1.0_wp / 15.0_wp * n3 &
                  - 437.0_wp / 1440.0_wp * n4 &
                  + 46.0_wp / 105.0_wp * n5 &
                  - 1118711.0_wp / 3870720.0_wp * n6

        beta(3) = 17.0_wp / 480.0_wp * n3 &
                  - 37.0_wp / 840.0_wp * n4 &
                  - 209.0_wp / 4480.0_wp * n5 &
                  + 5569.0_wp / 90720.0_wp * n6

        beta(4) = 4397.0_wp / 161280.0_wp * n4 &
                  - 11.0_wp / 504.0_wp * n5 &
                  - 830251.0_wp / 7257600.0_wp * n6

        beta(5) = 4583.0_wp / 161280.0_wp * n5 &
                  - 108847.0_wp / 3991680.0_wp * n6

        beta(6) = 20648693.0_wp / 638668800.0_wp * n6

        eta_s = eta
        nu_s = nu
        do j = 1, 6
            eta_s = eta_s - beta(j) * sin(2.0_wp * j * eta) * cosh(2.0_wp * j * nu)
            nu_s = nu_s - beta(j) * cos(2.0_wp * j * eta) * sinh(2.0_wp * j * nu)
        end do

        sinh_nu_s = sinh(nu_s)
        sin_eta_s = sin(eta_s)
        cos_eta_s = cos(eta_s)

        tau_s = sin_eta_s / sqrt(sinh_nu_s**2 + cos_eta_s**2)

        tau_i = tau_s
        delta_tau_i = 1.0_wp

        do while (abs(delta_tau_i) > 1.0e-12_wp)

            delta_i = sinh(e * ATANH(e * tau_i / sqrt(1.0_wp + tau_i**2)))

            tau_i_s = tau_i * sqrt(1.0_wp + delta_i**2) - delta_i * sqrt(1.0_wp + tau_i**2)

            delta_tau_i = (tau_s - tau_i_s) / sqrt(1.0_wp + tau_i_s**2) &
                          * (1.0_wp + (1.0_wp - e**2) * tau_i**2) &
                          / ((1.0_wp - e**2) * sqrt(1.0_wp + tau_i**2))

            tau_i = tau_i + delta_tau_i

        end do

        lat = atan(tau_i) / pi * 180.0_wp
        lon = atan2(sinh_nu_s, cos_eta_s) / pi * 180.0_wp + crs(4)

    end subroutine convert_utm_to_geographic

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> This function computes the magnus formula (Press et al., 1992).
    !> The magnus formula is needed to calculate the saturation vapor pressure.
    !--------------------------------------------------------------------------------------------------!
    function magnus_0d(t)
        !$ACC ROUTINE SEQ

        implicit none

        real(wp), intent(IN) :: t !< temperature (K)

        real(wp) :: magnus_0d

        !-- Saturation vapor pressure for a specific temperature:
        magnus_0d = 611.2_wp * exp(17.62_wp * (t - degc_to_k) / (t - 29.65_wp))

    end function magnus_0d

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> This function computes the magnus formula (Press et al., 1992).
    !> The magnus formula is needed to calculate the saturation vapor pressure.
    !--------------------------------------------------------------------------------------------------!
    function magnus_1d(t)

        implicit none

        real(wp), intent(IN), dimension(:) :: t !< temperature (K)

        real(wp), dimension(size(t)) :: magnus_1d

        !-- Saturation vapor pressure for a specific temperature:
        magnus_1d = 611.2_wp * exp(17.62_wp * (t - degc_to_k) / (t - 29.65_wp))

    end function magnus_1d

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> This function computes the magnus formula (Press et al., 1992) using the (ice-) liquid water
    !> potential temperature.
    !> The magnus formula is needed to calculate the saturation vapor pressure over a plane liquid water
    !> surface.
    !--------------------------------------------------------------------------------------------------!
    function magnus_tl_0d(t_l)

        implicit none

        real(wp), intent(IN) :: t_l !< liquid water temperature (K)

        real(wp) :: magnus_tl_0d

        !-- Saturation vapor pressure for a specific temperature:
        magnus_tl_0d = 610.78_wp * exp(17.269_wp * (t_l - 273.16_wp) / (t_l - 35.86_wp))

    end function magnus_tl_0d

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> This function computes the magnus formula (Press et al., 1992) using the (ice-) liquid water
    !> potential temperature.
    !> The magnus formula is needed to calculate the saturation vapor pressure over a plane liquid water
    !> surface.
    !--------------------------------------------------------------------------------------------------!
    function magnus_tl_1d(t_l)

        implicit none

        real(wp), intent(IN), dimension(:) :: t_l !< liquid water temperature (K)

        real(wp), dimension(size(t_l)) :: magnus_tl_1d
        !-- Saturation vapor pressure for a specific temperature:
        magnus_tl_1d = 610.78_wp * exp(17.269_wp * (t_l - 273.16_wp) / (t_l - 35.86_wp))

    end function magnus_tl_1d

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> This function computes the magnus formula (Press et al., 1992).
    !> The magnus formula is needed to calculate the saturation vapor pressure over a plane ice surface.
    !--------------------------------------------------------------------------------------------------!
    function magnus_0d_ice(t)

        implicit none

        real(wp), intent(IN) :: t !< temperature (K)

        real(wp) :: magnus_0d_ice

        !-- Saturation vapor pressure for a specific temperature:
        !magnus_0d_ice =  611.2_wp * EXP( 22.46_wp * ( t - degc_to_k ) / ( t - 0.53_wp  ) )
        magnus_0d_ice = 610.78_wp * exp(21.875_wp * (t - degc_to_k) / (t - 7.66_wp))

    end function magnus_0d_ice

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> This function computes the magnus formula (Press et al., 1992).
    !> The magnus formula is needed to calculate the saturation vapor pressure over a plane ice surface.
    !--------------------------------------------------------------------------------------------------!
    function magnus_1d_ice(t)

        implicit none

        real(wp), intent(IN), dimension(:) :: t !< temperature (K)

        real(wp), dimension(size(t)) :: magnus_1d_ice

        !-- Saturation vapor pressure for a specific temperature:
        !magnus_1d_ice =  611.2_wp * EXP( 22.46_wp * ( t - degc_to_k ) / ( t - 0.53_wp  ) )
        magnus_1d_ice = 610.78_wp * exp(21.875_wp * (t - degc_to_k) / (t - 7.66_wp))

    end function magnus_1d_ice

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Compute the ideal gas law for scalar arguments.
    !--------------------------------------------------------------------------------------------------!
    function ideal_gas_law_rho_0d(p, t)

        implicit none

        real(wp), intent(IN) :: p !< pressure (Pa)
        real(wp), intent(IN) :: t !< temperature (K)

        real(wp) :: ideal_gas_law_rho_0d

        !-- Compute density according to ideal gas law:
        ideal_gas_law_rho_0d = p / (r_d * t)

    end function ideal_gas_law_rho_0d

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Compute the ideal gas law for 1-D array arguments.
    !--------------------------------------------------------------------------------------------------!
    function ideal_gas_law_rho_1d(p, t)

        implicit none

        real(wp), intent(IN), dimension(:) :: p !< pressure (Pa)
        real(wp), intent(IN), dimension(:) :: t !< temperature (K)

        real(wp), dimension(size(p)) :: ideal_gas_law_rho_1d

        !-- Compute density according to ideal gas law:
        ideal_gas_law_rho_1d = p / (r_d * t)

    end function ideal_gas_law_rho_1d

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Compute the ideal gas law for scalar arguments.
    !--------------------------------------------------------------------------------------------------!
    function ideal_gas_law_rho_pt_0d(p, t)

        implicit none

        real(wp), intent(IN) :: p !< pressure (Pa)
        real(wp), intent(IN) :: t !< temperature (K)

        real(wp) :: ideal_gas_law_rho_pt_0d

        !-- Compute density according to ideal gas law:
        ideal_gas_law_rho_pt_0d = p / (r_d * exner_function(p) * t)

    end function ideal_gas_law_rho_pt_0d

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Compute the ideal gas law for 1-D array arguments.
    !--------------------------------------------------------------------------------------------------!
    function ideal_gas_law_rho_pt_1d(p, t)

        implicit none

        real(wp), intent(IN), dimension(:) :: p !< pressure (Pa)
        real(wp), intent(IN), dimension(:) :: t !< temperature (K)

        real(wp), dimension(size(p)) :: ideal_gas_law_rho_pt_1d

        !-- Compute density according to ideal gas law:
        ideal_gas_law_rho_pt_1d = p / (r_d * exner_function(p) * t)

    end function ideal_gas_law_rho_pt_1d

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Compute the exner function for scalar arguments.
    !--------------------------------------------------------------------------------------------------!
    function exner_function_0d(p)

        implicit none

        real(wp), intent(IN) :: p !< pressure (Pa)

        real(wp) :: exner_function_0d

        !-- Compute exner function:
        exner_function_0d = (p / p_0)**(rd_d_cp)

    end function exner_function_0d

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Compute the exner function for 1-D array arguments.
    !--------------------------------------------------------------------------------------------------!
    function exner_function_1d(p)

        implicit none

        real(wp), intent(IN), dimension(:) :: p !< pressure (Pa)

        real(wp), dimension(size(p)) :: exner_function_1d

        !-- Compute exner function:
        exner_function_1d = (p / p_0)**(rd_d_cp)

    end function exner_function_1d

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Compute the exner function for scalar arguments.
    !--------------------------------------------------------------------------------------------------!
    function exner_function_invers_0d(p)

        implicit none

        real(wp), intent(IN) :: p !< pressure (Pa)

        real(wp) :: exner_function_invers_0d

        !-- Compute exner function:
        exner_function_invers_0d = (p_0 / p)**(rd_d_cp)

    end function exner_function_invers_0d

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Compute the exner function for 1-D array arguments.
    !--------------------------------------------------------------------------------------------------!
    function exner_function_invers_1d(p)

        implicit none

        real(wp), intent(IN), dimension(:) :: p !< pressure (Pa)

        real(wp), dimension(size(p)) :: exner_function_invers_1d

        !-- Compute exner function:
        exner_function_invers_1d = (p_0 / p)**(rd_d_cp)

    end function exner_function_invers_1d

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Compute the barometric formula for scalar arguments. The calculation is based on the assumption
    !> of a polytropic atmosphere and neutral stratification, where the temperature lapse rate is g/cp.
    !--------------------------------------------------------------------------------------------------!
    function barometric_formula_0d(z, t_0, p_0)

        implicit none

        real(wp), intent(IN) :: z !< height (m)
        real(wp), intent(IN) :: t_0 !< temperature reference state (K)
        real(wp), intent(IN) :: p_0 !< surface pressure (Pa)

        real(wp) :: barometric_formula_0d

        !-- Compute barometric formula:
        barometric_formula_0d = p_0 * ((t_0 - g_d_cp * z) / t_0)**(cp_d_rd)

    end function barometric_formula_0d

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Compute the barometric formula for 1-D array arguments. The calculation is based on the
    !> assumption of a polytropic atmosphere and neutral stratification, where the temperature lapse
    !> rate is g/cp.
    !--------------------------------------------------------------------------------------------------!
    function barometric_formula_1d(z, t_0, p_0)

        implicit none

        real(wp), intent(IN), dimension(:) :: z !< height (m)
        real(wp), intent(IN) :: t_0 !< temperature reference state (K)
        real(wp), intent(IN) :: p_0 !< surface pressure (Pa)

        real(wp), dimension(size(z)) :: barometric_formula_1d

        !-- Compute barometric formula:
        barometric_formula_1d = p_0 * ((t_0 - g_d_cp * z) / t_0)**(cp_d_rd)

    end function barometric_formula_1d

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Compute relative humidty for an equilibrium thermodynamic state (i.e. bounds between 0 and 1)
    !> Scalar version
    !--------------------------------------------------------------------------------------------------!
    function get_relative_humidity_equilibrium_0d(vapor_content, temperature, pressure, air_density)

        implicit none

        real(wp), intent(IN) :: temperature !< thermodynamic temperature (i.e., exner x pt)
        real(wp), intent(IN) :: pressure !< pressure (hyp)
        real(wp), intent(IN) :: air_density !< rho_air
        real(wp), intent(IN) :: vapor_content !< q

        real(wp) :: get_relative_humidity_equilibrium_0d
        !-- Local variables
        real(wp) :: pressure_ratio_1 !< (hydrostatic / saturation vapor pressure) - 1
        real(wp) :: saturation_vapor_pressure !< saturation vapor pressure
        real(wp) :: vapor_content_1 !< 1- vapor content
        real(wp) :: vapor_ratio !< vapor content / (1 - vapor content)
        !-- Prevents singularity
        saturation_vapor_pressure = max(1.0e-12_wp, magnus(temperature))
        vapor_content_1 = max(1.0e-12_wp, (1.0_wp - vapor_content))
        !-- Vapor and pressure ratios
        vapor_ratio = vapor_content / vapor_content_1
        pressure_ratio_1 = (pressure / saturation_vapor_pressure) - 1.0_wp
        !-- Compute RH:
        get_relative_humidity_equilibrium_0d = max(0.0_wp, min(1.0_wp, &
                                                               (air_density / rd_d_rv) * &
                                                               pressure_ratio_1 * vapor_ratio))

    end function get_relative_humidity_equilibrium_0d

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Compute relative humidty for an equilibrium thermodynamic state (i.e. bounds between 0 and 1)
    !> 1D-array version
    !--------------------------------------------------------------------------------------------------!
    function get_relative_humidity_equilibrium_1d(vapor_content, temperature, pressure, air_density)

        implicit none

        real(wp), intent(IN), dimension(:) :: vapor_content !< q
        real(wp), intent(IN), dimension(:) :: temperature !< thermodynamic temperature (i.e., exner x pt)
        real(wp), intent(IN), dimension(:) :: pressure !< pressure (hyp)
        real(wp), intent(IN), dimension(:) :: air_density !< rho_air

        real(wp), dimension(size(vapor_content)) :: get_relative_humidity_equilibrium_1d
        !-- Local variables
        real(wp), dimension(size(vapor_content)) :: pressure_ratio_1 !< (p/es) - 1
        real(wp), dimension(size(vapor_content)) :: saturation_vapor_pressure !< es
        real(wp), dimension(size(vapor_content)) :: vapor_content_1 !< 1 - q
        real(wp), dimension(size(vapor_content)) :: vapor_ratio !< q / (1-q)
        !-- Prevent singularity
        saturation_vapor_pressure = max(1.0e-12_wp, magnus(temperature))
        vapor_content_1 = max(1.0e-12_wp, (1.0_wp - vapor_content))
        !-- Vapor and pressure ratios
        vapor_ratio = vapor_content / vapor_content_1
        pressure_ratio_1 = (pressure / saturation_vapor_pressure) - 1.0_wp
        !-- Compute RH:
        get_relative_humidity_equilibrium_1d = max(0.0_wp, min(1.0_wp, &
                                                               (air_density / rd_d_rv) * &
                                                               pressure_ratio_1 * vapor_ratio))

    end function get_relative_humidity_equilibrium_1d

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Compute relative humidty for a potentially supersaturated thermodynamic state (i.e., RH > 1)
    !> Scalar version
    !--------------------------------------------------------------------------------------------------!
    function get_relative_humidity_supersaturated_0d(vapor_content)

        implicit none

        real(wp), intent(IN) :: vapor_content !< q

        real(wp) :: get_relative_humidity_supersaturated_0d


        !-- placeholder for new implementation
        get_relative_humidity_supersaturated_0d = 0.0_wp * vapor_content

    end function get_relative_humidity_supersaturated_0d

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Compute relative humidty for a potentially supersaturated thermodynamic state (i.e., RH > 1)
    !> 1D array version
    !--------------------------------------------------------------------------------------------------!
    function get_relative_humidity_supersaturated_1d(vapor_content)

        implicit none

        real(wp), intent(IN), dimension(:) :: vapor_content !< q

        real(wp), dimension(size(vapor_content)) :: get_relative_humidity_supersaturated_1d

        !-- placeholder for new implementation
        get_relative_humidity_supersaturated_1d = 0.0_wp * vapor_content

    end function get_relative_humidity_supersaturated_1d

end module basic_constants_and_equations_mod
