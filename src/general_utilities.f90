!> @file general_utilities.f90
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
! Copyright 2021-2022 Leibniz Universitaet Hannover
!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> A collection of handy utilities.
!--------------------------------------------------------------------------------------------------!
module general_utilities

    use kinds

    use indices, &
        only: nx, &
              ny, &
              nz

    implicit none

    !-- Public functions
    public &
        cross_product, &
        gridpoint_id, &
        interpolate_linear, &
        normalize_vector

    interface cross_product
        module procedure cross_product
    end interface cross_product

    interface gridpoint_id
        module procedure gridpoint_id_2d
        module procedure gridpoint_id_3d
    end interface gridpoint_id

    interface interpolate_linear
        module procedure interpolate_linear_0d_wp
    end interface interpolate_linear

    interface normalize_vector
        module procedure normalize_vector
    end interface normalize_vector

contains

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! -------------------------------------------------------------------------------------------------!
    !> Compute the cross product of two vectors.
    !--------------------------------------------------------------------------------------------------!
    function cross_product(a, b) result(c)

        real(wp), dimension(3) :: a !< vector a
        real(wp), dimension(3) :: b !< vector b
        real(wp), dimension(3) :: c !< resulting vector c

        c(1) = a(2) * b(3) - a(3) * b(2)
        c(2) = a(3) * b(1) - a(1) * b(3)
        c(3) = a(1) * b(2) - a(2) * b(1)

    end function cross_product

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! -------------------------------------------------------------------------------------------------!
    !> This functions computes an unique ID for each given grid point (j,i).
    !--------------------------------------------------------------------------------------------------!
    function gridpoint_id_2d(j, i)

        integer(idp) :: gridpoint_id_2d !< grid point ID, unique for each (j,i)
        integer(iwp) :: i !< grid index in x-direction
        integer(iwp) :: j !< grid index in y-direction

        gridpoint_id_2d = i + (nx + 1) * j

    end function gridpoint_id_2d

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! -------------------------------------------------------------------------------------------------!
    !> This functions computes an unique ID for each given grid point (k,j,i).
    !--------------------------------------------------------------------------------------------------!
    function gridpoint_id_3d(k, j, i)

        integer(idp) :: gridpoint_id_3d !< grid point ID, unique for each (k,j,i)
        integer(iwp) :: i !< grid index in x-direction
        integer(iwp) :: j !< grid index in y-direction
        integer(iwp) :: k !< grid index in k-direction

        gridpoint_id_3d = i * (ny + 1) * (nz + 1) + j * (nz + 1) + k + 1

    end function gridpoint_id_3d

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    !--------------------------------------------------------------------------------------------------!
    !> Interpolation function, used to interpolate between two real-type scalar values, e.g. in space
    !> or time.
    !--------------------------------------------------------------------------------------------------!
    function interpolate_linear_0d_wp(var_x1, var_x2, fac)

        real(wp) :: fac !< interpolation factor
        real(wp) :: interpolate_linear_0d_wp !< interpolated value
        real(wp) :: var_x1 !< value at x1
        real(wp) :: var_x2 !< value at x2

        interpolate_linear_0d_wp = (1.0_wp - fac) * var_x1 + fac * var_x2

    end function interpolate_linear_0d_wp

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! -------------------------------------------------------------------------------------------------!
    !> Normalize a given vector of arbitrary dimension.
    !--------------------------------------------------------------------------------------------------!
    subroutine normalize_vector(a)

        integer(iwp) :: n !< running index

        real(wp) :: abs_value !< absolute value of given vector
        real(wp), dimension(:) :: a !< vector a

        abs_value = 0.0_wp
        do n = 1, size(a)
            abs_value = abs_value + a(n)**2
        end do
        abs_value = sqrt(abs_value)

        if (abs_value > 0.0_wp) a = a / abs_value

    end subroutine normalize_vector

end module general_utilities