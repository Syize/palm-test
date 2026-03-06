!> @file boundary_settings_mod.f90
!--------------------------------------------------------------------------------------------------!
! This file is part of the PALM model system.
!
! PALM is free software: you can redistribute it and/or modify it under the terms of the GNU General
! Public License as published by the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! PALM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
! implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
! Public License for more details.
!
! You should have received a copy of the GNU General Public License along with PALM. If not, see
! <http://www.gnu.org/licenses/>.
!
! Copyright 1997-2021 Leibniz Universitaet Hannover
!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> This module contains some general settings of specific boundary conditions used by various
!> modules.
!--------------------------------------------------------------------------------------------------!
module boundary_settings_mod

    use control_parameters, &
        only: bc_dirichlet_l, &
              bc_dirichlet_n, &
              bc_dirichlet_r, &
              bc_dirichlet_s, &
              bc_radiation_l, &
              bc_radiation_n, &
              bc_radiation_r, &
              bc_radiation_s

    use indices, &
        only: nbgp, &
              nxl, &
              nxlg, &
              nxr, &
              nxrg, &
              nyn, &
              nyng, &
              nys, &
              nysg

    use kinds

    implicit none

    interface set_lateral_neumann_bc
        module procedure set_lateral_neumann_bc_int1
        module procedure set_lateral_neumann_bc_int4
        module procedure set_lateral_neumann_bc_real
    end interface set_lateral_neumann_bc

!
!-- Public routines
    public set_lateral_neumann_bc

contains

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set Neumann boundary conditions at lateral boundaries for 1-byte 2D integer arrays.
!--------------------------------------------------------------------------------------------------!
    subroutine set_lateral_neumann_bc_int1(ar)

        integer(iwp) :: i !< running index over ghost points

        integer(ibp) :: ar(nysg:nyng, nxlg:nxrg) !< treated array

!
!-- Neumann-conditions at inflow/outflow/nested boundaries.
        if (bc_dirichlet_l .or. bc_radiation_l) then
            do i = nbgp, 1, -1
                ar(:, nxl - i) = ar(:, nxl)
            end do
        end if
        if (bc_dirichlet_r .or. bc_radiation_r) then
            do i = 1, nbgp
                ar(:, nxr + i) = ar(:, nxr)
            end do
        end if
        if (bc_dirichlet_s .or. bc_radiation_s) then
            do i = nbgp, 1, -1
                ar(nys - i, :) = ar(nys, :)
            end do
        end if
        if (bc_dirichlet_n .or. bc_radiation_n) then
            do i = 1, nbgp
                ar(nyn + i, :) = ar(nyn, :)
            end do
        end if

    end subroutine set_lateral_neumann_bc_int1

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set Neumann boundary conditions at lateral boundaries for 4-byte 2D integer arrays.
!--------------------------------------------------------------------------------------------------!
    subroutine set_lateral_neumann_bc_int4(ar)

        integer(iwp) :: i !< running index over ghost points

        integer(iwp) :: ar(nysg:nyng, nxlg:nxrg) !< treated array

!
!-- Neumann-conditions at inflow/outflow/nested boundaries.
        if (bc_dirichlet_l .or. bc_radiation_l) then
            do i = nbgp, 1, -1
                ar(:, nxl - i) = ar(:, nxl)
            end do
        end if
        if (bc_dirichlet_r .or. bc_radiation_r) then
            do i = 1, nbgp
                ar(:, nxr + i) = ar(:, nxr)
            end do
        end if
        if (bc_dirichlet_s .or. bc_radiation_s) then
            do i = nbgp, 1, -1
                ar(nys - i, :) = ar(nys, :)
            end do
        end if
        if (bc_dirichlet_n .or. bc_radiation_n) then
            do i = 1, nbgp
                ar(nyn + i, :) = ar(nyn, :)
            end do
        end if

    end subroutine set_lateral_neumann_bc_int4

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set Neumann boundary conditions at lateral boundaries for real-type 2D arrays.
!--------------------------------------------------------------------------------------------------!
    subroutine set_lateral_neumann_bc_real(ar)

        integer(iwp) :: i !< running index over ghost points

        real(wp) :: ar(nysg:nyng, nxlg:nxrg) !< treated array

!
!-- Neumann-conditions at inflow/outflow/nested boundaries.
        if (bc_dirichlet_l .or. bc_radiation_l) then
            do i = nbgp, 1, -1
                ar(:, nxl - i) = ar(:, nxl)
            end do
        end if
        if (bc_dirichlet_r .or. bc_radiation_r) then
            do i = 1, nbgp
                ar(:, nxr + i) = ar(:, nxr)
            end do
        end if
        if (bc_dirichlet_s .or. bc_radiation_s) then
            do i = nbgp, 1, -1
                ar(nys - i, :) = ar(nys, :)
            end do
        end if
        if (bc_dirichlet_n .or. bc_radiation_n) then
            do i = 1, nbgp
                ar(nyn + i, :) = ar(nyn, :)
            end do
        end if

    end subroutine set_lateral_neumann_bc_real

end module boundary_settings_mod
