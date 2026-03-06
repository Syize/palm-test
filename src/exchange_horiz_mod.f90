!> @file exchange_horiz.f90
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
!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Exchange of ghost point layers for subdomains (in parallel mode) and setting of cyclic lateral
!> boundary conditions for the total domain .
!--------------------------------------------------------------------------------------------------!
module exchange_horiz_mod
    use pegrid

    use kinds

    implicit none

    interface exchange_horiz
        module procedure exchange_horiz
    end interface exchange_horiz

    interface exchange_horiz_int
        module procedure exchange_horiz_int
    end interface exchange_horiz_int

    interface exchange_horiz_2d
        module procedure exchange_horiz_2d
    end interface exchange_horiz_2d

    interface exchange_horiz_2d_byte
        module procedure exchange_horiz_2d_byte
    end interface exchange_horiz_2d_byte

    interface exchange_horiz_2d_int
        module procedure exchange_horiz_2d_int
    end interface exchange_horiz_2d_int

    private

    public exchange_horiz, &
        exchange_horiz_int, &
        exchange_horiz_2d, &
        exchange_horiz_2d_byte, &
        exchange_horiz_2d_int


contains

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Exchange of ghost point layers for subdomains (in parallel mode) and setting of cyclic lateral
    !> boundary conditions for the total domain.
    !> This routine is for REAL 3d-arrays.
    !--------------------------------------------------------------------------------------------------!
    subroutine exchange_horiz(ar, nbgp_local, alternative_communicator, grid_level, mg_switch_to_pe0)

        use control_parameters, &
            ! use default value
            only: bc_lr_cyc, &
                  bc_ns_cyc


        use indices, &
            only: nxl, &
                  nxr, &
                  nyn, &
                  nys, &
                  nzb, &
                  nzt

        integer(iwp), optional, intent(IN) :: alternative_communicator !< alternative MPI communicator to be used
        integer(iwp), optional, intent(IN) :: grid_level !< multigrid grid level to be used

        ! #if defined( __parallel )
        !         integer(iwp) :: bufsize !< size of buffer for sending/receiving contiguous data
        ! #endif
        integer(iwp) :: communicator !< communicator that is used as argument in MPI calls
        integer(iwp) :: i !< loop index
        integer(iwp) :: j !< loop index
        integer(iwp) :: k !< loop index
        integer(iwp) :: l !< grid level to choose the dreived datatypes
        integer(iwp) :: left_pe !< id of left pe that is used as argument in MPI calls
        integer(iwp) :: nbgp_local !< number of ghost point layers
        integer(iwp) :: north_pe !< id of north pe that is used as argument in MPI calls
        integer(iwp) :: right_pe !< id of right pe that is used as argument in MPI calls
        integer(iwp) :: south_pe !< id of south pe that is used as argument in MPI calls

        logical, optional, intent(IN) :: mg_switch_to_pe0

        logical :: switch_to_pe0 !< local switch telling if total domain is gathered on one PE

        real(wp), dimension(nzb:nzt + 1, nys - nbgp_local:nyn + nbgp_local, &
                            nxl - nbgp_local:nxr + nbgp_local) :: ar !< 3d-array for which exchange is done


        !-- Set the communicator to be used
        if (present(alternative_communicator)) then
            !--    Alternative communicator is to be used
            communicator = communicator_configurations(alternative_communicator) % mpi_communicator
            left_pe = communicator_configurations(alternative_communicator) % pleft
            right_pe = communicator_configurations(alternative_communicator) % pright
            south_pe = communicator_configurations(alternative_communicator) % psouth
            north_pe = communicator_configurations(alternative_communicator) % pnorth

        else
            !--    Main communicator is to be used
            communicator = comm2d
            left_pe = pleft
            right_pe = pright
            south_pe = psouth
            north_pe = pnorth

        end if

        !-- Set the grid level.
        if (present(grid_level)) then
            l = grid_level
        else
            l = 0
        end if

        !-- Set the switch that tells if data of the toal domain is gathered on the PE
        if (present(mg_switch_to_pe0)) then
            switch_to_pe0 = mg_switch_to_pe0
        else
            switch_to_pe0 = .false.
        end if
        !-- Lateral boundary conditions in the non-parallel case.
        !-- Case dependent, because in GPU mode still not all arrays are on device. This workaround has to
        !-- be removed later. Also, since PGI compiler 12.5 has problems with array syntax, explicit loops
        !-- are used.
        if (present(alternative_communicator)) then
            if (alternative_communicator <= 2) then
                !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) DEFAULT(PRESENT) IF(lacc) IF(enable_openacc)
                do i = 1, nbgp_local
                    do j = nys - nbgp_local, nyn + nbgp_local
                        do k = nzb, nzt + 1
                            ar(k, j, nxl - nbgp_local - 1 + i) = ar(k, j, nxr - nbgp_local + i)
                        end do
                    end do
                end do
                !$ACC END PARALLEL LOOP
                !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) DEFAULT(PRESENT) IF(lacc) IF(enable_openacc)
                do i = 1, nbgp_local
                    do j = nys - nbgp_local, nyn + nbgp_local
                        do k = nzb, nzt + 1
                            ar(k, j, nxr + i) = ar(k, j, nxl - 1 + i)
                        end do
                    end do
                end do
                !$ACC END PARALLEL LOOP
            end if
        else
            if (bc_lr_cyc) then
                !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) DEFAULT(PRESENT) IF(lacc) IF(enable_openacc)
                do i = 1, nbgp_local
                    do j = nys - nbgp_local, nyn + nbgp_local
                        do k = nzb, nzt + 1
                            ar(k, j, nxl - nbgp_local - 1 + i) = ar(k, j, nxr - nbgp_local + i)
                        end do
                    end do
                end do
                !$ACC END PARALLEL LOOP
                !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) DEFAULT(PRESENT) IF(lacc) IF(enable_openacc)
                do i = 1, nbgp_local
                    do j = nys - nbgp_local, nyn + nbgp_local
                        do k = nzb, nzt + 1
                            ar(k, j, nxr + i) = ar(k, j, nxl - 1 + i)
                        end do
                    end do
                end do
                !$ACC END PARALLEL LOOP
            end if
        end if

        if (present(alternative_communicator)) then
            if (alternative_communicator == 1 .or. alternative_communicator == 3) then
                !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) DEFAULT(PRESENT) IF(lacc) IF(enable_openacc)
                do i = nxl - nbgp_local, nxr + nbgp_local
                    do j = 1, nbgp_local
                        do k = nzb, nzt + 1
                            ar(k, nys - nbgp_local - 1 + j, i) = ar(k, nyn - nbgp_local + j, i)
                        end do
                    end do
                end do
                !$ACC END PARALLEL LOOP
                !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) DEFAULT(PRESENT) IF(lacc) IF(enable_openacc)
                do i = nxl - nbgp_local, nxr + nbgp_local
                    do j = 1, nbgp_local
                        do k = nzb, nzt + 1
                            ar(k, nyn + j, i) = ar(k, nys - 1 + j, i)
                        end do
                    end do
                end do
                !$ACC END PARALLEL LOOP
            end if
        else
            if (bc_ns_cyc) then
                !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) DEFAULT(PRESENT) IF(lacc) IF(enable_openacc)
                do i = nxl - nbgp_local, nxr + nbgp_local
                    do j = 1, nbgp_local
                        do k = nzb, nzt + 1
                            ar(k, nys - nbgp_local - 1 + j, i) = ar(k, nyn - nbgp_local + j, i)
                        end do
                    end do
                end do
                !$ACC END PARALLEL LOOP
                !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) DEFAULT(PRESENT) IF(lacc) IF(enable_openacc)
                do i = nxl - nbgp_local, nxr + nbgp_local
                    do j = 1, nbgp_local
                        do k = nzb, nzt + 1
                            ar(k, nyn + j, i) = ar(k, nys - 1 + j, i)
                        end do
                    end do
                end do
                !$ACC END PARALLEL LOOP
            end if
        end if


        ! #endif
        ! call cpu_log(log_point_s(2), 'exchange_horiz', 'stop')

    end subroutine exchange_horiz
    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> @todo Missing subroutine description.
    !--------------------------------------------------------------------------------------------------!
    subroutine exchange_horiz_int(ar, nys_l, nyn_l, nxl_l, nxr_l, nzt_l, nbgp_local, type_xz_in, &
                                  type_yz_in, alternative_communicator)

        use control_parameters, &
            only: bc_lr_cyc, &
                  bc_ns_cyc

        use indices, &
            only: nzb

        integer(iwp), optional :: alternative_communicator !< alternative MPI communicator to be used

        integer(iwp) :: communicator !< communicator that is used as argument in MPI calls
        integer(iwp) :: left_pe !< id of left pe that is used as argument in MPI calls
        integer(iwp) :: nbgp_local !< number of ghost points
        integer(iwp) :: north_pe !< id of north pe that is used as argument in MPI calls
        integer(iwp) :: nxl_l !< local index bound at current grid level, left side
        integer(iwp) :: nxr_l !< local index bound at current grid level, right side
        integer(iwp) :: nyn_l !< local index bound at current grid level, north side
        integer(iwp) :: nys_l !< local index bound at current grid level, south side
        integer(iwp) :: nzt_l !< local index bound at current grid level, top
        integer(iwp) :: right_pe !< id of right pe that is used as argument in MPI calls
        integer(iwp) :: south_pe !< id of south pe that is used as argument in MPI calls
        integer(iwp) :: type_xz !< MPI datatype of exchange 3D data - left/right
        integer(iwp) :: type_yz !< MPI datatype of exchange 3D data - north/south

        integer(iwp), optional :: type_xz_in !< passed MPI datatype to exchange 3D data between left/right MPI ranks
        integer(iwp), optional :: type_yz_in !< passed MPI datatype to exchange 3D data between north/south MPI ranks

        integer(iwp), dimension(nzb:nzt_l + 1, nys_l - nbgp_local:nyn_l + nbgp_local, &
                                nxl_l - nbgp_local:nxr_l + nbgp_local) :: ar !< treated array

        !-- Set MPI datatype depending on the requested task.
        if (present(type_xz_in)) then
            type_xz = type_xz_in
        else
        ! #if defined( __parallel )
        !             type_xz = type_xz_int(0)
        ! #endif
        end if

        if (present(type_yz_in)) then
            type_yz = type_yz_in
        else
        ! #if defined( __parallel )
        !             type_yz = type_yz_int(0)
        ! #endif
        end if

        !-- Set the communicator to be used.
        if (present(alternative_communicator)) then
            !--    Alternative communicator is to be used.
            communicator = communicator_configurations(alternative_communicator) % mpi_communicator
            left_pe = communicator_configurations(alternative_communicator) % pleft
            right_pe = communicator_configurations(alternative_communicator) % pright
            south_pe = communicator_configurations(alternative_communicator) % psouth
            north_pe = communicator_configurations(alternative_communicator) % pnorth

        else
            !--    Main communicator is to be used.
            communicator = comm2d
            left_pe = pleft
            right_pe = pright
            south_pe = psouth
            north_pe = pnorth

        end if
        !-- Lateral boundary conditions in the non-parallel case.
        if (present(alternative_communicator)) then
            if (alternative_communicator <= 2) then
                ar(:, :, nxl_l - nbgp_local:nxl_l - 1) = ar(:, :, nxr_l - nbgp_local + 1:nxr_l)
                ar(:, :, nxr_l + 1:nxr_l + nbgp_local) = ar(:, :, nxl_l:nxl_l + nbgp_local - 1)
            end if
        else
            if (bc_lr_cyc) then
                ar(:, :, nxl_l - nbgp_local:nxl_l - 1) = ar(:, :, nxr_l - nbgp_local + 1:nxr_l)
                ar(:, :, nxr_l + 1:nxr_l + nbgp_local) = ar(:, :, nxl_l:nxl_l + nbgp_local - 1)
            end if
        end if

        if (present(alternative_communicator)) then
            if (alternative_communicator == 1 .or. alternative_communicator == 3) then
                ar(:, nys_l - nbgp_local:nys_l - 1, :) = ar(:, nyn_l - nbgp_local + 1:nyn_l, :)
                ar(:, nyn_l + 1:nyn_l + nbgp_local, :) = ar(:, nys_l:nys_l + nbgp_local - 1, :)
            end if
        else
            if (bc_ns_cyc) then
                ar(:, nys_l - nbgp_local:nys_l - 1, :) = ar(:, nyn_l - nbgp_local + 1:nyn_l, :)
                ar(:, nyn_l + 1:nyn_l + nbgp_local, :) = ar(:, nys_l:nys_l + nbgp_local - 1, :)
            end if
        end if


    ! #endif
    end subroutine exchange_horiz_int

    ! Description:
    ! ------------
    !> Exchange of lateral (ghost) boundaries (parallel computers) and cyclic boundary conditions,
    !> respectively, for 2D-arrays.
    !>
    !> TODO: This routine requires adjustments to be used for GPUs. The simple check further below
    !>       is commented, because the routine is called for output purposes in average_3d_data, so that
    !>       the check would stop the run since many of the output quantities do not yet exist on the
    !>       device.
    !--------------------------------------------------------------------------------------------------!
    subroutine exchange_horiz_2d(ar)

        ! #if ! defined( __parallel )
        use control_parameters, &
            only: bc_lr_cyc, &
                  bc_ns_cyc


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

        real(wp) :: ar(nysg:nyng, nxlg:nxrg) !<

        !-- Lateral boundary conditions in the non-parallel case
        if (bc_lr_cyc) then
            ar(:, nxlg:nxl - 1) = ar(:, nxr - nbgp + 1:nxr)
            ar(:, nxr + 1:nxrg) = ar(:, nxl:nxl + nbgp - 1)
        end if

        if (bc_ns_cyc) then
            ar(nysg:nys - 1, :) = ar(nyn - nbgp + 1:nyn, :)
            ar(nyn + 1:nyng, :) = ar(nys:nys + nbgp - 1, :)
        end if


        ! #endif
        ! call cpu_log(log_point_s(13), 'exchange_horiz_2d', 'stop')

    end subroutine exchange_horiz_2d

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Exchange of lateral (ghost) boundaries (parallel computers) and cyclic boundary conditions,
    !> respectively, for 2D 8-bit integer arrays.
    !--------------------------------------------------------------------------------------------------!
    subroutine exchange_horiz_2d_byte(ar, nys_l, nyn_l, nxl_l, nxr_l, nbgp_local)

        ! #if ! defined( __parallel )
        use control_parameters, &
            only: bc_lr_cyc, &
                  bc_ns_cyc


        integer(iwp) :: nbgp_local !< number of ghost layers to be exchanged
        integer(iwp) :: nxl_l !< local index bound at current grid level, left side
        integer(iwp) :: nxr_l !< local index bound at current grid level, right side
        integer(iwp) :: nyn_l !< local index bound at current grid level, north side
        integer(iwp) :: nys_l !< local index bound at current grid level, south side

        integer(ibp), dimension(nys_l - nbgp_local:nyn_l + nbgp_local, &
                                nxl_l - nbgp_local:nxr_l + nbgp_local) :: ar !< treated array

        !-- Lateral boundary conditions in the non-parallel case
        if (bc_lr_cyc) then
            ar(:, nxl_l - nbgp_local:nxl_l - 1) = ar(:, nxr_l - nbgp_local + 1:nxr_l)
            ar(:, nxr_l + 1:nxr_l + nbgp_local) = ar(:, nxl_l:nxl_l + nbgp_local - 1)
        end if

        if (bc_ns_cyc) then
            ar(nys_l - nbgp_local:nys_l - 1, :) = ar(nyn_l + 1 - nbgp_local:nyn_l, :)
            ar(nyn_l + 1:nyn_l + nbgp_local, :) = ar(nys_l:nys_l - 1 + nbgp_local, :)
        end if


        ! #endif
        ! call cpu_log(log_point_s(13), 'exchange_horiz_2d', 'stop')

    end subroutine exchange_horiz_2d_byte

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Exchange of lateral (ghost) boundaries (parallel computers) and cyclic boundary conditions,
    !> respectively, for 2D 32-bit integer arrays.
    !--------------------------------------------------------------------------------------------------!
    subroutine exchange_horiz_2d_int(ar, nys_l, nyn_l, nxl_l, nxr_l, nbgp_local)

        ! #if ! defined( __parallel )
        use control_parameters, &
            only: bc_lr_cyc, &
                  bc_ns_cyc


        integer(iwp) :: nbgp_local !< number of ghost layers to be exchanged
        integer(iwp) :: nxl_l !< local index bound at current grid level, left side
        integer(iwp) :: nxr_l !< local index bound at current grid level, right side
        integer(iwp) :: nyn_l !< local index bound at current grid level, north side
        integer(iwp) :: nys_l !< local index bound at current grid level, south side

        integer(iwp), dimension(nys_l - nbgp_local:nyn_l + nbgp_local, &
                                nxl_l - nbgp_local:nxr_l + nbgp_local) :: ar !< treated array

        !-- Lateral boundary conditions in the non-parallel case
        if (bc_lr_cyc) then
            ar(:, nxl_l - nbgp_local:nxl_l - 1) = ar(:, nxr_l - nbgp_local + 1:nxr_l)
            ar(:, nxr_l + 1:nxr_l + nbgp_local) = ar(:, nxl_l:nxl_l + nbgp_local - 1)
        end if

        if (bc_ns_cyc) then
            ar(nys_l - nbgp_local:nys_l - 1, :) = ar(nyn_l + 1 - nbgp_local:nyn_l, :)
            ar(nyn_l + 1:nyn_l + nbgp_local, :) = ar(nys_l:nys_l - 1 + nbgp_local, :)
        end if


        ! #endif
        ! call cpu_log(log_point_s(13), 'exchange_horiz_2d', 'stop')

    end subroutine exchange_horiz_2d_int

end module exchange_horiz_mod
