!> @file topography_mod.f90
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
!> Setup of PALM's topography representation
!> @todo: Rearrange topo flag list
!> @todo: reference 3D buildings on top of orography is not tested and may need further improvement
!>        for steep slopes
!> @todo: Use more advanced setting of building type at filled holes
!--------------------------------------------------------------------------------------------------!
module topography_mod


    ! #if defined( __parallel )
    !     use MPI
    ! #endif
    use arrays_3d, &
        only: dzw, &
              dzu, &
              zu, &
              zw, ddzu, ddzw

    use boundary_settings_mod, &
        only: set_lateral_neumann_bc

    use control_parameters

    use exchange_horiz_mod, &
        only: exchange_horiz_2d, &
              exchange_horiz_2d_byte, &
              exchange_horiz_2d_int, &
              exchange_horiz_int

    use indices, &
        only: nbgp, &
              nx, &
              nxl, &
              nxlg, &
              nxr, &
              nxrg, &
              ny, &
              nys, &
              nysg, &
              nyn, &
              nyng, &
              nz, &
              nzb, &
              nzb_max, &
              nzt, &
              topo_min_level, &
              topo_top_ind, &
              topo_flags

    use general_utilities, &
        only: gridpoint_id

    use grid_variables, &
        only: dx, &
              dy, &
              zu_s_inner, &
              zw_w_inner

    use kinds

    use, intrinsic :: iso_fortran_env, only: int32, real32

    use netcdf_data_input_mod, &
        only: buildings_f, &
              building_id_f, &
              building_type_f, &
              char_fill, &
              char_lod, &
              check_existence, &
              close_input_file, &
              dims_xy, &
              get_attribute, &
              get_dimension_length, &
              get_variable, &
              init_model, &
              input_file_static, &
              input_pids_static, &
              inquire_num_variables, &
              inquire_variable_names, &
              list_building_ids, &
              open_read_file, &
              terrain_height_f

    use pegrid


    ! #if defined( __parallel )
    !     use pmc_handle_communicator, &
    !         only: pmc_get_model_info
    ! #endif
    ! use pmc_interface, &
    !     only: atmosphere_ocean_coupled_run, &
    !           nest_shift_z

    implicit none

    integer(int32), dimension(:), allocatable :: build_ids_unique !< list of building IDs in the domain
    integer(int32), dimension(:, :, :), allocatable :: topo !< temporary array used to initialize topography

    logical :: topo_read_all_domains !< control flag indicating whether topography is read from file in all domains

    save

    private
    !-- Public subroutines
    public init_topography

    interface init_topography
        module procedure init_topography
    end interface init_topography

contains

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! -------------------------------------------------------------------------------------------------!
    !> Routine that controls the basic topography initialization.
    !--------------------------------------------------------------------------------------------------!
    subroutine init_topography

        ! call location_message('Setup topography', 'start')
        !-- If a topography input file is available (static input file or ASCII), read the topography
        !-- information. Note, this is already done here in order to check for correct parameter settings
        !-- and file availability.
        call topography_input
        !-- Perform general initialization tasks and consistency checks.
        call topography_prepare_and_check
        !-- In case of cut-cell topography, further information is read from the static input file and
        !-- consistency checks are carried out. This is done in the cut-cell topography module.
        ! if (cut_cell_topography) call cct_input
        !-- Define the topography. Either generic topography is defined, or file topography is mapped onto
        !-- the numeric grid. In the later case, further processing steps are performed, e.g. small
        !-- cavities on the grid scale are filtered, buildings are mapped onto the underlying terrain if
        !-- required, etc.. In case of cut-cell topography, topography grid points are still defined
        !-- based on pre-processed rastered topography.
        call define_topography
        !-- Further define cut-cell topography on top of the already existing grid-cell topography.
        !-- At the moment, cut-cell topography simply replaces rastered topography, but this
        !-- might change in future.
        ! if (cut_cell_topography) call cct_define_topography(topo)
        !-- Set flags to mark the topography features on the grid.
        call topography_set_flags
        !-- Determine further topography grid indices, e.g. to output profile data, or to control
        !-- the degradation of the advection terms.
        call topography_set_indices

        ! call location_message('Setup topography', 'finished')

    end subroutine init_topography

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! -------------------------------------------------------------------------------------------------!
    !> Reads topography information from file or sets generic topography. Moreover, all
    !> topography-relevant topography arrays are initialized, and grid flags are set.
    !--------------------------------------------------------------------------------------------------!
    subroutine define_topography

        integer(int32) :: bh !< temporary vertical index of building height
        integer(int32) :: ch !< temporary vertical index for canyon height
        integer(int32) :: hv_in !< heavyside function to model inner tunnel surface
        integer(int32) :: i !< index variable along x
        integer(int32) :: index_left_bwall !< index for left building wall
        integer(int32) :: index_north_bwall !< index for north building wall
        integer(int32) :: index_right_bwall !< index for right building wall
        integer(int32) :: index_south_bwall !< index for south building wall
        integer(int32) :: index_left_cwall !< index for left canyon wall
        integer(int32) :: index_north_cwall !< index for north canyon wall
        integer(int32) :: index_right_cwall !< index for right canyon wall
        integer(int32) :: index_south_cwall !< index for south canyon wall
        integer(int32) :: j !< index variable along y
        integer(int32) :: k !< index variable along z
        integer(int32) :: ngp_bx !< grid point number of building size along x
        integer(int32) :: ngp_by !< grid point number of building size along y
        integer(int32) :: ngp_cx !< grid point number of canyon size along x
        integer(int32) :: ngp_cy !< grid point number of canyon size along y
        integer(int32) :: hv_out !< heavyside function to model outer tunnel surface
        integer(int32) :: td !< tunnel wall depth
        integer(int32) :: th !< height of outer tunnel wall
        integer(int32) :: txe_in !< end position of inner tunnel wall in x
        integer(int32) :: txe_out !< end position of outer tunnel wall in x
        integer(int32) :: txs_in !< start position of inner tunnel wall in x
        integer(int32) :: txs_out !< start position of outer tunnel wall in x
        integer(int32) :: tye_in !< end position of inner tunnel wall in y
        integer(int32) :: tye_out !< end position of outer tunnel wall in y
        integer(int32) :: tys_in !< start position of inner tunnel wall in y
        integer(int32) :: tys_out !< start position of outer tunnel wall in y

        integer(int32), dimension(:, :), allocatable :: nzb_local !< index for topography top at cell-center

        logical :: root_model !< flag to indicate if root or child

        !-- Check, if 'read_from_file' has been chosen for all domains (parents/childs) and set a
        !-- respective flag. Only then lateron global reduction operations (e.g. for getting the lowest
        !-- topography throughout all domains) will be allowed.
        topo_read_all_domains = (topography == 'read_from_file')
        ! #if defined( __parallel )
        !         call MPI_ALLREDUCE(MPI_IN_PLACE, topo_read_all_domains, 1, MPI_LOGICAL, MPI_LAND, &
        !                            MPI_COMM_WORLD, ierr)
        ! #endif
        !-- Define topography and set respective flags. Topography is either flat (only grid point at
        !-- k=nzb is flagged), generically defined (block, canyon or tunnel), or it can be read from file.
        !-- In case of elevated childs, the lowest grid point at k=nzb may belong to the atmosphere, so
        !-- this point is not flagged.
        select case (trim(topography))

        ! case ('flat')
        !     !--       Initialize 3D topography array, used later for initializing flags.
        !     !--       In case of vertically shifted nests, the boundary vertical coordinate does not have its
        !     !--       standard location 0.0, meaning that the respective boundary is "open".
        !     !--       Childs in ocean_mode always have an open bottom boundary.
        !     root_model = .true.
        !     ! #if defined( __parallel )
        !     !             if (nested_run) call pmc_get_model_info(root_model=root_model)
        !     ! #endif
        !     if (nest_shift_z == 0.0_real32 .and. .not. (ocean_mode .and. .not. root_model)) then
        !         topo(nzb, :, :) = ibset(topo(nzb, :, :), 0)
        !     end if

        case ('closed_channel')
            !--       Initialilize 3D topography array, used later for initializing flags.
            topo(nzb, :, :) = ibset(topo(nzb, :, :), 0)

        case ('single_building')
            !--       Single rectangular building, by default centered in the middle of the
            !--       total domain.
            ngp_bx = nint(building_length_x / dx)
            ngp_by = nint(building_length_y / dy)
            bh = minloc(abs(zw - building_height), 1) - 1
            if (abs(zw(bh) - building_height) == abs(zw(bh + 1) - building_height)) bh = bh + 1
            if (building_wall_left == 9999999.9_real32) then
                building_wall_left = (nx + 1 - ngp_bx) / 2 * dx
            end if
            index_left_bwall = nint(building_wall_left / dx)
            index_right_bwall = index_left_bwall + ngp_bx

            if (building_wall_south == 9999999.9_real32) then
                building_wall_south = (ny + 1 - ngp_by) / 2 * dy
            end if
            index_south_bwall = nint(building_wall_south / dy)
            index_north_bwall = index_south_bwall + ngp_by

            !--       Building size has to meet some requirements.
            if ((index_left_bwall < 1) .or. (index_right_bwall > nx - 1) .or. &
                (index_right_bwall < index_left_bwall + 3) .or. &
                (index_south_bwall < 1) .or. (index_north_bwall > ny - 1) .or. &
                (index_north_bwall < index_south_bwall + 3)) &
                then
                write (*, *) 'inconsistent building parameters:', &
                    '&index_left_bwall=', index_left_bwall, &
                    'index_right_bwall=', index_right_bwall, &
                    'index_south_bwall=', index_south_bwall, &
                    'index_north_bwall=', index_north_bwall, &
                    'nx=', nx, 'ny=', ny
                ! call message('topography_mod', 'PAC0320', 1, 2, 0, 6, 0)
            end if

            allocate (nzb_local(nysg:nyng, nxlg:nxrg))
            nzb_local = 0
            !--       Define the building.
            if (index_left_bwall <= nxr .and. index_right_bwall >= nxl .and. &
                index_south_bwall <= nyn .and. index_north_bwall >= nys) &
                then
                nzb_local(max(nys, index_south_bwall):min(nyn, index_north_bwall), &
                          max(nxl, index_left_bwall):min(nxr, index_right_bwall)) = bh
            end if
            !--       Set bit array on basis of nzb_local.
            do i = nxl, nxr
                do j = nys, nyn
                    topo(nzb:nzb_local(j, i), j, i) = ibset(topo(nzb:nzb_local(j, i), j, i), 0)
                end do
            end do

            deallocate (nzb_local)
            !--       Exchange ghost points. Further, in case of non-cyclic boundary conditions Neumann BC
            !--       are set for the topography, i.e. it is assumed that buildings continue with same height
            !--       outside the gloabl domain.
            call exchange_horiz_int(topo, nys, nyn, nxl, nxr, nzt, nbgp)
            call topography_set_non_cyc_bc(topo)

        case ('single_street_canyon')
            !--       Single quasi-2D street canyon of infinite length in x or y direction.
            !--       The canyon is centered in the other direction by default.
            if (canyon_width_x /= 9999999.9_real32) then
                !--          Street canyon in y direction
                ngp_cx = nint(canyon_width_x / dx)
                if (canyon_wall_left == 9999999.9_real32) then
                    canyon_wall_left = (nx + 1 - ngp_cx) / 2 * dx
                end if
                index_left_cwall = nint(canyon_wall_left / dx)
                index_right_cwall = index_left_cwall + ngp_cx
            elseif (canyon_width_y /= 9999999.9_real32) then
                !--          Street canyon in x direction
                ngp_cy = nint(canyon_width_y / dy)
                if (canyon_wall_south == 9999999.9_real32) then
                    canyon_wall_south = (ny + 1 - ngp_cy) / 2 * dy
                end if
                index_south_cwall = nint(canyon_wall_south / dy)
                index_north_cwall = index_south_cwall + ngp_cy

            else

                print *, 'no street canyon width given'
                ! call message('topography_mod', 'PAC0321', 1, 2, 0, 6, 0)

            end if

            ch = minloc(abs(zw - canyon_height), 1) - 1
            if (abs(zw(ch) - canyon_height) == abs(zw(ch + 1) - canyon_height)) ch = ch + 1
            dp_level_ind_b = ch
            !--       Street canyon size has to meet some requirements.
            if (canyon_width_x /= 9999999.9_real32) then
                if ((index_left_cwall < 1) .or. (index_right_cwall > nx - 1) .or. (ngp_cx < 3)) &
                    then
                    write (*, *) 'inconsistent canyon parameters:', &
                        '&index_left_cwall=', index_left_cwall, &
                        ' index_right_cwall=', index_right_cwall, &
                        ' ngp_cx=', ngp_cx, ' ch=', ch, ' nx=', nx, ' ny=', ny
                    ! call message('topography_mod', 'PAC0322', 1, 2, 0, 6, 0)
                end if
            elseif (canyon_width_y /= 9999999.9_real32) then
                if ((index_south_cwall < 1) .or. (index_north_cwall > ny - 1) .or. (ngp_cy < 3)) &
                    then
                    write (*, *) 'inconsistent canyon parameters:', &
                        '&index_south_cwall=', index_south_cwall, &
                        ' index_north_cwall=', index_north_cwall, &
                        ' ngp_cy=', ngp_cy, ' ch=', ch, ' nx=', nx, ' ny=', ny
                    ! call message('topography_mod', 'PAC0323', 1, 2, 0, 6, 0)
                end if
            end if
            if (canyon_width_x /= 9999999.9_real32 .and. canyon_width_y /= 9999999.9_real32) then
                print *, 'inconsistent canyon parameters:' // &
                                 '&street canyon can only be oriented either in x- or in y-direction'
                ! call message('topography_mod', 'PAC0324', 1, 2, 0, 6, 0)
            end if

            allocate (nzb_local(nysg:nyng, nxlg:nxrg))
            nzb_local = ch
            if (canyon_width_x /= 9999999.9_real32) then
                if (index_left_cwall <= nxr .and. index_right_cwall >= nxl) &
                    nzb_local(:, max(nxl, index_left_cwall + 1):min(nxr, index_right_cwall - 1)) = 0
            elseif (canyon_width_y /= 9999999.9_real32) then
                if (index_south_cwall <= nyn .and. index_north_cwall >= nys) &
                    nzb_local(max(nys, index_south_cwall + 1):min(nyn, index_north_cwall - 1), :) = 0
            end if
            !--       Set bit array on basis of nzb_local
            do i = nxl, nxr
                do j = nys, nyn
                    topo(nzb:nzb_local(j, i), j, i) = ibset(topo(nzb:nzb_local(j, i), j, i), 0)
                end do
            end do
            deallocate (nzb_local)
            !--       Exchange ghost points. Further, in case of non-cyclic boundary conditions Neumann BC
            !--       are set for the topography, i.e. it is assumed that buildings continue with same height
            !--       outside the gloabl domain.
            call exchange_horiz_int(topo, nys, nyn, nxl, nxr, nzt, nbgp)
            call topography_set_non_cyc_bc(topo)

        case ('tunnel')
            !--       Initialize surface height with zero and set lowest model grid point to topography.
            topo(nzb, :, :) = ibset(topo(nzb, :, :), 0)
            !--       Tunnel height.
            if (tunnel_height == 9999999.9_real32) then
                th = zw(int(0.2 * nz))
            else
                th = tunnel_height
            end if
            !--       Tunnel-wall depth.
            if (tunnel_wall_depth == 9999999.9_real32) then
                td = max(dx, dy, dz(1))
            else
                td = tunnel_wall_depth
            end if
            !--       Check for tunnel width
            if (tunnel_width_x == 9999999.9_real32 .and. tunnel_width_y == 9999999.9_real32) then
                print *, 'No tunnel width is given'
                ! call message('topography_mod', 'PAC0325', 1, 2, 0, 6, 0)
            end if
            if (tunnel_width_x /= 9999999.9_real32 .and. tunnel_width_y /= 9999999.9_real32) then
                print *, 'inconsistent tunnel parameters tunnel_width_x / tunnel_width_y'
                ! call message('topography_mod', 'PAC0326', 1, 2, 0, 6, 0)
            end if
            !--       Check for too small tunnel width in x- and y-direction
            if (tunnel_width_x /= 9999999.9_real32 .and. &
                tunnel_width_x - 2.0_real32 * td <= 2.0_real32 * dx) then
                print *, 'tunnel_width_x too small'
                ! call message('topography_mod', 'PAC0327', 1, 2, 0, 6, 0)
            end if
            if (tunnel_width_y /= 9999999.9_real32 .and. &
                tunnel_width_y - 2.0_real32 * td <= 2.0_real32 * dy) then
                print *, 'tunnel_width_y too small'
                ! call message('topography_mod', 'PAC0328', 1, 2, 0, 6, 0)
            end if
            !--       Check for too large tunnel width. Tunnel axis along y.
            if (tunnel_width_x /= 9999999.9_real32) then
                if (tunnel_width_x > (nx + 1) * dx) then
                    print *, 'tunnel_width_x too large'
                    ! call message('topography_mod', 'PAC0329', 1, 2, 0, 6, 0)
                end if

                txs_out = int((nx + 1) * 0.5_real32 * dx - tunnel_width_x * 0.5_real32)
                txe_out = int((nx + 1) * 0.5_real32 * dx + tunnel_width_x * 0.5_real32)
                txs_in = int((nx + 1) * 0.5_real32 * dx - (tunnel_width_x * 0.5_real32 - td))
                txe_in = int((nx + 1) * 0.5_real32 * dx + (tunnel_width_x * 0.5_real32 - td))

                tys_out = int((ny + 1) * 0.5_real32 * dy - tunnel_length * 0.5_real32)
                tye_out = int((ny + 1) * 0.5_real32 * dy + tunnel_length * 0.5_real32)
                tys_in = tys_out
                tye_in = tye_out
            end if
            !--       Tunnel axis along x.
            if (tunnel_width_y /= 9999999.9_real32) then
                if (tunnel_width_y > (ny + 1) * dy) then
                    print *, 'tunnel_width_y too large'
                    ! call message('topography_mod', 'PAC0330', 1, 2, 0, 6, 0)
                end if

                txs_out = int((nx + 1) * 0.5_real32 * dx - tunnel_length * 0.5_real32)
                txe_out = int((nx + 1) * 0.5_real32 * dx + tunnel_length * 0.5_real32)
                txs_in = txs_out
                txe_in = txe_out

                tys_out = int((ny + 1) * 0.5_real32 * dy - tunnel_width_y * 0.5_real32)
                tye_out = int((ny + 1) * 0.5_real32 * dy + tunnel_width_y * 0.5_real32)
                tys_in = int((ny + 1) * 0.5_real32 * dy - (tunnel_width_y * 0.5_real32 - td))
                tye_in = int((ny + 1) * 0.5_real32 * dy + (tunnel_width_y * 0.5_real32 - td))
            end if

            do i = nxl, nxr
                do j = nys, nyn
                    !--             Use heaviside function to model outer tunnel surface.
                    hv_out = th * 0.5_real32 * ((sign(1.0_real32, i * dx - txs_out) + 1.0_real32) &
                                            - (sign(1.0_real32, i * dx - txe_out) + 1.0_real32))

                    hv_out = hv_out * 0.5_real32 * ((sign(1.0_real32, j * dy - tys_out) + 1.0_real32) &
                                                - (sign(1.0_real32, j * dy - tye_out) + 1.0_real32))
                    !--             Use heaviside function to model inner tunnel surface.
                    hv_in = (th - td) * 0.5_real32 * ((sign(1.0_real32, i * dx - txs_in) + 1.0_real32) &
                                                  - (sign(1.0_real32, i * dx - txe_in) + 1.0_real32))

                    hv_in = hv_in * 0.5_real32 * ((sign(1.0_real32, j * dy - tys_in) + 1.0_real32) &
                                              - (sign(1.0_real32, j * dy - tye_in) + 1.0_real32))

                    if (hv_out - hv_in == 0.0_real32) then
                        !--                Set flags at x-y-positions without any tunnel surface.
                        topo(nzb + 1:nzt + 1, j, i) = ibclr(topo(nzb + 1:nzt + 1, j, i), 0)
                    else
                        !--                Set flags at x-y-positions with tunnel surfaces.
                        do k = nzb + 1, nzt + 1
                            !--                   Inner tunnel.
                            if (hv_out - hv_in == th) then
                                if (zw(k) <= hv_out) then
                                    topo(k, j, i) = ibset(topo(k, j, i), 0)
                                else
                                    topo(k, j, i) = ibclr(topo(k, j, i), 0)
                                end if
                            end if
                            !--                   Lateral tunnel walls
                            if (hv_out - hv_in == td) then
                                if (zw(k) <= hv_in) then
                                    topo(k, j, i) = ibclr(topo(k, j, i), 0)
                                elseif (zw(k) > hv_in .and. zw(k) <= hv_out) then
                                    topo(k, j, i) = ibset(topo(k, j, i), 0)
                                elseif (zw(k) > hv_out) then
                                    topo(k, j, i) = ibclr(topo(k, j, i), 0)
                                end if
                            end if
                        end do
                    end if
                end do
            end do
            !--       Exchange ghost points. Further, in case of non-cyclic boundary conditions Neumann BC
            !--       are set for the topography, i.e. it is assumed that buildings continue with same height
            !--       outside the gloabl domain.
            call exchange_horiz_int(topo, nys, nyn, nxl, nxr, nzt, nbgp)
            call topography_set_non_cyc_bc(topo)

        case ('read_from_file')
            !--       Note, topography information has already been read.
            !--       If required, further process topography, i.e. reference buildings on top of orography and
            !--       set temporary 3D topography array, which is used later to set grid flags. Calling of this
            !--       routine is also required in case of ASCII input, even though no distinction between
            !--       terrain- and building height is made in this case.
            call process_topography
            !--       Filter holes resolved by only one grid-point.
            ! if (.not. cut_cell_topography) call filter_topography
            call filter_topography
            !--       Exchange ghost points. Further, in case of non-cyclic boundary conditions Neumann BC
            !--       are set for the topography, i.e. it is assumed that buildings continue with same height
            !--       outside the gloabl domain.
            call exchange_horiz_int(topo, nys, nyn, nxl, nxr, nzt, nbgp)
            call topography_set_non_cyc_bc(topo)

        case DEFAULT
            !--       The DEFAULT case is reached either if the parameter topography contains a wrong character
            !--       string or if the user has defined a special case in the user interface. There, the
            !--       subroutine user_init_grid checks which of these two conditions applies.
            ! call user_init_grid(topo)
            ! if (.not. cut_cell_topography) call filter_topography
            call filter_topography
            !--       Exchange ghost points. Further, in case of non-cyclic boundary conditions Neumann BC
            !--       are set for the topography, i.e. it is assumed that buildings continue with same height
            !--       outside the gloabl domain.
            call exchange_horiz_int(topo, nys, nyn, nxl, nxr, nzt, nbgp)
            call topography_set_non_cyc_bc(topo)

        end select

        !-- Consistency checks and index array initialization are only required for non-flat topography.
        ! if (trim(topography) /= 'flat' .and. .not. cut_cell_topography) then
        if (trim(topography) /= 'flat') then
            !--    In case of non-flat topography, check whether the convention how to define the topography
            !--    grid has been set correctly, or whether the default is applicable. If this is not possible,
            !--    abort.
            if (trim(topography_grid_convention) == ' ') then
                if (trim(topography) /= 'closed_channel' .and. &
                    trim(topography) /= 'single_building' .and. &
                    trim(topography) /= 'single_street_canyon' .and. &
                    trim(topography) /= 'tunnel' .and. &
                    trim(topography) /= 'read_from_file') &
                    then
                    !--          The default value is not applicable here, because it is only valid for the four
                    !--          standard cases 'single_building', 'single_street_canyon', 'tunnel' and 'read_from_file'.
                    print *, 'missing value for topography_grid_convention'
                    ! call message('topography_mod', 'PAC0331', 1, 2, 0, 6, 0)
                else
                    !--          The default value is applicable here. Set convention according to topography.
                    if (trim(topography) == 'single_building' .or. &
                        trim(topography) == 'single_street_canyon') &
                        then
                        topography_grid_convention = 'cell_edge'
                    elseif (trim(topography) == 'read_from_file' .or. trim(topography) == 'tunnel') &
                        then
                        topography_grid_convention = 'cell_center'
                    end if
                end if
            elseif (trim(topography_grid_convention) /= 'cell_edge' .and. &
                    trim(topography_grid_convention) /= 'cell_center') &
                then
                write (*, *) 'illegal value for topography_grid_convention: "' // &
                    trim(topography_grid_convention) // '"'
                ! call message('topography_mod', 'PAC0332', 1, 2, 0, 6, 0)
            end if

            if (topography_grid_convention == 'cell_edge') then
                !--       The topography as defined using the 'cell_edge' convention describes the actual total size
                !--       of topography which is defined at the cell edges where u=0 on the topography walls in
                !--       x-direction and v=0 on the topography walls in y-direction.
                !--       Therefore, the existence of topography in the grid center is now reduced by 1dx at the
                !--       left topography walls and by 1dy at the north topography walls to form the basis for
                !--       the grid-center flag.
                !--       Note, the reverse memory access (first j loop, then i loop) is absolutely required at
                !--       this point.
                do j = nys + 1, nyn + 1
                    do i = nxl - 1, nxr
                        do k = nzb, nzt + 1
                            if (.not. btest(topo(k, j, i), 0) .or. .not. btest(topo(k, j, i + 1), 0)) &
                                topo(k, j, i) = ibclr(topo(k, j, i), 0)
                        end do
                    end do
                end do
                call exchange_horiz_int(topo, nys, nyn, nxl, nxr, nzt, nbgp)

                do i = nxl, nxr + 1
                    do j = nys - 1, nyn
                        do k = nzb, nzt + 1
                            if (.not. btest(topo(k, j, i), 0) .or. .not. btest(topo(k, j + 1, i), 0)) &
                                topo(k, j, i) = ibclr(topo(k, j, i), 0)
                        end do
                    end do
                end do
                call exchange_horiz_int(topo, nys, nyn, nxl, nxr, nzt, nbgp)

            end if
        end if

    end subroutine define_topography

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! -------------------------------------------------------------------------------------------------!
    !> Reads orography and building information.
    !--------------------------------------------------------------------------------------------------!
    subroutine topography_input

        character(LEN=100), dimension(:), allocatable :: var_names !< variable names in static input file

        integer(int32) :: i !< running index along x-direction
        integer(int32) :: id_topo !< NetCDF id of topograhy input file
        integer(int32) :: ii !< running index for IO blocks
        integer(int32) :: io_status !< status after reading the ascii topo file
        integer(int32) :: j !< running index along y-direction
        integer(int32) :: k !< running index along z-direction
        integer(int32) :: num_vars !< number of variables in netcdf input file
        integer(int32) :: skip_n_rows !< counting variable to skip rows while reading topography file

        real(real32) :: dum !< dummy variable to skip columns while reading topography file

        type(dims_xy) :: dim_static !< data structure for x, y-dimension in static input file

        !-- CPU measurement.
        ! call cpu_log(log_point_s(83), 'NetCDF/ASCII input topo', 'start')
        !-- Input via palm-input data standard.
        if (input_pids_static) then
            ! #if defined ( __netcdf )
            !--    Open file in read-only mode.
            call open_read_file(trim(input_file_static)//trim(coupling_char), id_topo)
            !--    At first, inquire all variable names.
            !--    This will be used to check whether an input variable exists or not.
            call inquire_num_variables(id_topo, num_vars)
            !--    Allocate memory to store variable names and inquire them.
            allocate (var_names(1:num_vars))
            call inquire_variable_names(id_topo, var_names)
            !--    Read x, y - dimensions. Only required for consistency checks.
            call get_dimension_length(id_topo, dim_static % nx, 'x')
            call get_dimension_length(id_topo, dim_static % ny, 'y')
            allocate (dim_static % x(0:dim_static % nx - 1))
            allocate (dim_static % y(0:dim_static % ny - 1))
            call get_variable(id_topo, 'x', dim_static % x)
            call get_variable(id_topo, 'y', dim_static % y)
            !--    Check whether dimension size in input file matches the model dimensions.
            if (dim_static % nx - 1 /= nx) then
                write (*, *) 'static driver: horizontal dimension in x-direction (=', &
                    dim_static % nx - 1, ') does not match the respective model ', &
                    'dimension (=', nx, ')'
                ! ! call message('topography_mod', 'PAC0333', 1, 2, 0, 6, 0)
            end if
            if (dim_static % ny - 1 /= ny) then
                write (*, *) 'static driver: horizontal dimension in y-direction (=', &
                    dim_static % ny - 1, ') does not match the respective model ', &
                    'dimension (=', ny, ')'
                ! ! call message('topography_mod', 'PAC0333', 1, 2, 0, 6, 0)
            end if
            !--    Check if grid spacing of provided input data matches the respective grid spacing in the
            !--    model. The allowed tolerance is 0.1% of the respective model grid spacing.
            if (abs(dim_static % x(1) - dim_static % x(0) - dx) > 0.001_real32 * dx) then
                write (*, *) 'static driver: horizontal grid spacing in x-direction (=', &
                    dim_static % x(1) - dim_static % x(0), ') does not match the ', &
                    'respective model model grid spacing (=', dx, ')'
                ! ! call message('topography_mod', 'PAC0334', 1, 2, 0, 6, 0)
            end if
            if (abs(dim_static % y(1) - dim_static % y(0) - dy) > 0.001_real32 * dy) then
                write (*, *) 'static driver: horizontal grid spacing in y-direction (=', &
                    dim_static % y(1) - dim_static % y(0), ') does not match the ', &
                    'respective model model grid spacing (=', dy, ')'
                ! ! call message('topography_mod', 'PAC0334', 1, 2, 0, 6, 0)
            end if
            !--    Terrain height. First, get variable-related _fillvalue attribute.
            if (check_existence(var_names, 'zt')) then
                terrain_height_f % from_file = .true.
                call get_attribute(id_topo, char_fill, terrain_height_f % fill, .false., 'zt')
                !--       Input 2D terrain height.
                allocate (terrain_height_f % var(nys:nyn, nxl:nxr))

                call get_variable(id_topo, 'zt', terrain_height_f % var, nxl, nxr, nys, nyn, nbgp=0)
            else
                terrain_height_f % from_file = .false.
            end if

            !--    Read building height. First, read its _fillvalue attribute, as well as lod attribute.
            buildings_f % from_file = .false.
            if (check_existence(var_names, 'buildings_2d')) then
                buildings_f % from_file = .true.
                call get_attribute(id_topo, char_lod, buildings_f % lod, .false., 'buildings_2d')
                call get_attribute(id_topo, char_fill, buildings_f % fill1, .false., 'buildings_2d')

                !--       Read 2D buildings.
                if (buildings_f % lod == 1) then
                    allocate (buildings_f % var_2d(nys:nyn, nxl:nxr))
                    call get_variable(id_topo, 'buildings_2d', buildings_f % var_2d, nxl, nxr, nys, nyn, &
                                      nbgp=0)
                else
                    write (*, *) 'static driver: wrong netCDF attribute lod (=', &
                        buildings_f % lod, ') for buildings_2d'
                    ! ! call message('topography_mod', 'PAC0335', 1, 2, 0, 6, 0)
                end if
            end if
            !--    If available, also read 3D building information. If both are available, use 3D information.
            if (check_existence(var_names, 'buildings_3d')) then
                buildings_f % from_file = .true.
                call get_attribute(id_topo, char_lod, buildings_f % lod, .false., 'buildings_3d')
                call get_attribute(id_topo, char_fill, buildings_f % fill2, .false., 'buildings_3d')
                call get_dimension_length(id_topo, buildings_f % nz, 'z')
                !--       Read 3D buildings
                if (buildings_f % lod == 2) then
                    allocate (buildings_f % z(nzb:buildings_f % nz - 1))
                    call get_variable(id_topo, 'z', buildings_f % z)
                    !--          Check if building information is consistent to numeric grid.
                    if (buildings_f % nz > size(zu)) then
                        write (*, *) 'static driver: too much data points (=', &
                            buildings_f % nz, ') along the vertical coordinate for', &
                            ' 3d building data (maximum allowed=', size(zu), ')'
                        ! ! call message('topography_mod', 'PAC0336', 2, 2, 0, 6, 0)
                    end if
                    if (any(abs(buildings_f % z(0:buildings_f % nz - 1) - zu(0:buildings_f % nz - 1)) > &
                            0.001_real32 * minval(dz(1:number_dz)))) then
                        k = 0
                        do while (k <= buildings_f % nz - 1)
                            if (abs(buildings_f % z(k) - zu(k)) > &
                                0.001_real32 * minval(dz(1:number_dz))) exit
                            k = k + 1
                        end do
                        write (*, *) 'static driver: vertical coordinate do not match ', &
                            'numeric grid at z(', k, ') for 3d building data'
                        ! ! call message('topography_mod', 'PAC0337', 2, 2, 0, 6, 0)
                    end if

                    allocate (buildings_f % var_3d(nzb:buildings_f % nz - 1, nys:nyn, nxl:nxr))
                    buildings_f % var_3d = 0
                    call get_variable(id_topo, 'buildings_3d', buildings_f % var_3d, nxl, nxr, nys, nyn, 0, &
                                      buildings_f % nz - 1, nbgp=0)
                else
                    write (*, *) 'static driver: wrong netCDF attribute lod (=', &
                        buildings_f % lod, ') for buildings_3d'
                    ! ! call message('topography_mod', 'PAC0338', 1, 2, 0, 6, 0)
                end if
            end if
            !--    Read building IDs and its FillValue attribute. Further required for mapping buildings on top
            !--    of orography.
            if (check_existence(var_names, 'building_id')) then
                building_id_f % from_file = .true.
                call get_attribute(id_topo, char_fill, building_id_f % fill, .false., 'building_id')
                allocate (building_id_f % var(nysg:nyng, nxlg:nxrg))
                call get_variable(id_topo, 'building_id', building_id_f % var, nxl, nxr, nys, nyn, &
                                  nbgp=nbgp)
            else
                building_id_f % from_file = .false.
            end if
            !--    Read building_type and required attributes.
            if (check_existence(var_names, 'building_type')) then
                building_type_f % from_file = .true.
                call get_attribute(id_topo, char_fill, building_type_f % fill, .false., 'building_type')
                allocate (building_type_f % var(nysg:nyng, nxlg:nxrg))
                call get_variable(id_topo, 'building_type', building_type_f % var, nxl, nxr, nys, nyn, &
                                  nbgp=nbgp)
            else
                building_type_f % from_file = .false.
            end if
            !--    Close topography input file.
            call close_input_file(id_topo)
        ! #else
        !             continue
        ! #endif
        !-- ASCII input
        end if
        !-- End of CPU measurement.
        ! call cpu_log(log_point_s(83), 'NetCDF/ASCII input topo', 'stop')

    end subroutine topography_input

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! -------------------------------------------------------------------------------------------------!
    !> Performs preparatory tasks, such as computing a list of building IDs, general checks, etc..
    !--------------------------------------------------------------------------------------------------!
    subroutine topography_prepare_and_check

        real(real32) :: oro_min = 0.0_real32 !< minimum terrain height in entire model domain, used to reference terrain to zero

        !-- Check for minimum requirement to setup building topography. If buildings are provided, also an
        !-- ID and a type are required.
        !-- Note that performing this check in check_parameters will be too late (data will be used for
        !-- grid initialization before).
        if (input_pids_static) then
            if (buildings_f % from_file .and. .not. building_id_f % from_file) then
                print *, 'static driver: building ID is missing'
                ! ! call message('topography_mod', 'PAC0342', 1, 2, 0, 6, 0)
            end if
        end if

        if (terrain_height_f % from_file) then
            !--    Check orography for fill-values.
            !--    For the moment, give an error message. More advanced methods, e.g. a nearest neighbor
            !--    algorithm as used in GIS systems might be implemented later.
            !--    Note: This check must be placed here as terrain_height_f is altered later.
            if (any(terrain_height_f % var == terrain_height_f % fill)) then
                print *, 'static driver: fill value for variable zt found'
                ! ! call message('topography_mod', 'PAC0343', 2, 2, myid, 6, 0)
            end if
        else
            !--    In case no terrain height is provided by static input file, allocate array nevertheless and
            !--    set terrain height to 0, which simplifies topography initialization.
            allocate (terrain_height_f % var(nys:nyn, nxl:nxr))
            terrain_height_f % var = 0.0_real32
        end if
        !-- Finally, exchange 1 ghost point for building ID and type.
        !-- In case of non-cyclic boundary conditions set Neumann conditions at the lateral boundaries.
        if (building_id_f % from_file) then
            call exchange_horiz_2d_int(building_id_f % var, nys, nyn, nxl, nxr, nbgp)
            call set_lateral_neumann_bc(building_id_f % var)
        end if

        if (building_type_f % from_file) then
            call exchange_horiz_2d_byte(building_type_f % var, nys, nyn, nxl, nxr, nbgp)
            call set_lateral_neumann_bc(building_type_f % var)
        end if
        !-- Check for correct setting of the namelist parameter topography. If topography information is
        !-- read from file but topography = 'flat', initialization does not work properly.
        if ((buildings_f % from_file .or. terrain_height_f % from_file) .and. &
            trim(topography) /= 'read_from_file') &
            then
            print *, 'wrong setting topography = "' // trim(topography) // '"'
            ! call message('topography_mod', 'PAC0319', 1, 2, 0, 6, 0)
        end if
        !-- Check, if 'read_from_file' has been chosen for all domains (parents/childs) and set a
        !-- respective flag. Only then lateron global reduction operations (e.g. for getting the lowest
        !-- topography throughout all domains) will be allowed.
        topo_read_all_domains = (topography == 'read_from_file')
        ! #if defined( __parallel )
        !         call MPI_ALLREDUCE(MPI_IN_PLACE, topo_read_all_domains, 1, MPI_LOGICAL, MPI_LAND, &
        !                            MPI_COMM_WORLD, ierr)
        ! #endif
        !-- Determine array of all building IDs within the domain where each ID occurs only once.
        call list_building_ids(build_ids_unique)
        !-- Allocate 3D array to set topography.
        allocate (topo(nzb:nzt + 1, nysg:nyng, nxlg:nxrg))
        topo = 0
        !-- Reference lowest terrain height to zero. This ensures that first, non-required gird levels
        !-- (those which lie entirely below the minimum orography) are avoided, and second, that also
        !-- negative orography can be used within the input file.
        !-- Please note, in case of a nested run, the global minimum from all parent and childs needs to be
        !-- removed to avoid steep edges at the child-domain boundaries.
        !-- Moreover, please note, the global minimum topography is only substracted if topography is
        !-- defined in all coupled models. If it is defined only in some of the models, this step is
        !-- skipped, same as in coupled atmosphere-ocean runs where topography height relates to
        !-- different reference levels.
        ! if (topo_read_all_domains .and. .not. atmosphere_ocean_coupled_run) then
        if (topo_read_all_domains) then

            if (input_pids_static) then
                ! #if defined( __parallel )
                !                 call MPI_ALLREDUCE(minval(terrain_height_f % var), oro_min, 1, MPI_REAL, MPI_MIN, &
                !                                    MPI_COMM_WORLD, ierr)
                ! #else
                oro_min = minval(terrain_height_f % var)
                ! #endif
                terrain_height_f % var = terrain_height_f % var - oro_min
                !--       Update reference height used within output files
                init_model % origin_z = init_model % origin_z + oro_min
            !--    ASCII topography branch. In this case, in contrast to the static driver input, topography is
            !--    input via the variable buildings_2d, which is used as a dummy here. This case, the minimum
            !--    building height is subtracted from the building array.
            else
                ! #if defined( __parallel )
                !                 call MPI_ALLREDUCE(minval(buildings_f % var_2d), oro_min, 1, MPI_REAL, MPI_MIN, &
                !                                    MPI_COMM_WORLD, ierr)
                ! #else
                oro_min = minval(buildings_f % var_2d)
                ! #endif
                buildings_f % var_2d = buildings_f % var_2d - oro_min
                !--       Update reference height used within output files.
                init_model % origin_z = init_model % origin_z + oro_min
            end if

        end if

    end subroutine topography_prepare_and_check

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! -------------------------------------------------------------------------------------------------!
    !> Filter topography. This subroutine contains two filter steps. First, one-grid point wide
    !> structures are filled. More precisely, all fluid grid points that are surrounded by more than 4
    !> topography grid points in the x-,y-,z-direction are filled. This filtering is applied
    !> iteratively until no grid-point wide holes exit any more.
    !> Such holes are suspected to lead to velocity blow-ups as the continuity equation on a discrete
    !> grid cannot be fulfilled in such case.
    !> In a second step, enclosed narrow cavities are filled. These are also suspected to lead to
    !> numerical instabilities, in particular when surface scalar fluxes are not zero. This case, scalar
    !> values can increase to unrealistic levels since almost no mixing between the narrow cavity and
    !> its surroundings takes place. At the moment, enclosed cavities of up to 9 grid points are
    !> filtered.
    !> Attention: The first filter step removes narrow elongated structures like streets, while the
    !> second step only filters small fluid cavities that are completely surrounded by topography/
    !> buildings.
    !--------------------------------------------------------------------------------------------------!
    subroutine filter_topography

        integer(int32) :: bl !< search bound at left subdomain boundary
        integer(int32) :: bn !< search bound at north subdomain boundary
        integer(int32) :: br !< search bound at right subdomain boundary
        integer(int32) :: bs !< search bound at south subdomain boundary
        integer(int32) :: f !< running index over all fluid indices flagged to be filtered
        integer(int32) :: i !< running index along x-direction
        integer(int32) :: i_f !< grid index of filtered grid point in x-direction
        integer(int32) :: j !< running index along y-direction
        integer(int32) :: j_f !< grid index of filtered grid point in y-direction
        integer(int32) :: k !< running index along z-direction
        ! #if defined( __parallel )
        !         integer(int32) :: ngp_yz !< number of extended ghost points
        ! #endif
        integer(int32) :: num_cavity !< number of narrow cavities that have been filled
        integer(int32) :: num_hole !< number of holes (in topography) resolved by only one grid point detected during one sweep
        integer(int32) :: num_hole_total
        !< total number of holes (in topography) resolved by only one grid point detected during all sweeps
        integer(int32) :: num_wall !< number of surrounding vertical walls for a single grid point
        integer(int32) :: sweep !< counts how often the topography filter is apllied
        integer(int32) :: type_xz_ext_int
        !< derived MPI datatype for 3-D integer ghost-point exchange with extended number of ghost points - left / right
        integer(int32) :: type_yz_ext_int
        !< derived MPI datatype for 3-D integer ghost-point exchange with extended number of ghost points - south / north

        integer(int32), dimension(nys:nyn, nxl:nxr) :: cav_filled !< flag indiating whether a cavity has been filled or not

        integer(int32), dimension(:, :, :), allocatable :: topo_tmp !< temporary 3D-topography used to fill holes

        type filter_type
            integer(int32) :: num_gp !< number of fluid grid point in current trace
            integer(int32) :: num_thresh = 9 !< maximum number of grid points forming a cavity that will be filtered

            integer(int32), dimension(:), allocatable :: i !< grid index in x-direction indicating fluid grid point in current trace
            integer(int32), dimension(:), allocatable :: j !< grid index in y-direction indicating fluid grid point in current trace
            integer(idp), dimension(:), allocatable :: ij_id !< unique ID for each (ji)-pair
        end type filter_type

        type(filter_type) :: filter !< derived structure summarizing all fluid-grid point in current trace

        !-- Before checking for holes, set lateral boundary conditions for topography. After hole-filling,
        !-- boundary conditions must be set again. Several iterations are performed, in order to fill
        !-- holes which might emerge by applying the filling-algorithm itself.
        !-- Attention: Beside small holes, also narrow elongated structures without a length limit (like
        !-- streets with a width of only one gridpoint) will be filtered!
        allocate (topo_tmp(nzb:nzt + 1, nysg:nyng, nxlg:nxrg))
        topo_tmp = 0

        num_hole = 99999
        num_hole_total = 0
        sweep = 0
        do while (num_hole > 0)

            sweep = sweep + 1
            num_hole = 0
            call exchange_horiz_int(topo, nys, nyn, nxl, nxr, nzt, nbgp)
            !--    Exchange also building ID and type. Note, building_type is a one-byte variable.
            if (building_id_f % from_file) then
                call exchange_horiz_2d_int(building_id_f % var, nys, nyn, nxl, nxr, nbgp)
                call set_lateral_neumann_bc(building_id_f % var)
            end if
            if (building_type_f % from_file) then
                call exchange_horiz_2d_byte(building_type_f % var, nys, nyn, nxl, nxr, nbgp)
                call set_lateral_neumann_bc(building_type_f % var)
            end if

            topo_tmp = topo
            !--    In case of non-cyclic lateral boundaries, assume lateral boundary to be a solid wall. Thus,
            !--    intermediate spaces of one grid point between boundary and some topographic structure will be
            !--    filled.
            if (.not. bc_ns_cyc) then
                if (nys == 0) topo_tmp(:, -1, :) = ibset(topo_tmp(:, 0, :), 0)
                if (nyn == ny) topo_tmp(:, ny + 1, :) = ibset(topo_tmp(:, ny, :), 0)
            end if

            if (.not. bc_lr_cyc) then
                if (nxl == 0) topo_tmp(:, :, -1) = ibset(topo_tmp(:, :, 0), 0)
                if (nxr == nx) topo_tmp(:, :, nx + 1) = ibset(topo_tmp(:, :, nx), 0)
            end if

            num_hole = 0
            do i = nxl, nxr
                do j = nys, nyn
                    do k = nzb + 1, nzt
                        if (.not. btest(topo_tmp(k, j, i), 0)) then
                            num_wall = 0
                            if (btest(topo_tmp(k, j - 1, i), 0)) num_wall = num_wall + 1
                            if (btest(topo_tmp(k, j + 1, i), 0)) num_wall = num_wall + 1
                            if (btest(topo_tmp(k, j, i - 1), 0)) num_wall = num_wall + 1
                            if (btest(topo_tmp(k, j, i + 1), 0)) num_wall = num_wall + 1
                            if (btest(topo_tmp(k - 1, j, i), 0)) num_wall = num_wall + 1
                            if (btest(topo_tmp(k + 1, j, i), 0)) num_wall = num_wall + 1

                            if (num_wall >= 4) then
                                num_hole = num_hole + 1
                                !--                   Clear flag 0 and set special flag ( bit 15) to indicate that new topography
                                !--                   point is a result of filtering process.
                                topo(k, j, i) = ibset(topo(k, j, i), 0)
                                topo(k, j, i) = ibset(topo(k, j, i), 15)
                                !--                   If filled grid point is occupied by a building, classify it as building grid
                                !--                   point.
                                if (building_type_f % from_file) then
                                    if (building_type_f % var(j, i) /= building_type_f % fill .or. &
                                        building_type_f % var(j + 1, i) /= building_type_f % fill .or. &
                                        building_type_f % var(j - 1, i) /= building_type_f % fill .or. &
                                        building_type_f % var(j, i + 1) /= building_type_f % fill .or. &
                                        building_type_f % var(j, i - 1) /= building_type_f % fill) &
                                        then
                                        !--                         Set flag indicating building surfaces
                                        topo(k, j, i) = ibset(topo(k, j, i), 12)
                                        !--                         Set building_type and ID at this position if not already set. This is
                                        !--                         required for proper initialization of urban-surface energy balance
                                        !--                         solver.
                                        if (building_type_f % var(j, i) == building_type_f % fill) then

                                            if (building_type_f % var(j + 1, i) /= building_type_f % fill) then
                                                building_type_f % var(j, i) = building_type_f % var(j + 1, i)
                                                building_id_f % var(j, i) = building_id_f % var(j + 1, i)
                                            elseif (building_type_f % var(j - 1, i) /= building_type_f % fill) then
                                                building_type_f % var(j, i) = building_type_f % var(j - 1, i)
                                                building_id_f % var(j, i) = building_id_f % var(j - 1, i)
                                            elseif (building_type_f % var(j, i + 1) /= building_type_f % fill) then
                                                building_type_f % var(j, i) = building_type_f % var(j, i + 1)
                                                building_id_f % var(j, i) = building_id_f % var(j, i + 1)
                                            elseif (building_type_f % var(j, i - 1) /= building_type_f % fill) then
                                                building_type_f % var(j, i) = building_type_f % var(j, i - 1)
                                                building_id_f % var(j, i) = building_id_f % var(j, i - 1)
                                            end if
                                        end if
                                    end if
                                end if
                                !--                   If filled grid point is already classified as building everything is fine,
                                !--                   else classify this grid point as natural type grid point. This case, values
                                !--                   for the surface type are already set.
                                if (.not. btest(topo(k, j, i), 12)) then
                                    topo(k, j, i) = ibset(topo(k, j, i), 11)
                                end if
                            end if
                        end if
                    end do
                end do
            end do
            !--    Count the total number of holes, required for informative message.
            ! #if defined( __parallel )
            !             call MPI_ALLREDUCE(MPI_IN_PLACE, num_hole, 1, MPI_INTEGER, MPI_SUM, comm2d, ierr)
            ! #endif
            if (num_hole > 0) then
                num_hole_total = num_hole_total + num_hole
            else
                sweep = sweep - 1
            end if

        end do

        deallocate (topo_tmp)
        !-- Issue an informative message if any 1-grid point wide holes were filled.
        if (num_hole_total > 0) then
            write (*, '(A,I6,A,A,I2,A)') 'topography was filtered: ', num_hole_total, &
                ' hole(s) resolved by only one grid point were filled during ', &
                'initialization in ', sweep, ' sweep(s)'
            ! call message('topography_mod', 'PAC0344', 0, 0, 0, 6, 0)
        end if
        !-- Finally, exchange topo array again and if necessary set Neumann boundary condition in case of
        !-- non-cyclic lateral boundaries.
        call exchange_horiz_int(topo, nys, nyn, nxl, nxr, nzt, nbgp)
        call topography_set_non_cyc_bc(topo)
        !-- Exchange building ID and type. Note, building_type is an one-byte variable.
        if (building_id_f % from_file) then
            call exchange_horiz_2d_int(building_id_f % var, nys, nyn, nxl, nxr, nbgp)
            call set_lateral_neumann_bc(building_id_f % var)
        end if
        if (building_type_f % from_file) then
            call exchange_horiz_2d_byte(building_type_f % var, nys, nyn, nxl, nxr, nbgp)
            call set_lateral_neumann_bc(building_type_f % var)
        end if

        !-- On top of this 1-gridpoint cavity filtering, also larger structures are filtered, i.e.
        !-- enclosed cavities of less than 9 gridpoints.
        !-- Note, at the moment only courtyard-like cavities in the xy-plance are filtered, while
        !-- cavities in the yz- and xz-plane are not treated at the moment.
        !-- In a first step, search for all grid points that are somehow topography bounded in their
        !-- vicinity and store their indices. Do this layer by layer, which is the reason for the reverse
        !-- loop structure.
        !-- First of all, determine the search bounds at the subdomain boundaries. Note, for the filter
        !-- algorithm extended ghost layers are requied.
        bl = merge(nxl - filter % num_thresh, 0,.not. (bc_dirichlet_l .or. bc_radiation_l))
        br = merge(nxr + filter % num_thresh, nx + 1,.not. (bc_dirichlet_r .or. bc_radiation_r))
        bn = merge(nyn + filter % num_thresh, ny + 1,.not. (bc_dirichlet_n .or. bc_radiation_n))
        bs = merge(nys - filter % num_thresh, 0,.not. (bc_dirichlet_s .or. bc_radiation_s))

        !-- Define MPI-datatypes for extended ghost-point exchange. This is necessary to detect also
        !-- elongated cavities with length of up to 9 grid points.
        !-- However, extended ghost layers require the subdomains in x- and y-direction to hold at least
        !-- the same number of grid points. Else (at least at the moment), ghost point exchange will not
        !-- work. Hence, in case of smaller subdomains, the number of ghost points and the filter
        !-- threshold must be reduced. In unlucky cases this can result in the situation that elongated
        !-- cavities are not fully filtered.
        filter % num_thresh = min(filter % num_thresh, nyn - nys + 1, nxr - nxl + 1)




        ! #if defined( __parallel )
        !         ngp_yz = (nzt - nzb + 2) * (nyn - nys + 1 + 2 * filter % num_thresh)
        !         call MPI_TYPE_VECTOR(nxr - nxl + 1 + 2 * filter % num_thresh, filter % num_thresh * (nzt - nzb + 2), ngp_yz, &
        !                              MPI_INTEGER, type_xz_ext_int, ierr)
        !         call MPI_TYPE_COMMIT(type_xz_ext_int, ierr)
        !         call MPI_TYPE_VECTOR(filter % num_thresh, ngp_yz, ngp_yz, MPI_INTEGER, type_yz_ext_int, ierr)
        !         call MPI_TYPE_COMMIT(type_yz_ext_int, ierr)
        !         !-- Set a barrier so that all MPI datatypes are defined before the ghost-point exchange starts.
        !         !-- Without such a barrier, a MPI error concerning insufficiently large buffer size may occur.
        !         call MPI_BARRIER(comm2d, ierr)
        ! #endif
        allocate (topo_tmp(nzb:nzt + 1, nys - filter % num_thresh:nyn + filter % num_thresh, &
                           nxl - filter % num_thresh:nxr + filter % num_thresh))

        topo_tmp = 0
        topo_tmp(:, nys:nyn, nxl:nxr) = topo(:, nys:nyn, nxl:nxr)
        call exchange_horiz_int(topo_tmp, nys, nyn, nxl, nxr, nzt, filter % num_thresh, &
                                type_xz_ext_int, type_yz_ext_int)
        !-- Allocate arrays containing the flagged grid-indices and their corresponding grid-point ID
        !-- that are flagged to be potentially filtered.
        allocate (filter % i(filter % num_thresh))
        allocate (filter % j(filter % num_thresh))
        allocate (filter % ij_id(filter % num_thresh))

        cav_filled = 0
        do k = nzb + 1, nzt
            do i = nxl, nxr
                do j = nys, nyn
                    !--          Only employ filter algorithm to non-topography grid points..
                    if (.not. btest(topo_tmp(k, j, i), 0)) then
                        !--             Starting from each grid point, search
                        filter % num_gp = 1
                        filter % i(:) = -huge(1)
                        filter % j(:) = -huge(1)
                        filter % ij_id(:) = -huge(1)

                        filter % i(1) = i
                        filter % j(1) = j
                        filter % ij_id(1) = gridpoint_id(j, i)
                        !--             Search for fluid grid points in the surroundings of (j,i). Note, the search is
                        !--             done recursively. If a surrounding grid point is flagged as fluid, further extend
                        !--             the search in the surrounding of this flagged grid point, and so on. The search
                        !--             for fluid grid points will be stopped if a treshold value of fluid grid point
                        !--             is reached. This case, no topography filter is employed. However, if the number
                        !--             of continuously connected fluid grid points in the farther surrounding of (j,i) is
                        !--             less than this threshold number, topography will be filtered and the cavity is
                        !--             filled.
                        call trace_surrounding_gridpoints(k, j, i)

                        if (filter % num_gp < filter % num_thresh) then
                            cav_filled(j, i) = 1

                            do f = 1, filter % num_gp
                                i_f = filter % i(f)
                                j_f = filter % j(f)
                                !--                   Set topography bit at flagged grid points. Only set topography if the
                                !--                   k-1 level is also topography. By this condition, filling
                                !--                   of courtyards with lateral openings, urn-shaped cavities, or elevated
                                !--                   "bottlenecks", which are narrower in their upper levels should be avoided.
                                !--                   Same as for the 1 grid point hole filling algorithm above, set bit 4 to
                                !--                   indicate filtered grid points.
                                if (btest(topo_tmp(k - 1, j_f, i_f), 0)) then
                                    topo_tmp(k, j_f, i_f) = ibset(topo_tmp(k, j_f, i_f), 0)
                                    topo_tmp(k, j_f, i_f) = ibset(topo_tmp(k, j_f, i_f), 15)

                                    !--                      If filled grid point is occupied by a building, classify it as building
                                    !--                      grid point.
                                    if (building_type_f % from_file) then
                                        if (j_f >= nysg .and. j_f <= nyng .and. &
                                            i_f >= nxlg .and. i_f <= nxrg) then
                                            if (building_type_f % var(j_f, i_f) /= building_type_f % fill) then
                                                !--                               Set flag indicating building surfaces.
                                                topo_tmp(k, j_f, i_f) = ibset(topo_tmp(k, j_f, i_f), 12)
                                            end if
                                        end if
                                    end if
                                    !--                      If filled grid point is already classified as building everything is fine,
                                    !--                      else classify this grid point as natural type grid point. This case,
                                    !--                      values for the surface type are already set.
                                    if (.not. btest(topo_tmp(k, j_f, i_f), 12)) then
                                        topo_tmp(k, j_f, i_f) = ibset(topo_tmp(k, j_f, i_f), 11)
                                    end if

                                end if
                            end do
                        end if
                    end if
                end do
            end do
        end do
        deallocate (filter % i)
        deallocate (filter % j)
        deallocate (filter % ij_id)
        !-- Count the total number of filled cavities, required for informative message.
        num_cavity = sum(cav_filled)
        ! #if defined( __parallel )
        !         call MPI_ALLREDUCE(MPI_IN_PLACE, num_cavity, 1, MPI_INTEGER, MPI_SUM, comm2d, ierr)
        ! #endif
        !-- Create an informative message if any narrow cavities were filled.
        if (num_cavity > 0) then
            write (*, '(I6,A,I2,A)') num_cavity, ' narrow cavities of less than ', filter % num_thresh, &
                ' horizontal grid points were filled during initialization'
            ! call message('topography_mod', 'PAC0345', 0, 0, 0, 6, 0)
        end if

        topo(:, nys:nyn, nxl:nxr) = topo_tmp(:, nys:nyn, nxl:nxr)
        !-- Finally, exchange topo array again and if necessary set Neumann boundary condition in case of
        !-- non-cyclic lateral boundaries.
        call exchange_horiz_int(topo, nys, nyn, nxl, nxr, nzt, nbgp)
        call topography_set_non_cyc_bc(topo)
        !-- Exchange building ID and type. Note, building_type is an one-byte variable.
        if (building_id_f % from_file) then
            call exchange_horiz_2d_int(building_id_f % var, nys, nyn, nxl, nxr, nbgp)
            call set_lateral_neumann_bc(building_id_f % var)
        end if
        if (building_type_f % from_file) then
            call exchange_horiz_2d_byte(building_type_f % var, nys, nyn, nxl, nxr, nbgp)
            call set_lateral_neumann_bc(building_type_f % var)
        end if

    contains

        !--------------------------------------------------------------------------------------------------!
        ! Description:
        ! -------------------------------------------------------------------------------------------------!
        !> Core of the cavity filter algorithm. This routine checks if horizontally surrounding grid points
        !> of (j,i) belong to fluid or topography. If they belong to fluid, recursively extend the search
        !> until a threshold value is reached (meaning that the cavity is wide enough), or until all farther
        !> surrounding fluid grid points are counted.
        !--------------------------------------------------------------------------------------------------!
        recursive subroutine trace_surrounding_gridpoints(k, j, i)

            integer(int32) :: i !< input grid index in x-direction
            integer(int32) :: i_trace !< grid index in x-direction for next trace grid cell
            integer(int32) :: j !< input grid index in y-direction
            integer(int32) :: j_trace !< grid index in y-direction for next trace grid cell
            integer(int32) :: k !< input grid index in z-direction
            integer(int32) :: n !< loop variable for the 4 grid-line directions

            integer(int32), dimension(4) :: off_x = (/-1, 1, 0, 0/)
            integer(int32), dimension(4) :: off_y = (/0, 0, -1, 1/)

            !--    From the input grid point, start to look around the 4 grid-line surrounding grid points.
            do n = 1, size(off_x)
                !--       Exit subroutine if more than filter%num_thresh (currently 9) fluid grid points are found.
                if (filter % num_gp >= filter % num_thresh) return

                i_trace = i + off_x(n)
                j_trace = j + off_y(n)
                !--       Skip trace coordinates that are out of the search- or subdomain boundaries.
                if (j_trace < bs .or. j_trace > bn .or. i_trace < bl .or. i_trace > br) cycle
                !--       Check if the grid point has been already counted. This is identified by a unique ID each
                !--       (j,i) is given.
                if (.not. any(gridpoint_id(j_trace, i_trace) == filter % ij_id)) then
                    !--          If the trace grid point (j_trace,i_trace) is a fluid grid point, add its
                    !--          coordinates to the already existing list of coordinates.
                    if (.not. btest(topo_tmp(k, j_trace, i_trace), 0)) then
                        !--             Increment the fluid grid point counter.
                        filter % num_gp = filter % num_gp + 1
                        !--             Add new index pair and a corresponding unique ID.
                        filter % i(filter % num_gp) = i_trace
                        filter % j(filter % num_gp) = j_trace
                        filter % ij_id(filter % num_gp) = gridpoint_id(j_trace, i_trace)
                        !--             Based on the updated (j_trace,i_trace) coordinate, further search for
                        !--             fluid grid points in its vicinity.
                        if (filter % num_gp < filter % num_thresh) &
                            call trace_surrounding_gridpoints(k, j_trace, i_trace)
                    end if
                end if
            end do

        end subroutine trace_surrounding_gridpoints

    end subroutine filter_topography

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! -------------------------------------------------------------------------------------------------!
    !> Set static topography flags
    !--------------------------------------------------------------------------------------------------!
    subroutine topography_set_flags

        integer(int32) :: i !< index variable along x
        integer(int32) :: ibit !< integer bit position of topgraphy masking array
        integer(int32) :: j !< index variable along y
        integer(int32) :: k !< index variable along z

        allocate (topo_flags(nzb:nzt + 1, nysg:nyng, nxlg:nxrg))
        topo_flags = 0
        !-- Set-up topography flags. First, set flags only for s, u, v and w-grid.
        !-- Further special flags will be set in following loops. Note, topography is defined when the
        !-- temporary array "topo", bit 0 is one, else it is atmosphere.
        do i = nxl, nxr
            do j = nys, nyn
                do k = nzb, nzt + 1
                    !--          scalar grid
                    if (.not. btest(topo(k, j, i), 0)) topo_flags(k, j, i) = ibset(topo_flags(k, j, i), 0)
                end do
                !--       Cartesian topography.
                ! if (.not. cut_cell_topography) then
                if (.true.) then
                    do k = nzb, nzt + 1
                        !--             u grid.
                        if (.not. btest(topo(k, j, i), 0) .and. .not. btest(topo(k, j, i - 1), 0)) &
                            topo_flags(k, j, i) = ibset(topo_flags(k, j, i), 1)
                        !--             v grid.
                        if (.not. btest(topo(k, j, i), 0) .and. .not. btest(topo(k, j - 1, i), 0)) &
                            topo_flags(k, j, i) = ibset(topo_flags(k, j, i), 2)
                    end do

                    do k = nzb, nzt
                        !--             w grid.
                        if (.not. btest(topo(k, j, i), 0) .and. .not. btest(topo(k + 1, j, i), 0)) &
                            topo_flags(k, j, i) = ibset(topo_flags(k, j, i), 3)
                    end do
                !--       Slanted surface approach´.
                else
                    do k = nzb, nzt + 1
                        !--             u grid.
                        if (.not. btest(topo(k, j, i), 1)) &
                            topo_flags(k, j, i) = ibset(topo_flags(k, j, i), 1)
                        !--             v grid.
                        if (.not. btest(topo(k, j, i), 2)) &
                            topo_flags(k, j, i) = ibset(topo_flags(k, j, i), 2)

                    end do

                    do k = nzb, nzt
                        !--             w grid.
                        if (.not. btest(topo(k, j, i), 3)) &
                            topo_flags(k, j, i) = ibset(topo_flags(k, j, i), 3)
                    end do
                end if

                if (topography /= 'closed_channel') then
                    topo_flags(nzt + 1, j, i) = ibset(topo_flags(nzt + 1, j, i), 3)
                end if

            end do
        end do

        call exchange_horiz_int(topo_flags, nys, nyn, nxl, nxr, nzt, nbgp)
        !-- Set outer array for scalars to mask near-surface grid points. Note, on basis of flag 24 futher
        !-- flags will be derived which are used to control production of subgrid TKE production near walls.
        do i = nxl, nxr
            do j = nys, nyn
                do k = nzb, nzt + 1
                    if (btest(topo_flags(k, j - 1, i), 0) .and. &
                        btest(topo_flags(k, j + 1, i), 0) .and. &
                        btest(topo_flags(k, j, i - 1), 0) .and. &
                        btest(topo_flags(k, j, i + 1), 0) .and. &
                        btest(topo_flags(k, j - 1, i - 1), 0) .and. &
                        btest(topo_flags(k, j + 1, i - 1), 0) .and. &
                        btest(topo_flags(k, j - 1, i + 1), 0) .and. &
                        btest(topo_flags(k, j + 1, i + 1), 0)) &
                        then
                        topo_flags(k, j, i) = ibset(topo_flags(k, j, i), 24)
                    end if
                end do
            end do
        end do
        !-- Set further special flags.
        do i = nxl, nxr
            do j = nys, nyn
                do k = nzb, nzt + 1
                    !--          Scalar grid, former nzb_diff_s_inner.
                    !--          Note, use this flag also to mask topography in diffusion_u and diffusion_v along the
                    !--          vertical direction. In case of use_surface_fluxes, fluxes are calculated via MOST,
                    !--          else, simple gradient approach is applied. Please note, in case of u- and v-diffuison,
                    !--          a small error is made at edges (on the east side for u, at the north side for v), since
                    !--          topography on scalar grid point is used instead of topography on u/v-grid. As number of
                    !--          topography grid points on uv-grid is different than s-grid, different number of surface
                    !--          elements would be required. In order to avoid this, treat edges (u(k,j,i+1)) simply by
                    !--          a gradient approach, i.e. these points are not masked within diffusion_u. Tests had
                    !--          shown that the effect on the flow is negligible.
                    if (constant_flux_layer .or. use_surface_fluxes) then
                        if (btest(topo_flags(k, j, i), 0)) &
                            topo_flags(k, j, i) = ibset(topo_flags(k, j, i), 8)
                    else
                        topo_flags(k, j, i) = ibset(topo_flags(k, j, i), 8)
                    end if

                end do
                !--       Special flag to control vertical diffusion at model top - former nzt_diff.
                topo_flags(:, j, i) = ibset(topo_flags(:, j, i), 9)
                if (use_top_fluxes) topo_flags(nzt + 1, j, i) = ibclr(topo_flags(nzt + 1, j, i), 9)

                do k = nzb + 1, nzt
                    !--          Special flag on u grid, former nzb_u_inner + 1, required for disturb_field and
                    !--          initialization. Do not disturb directly at topography, as well as initialize u with
                    !--          zero one grid point outside of topography.
                    if (btest(topo_flags(k - 1, j, i), 1) .and. &
                        btest(topo_flags(k, j, i), 1) .and. &
                        btest(topo_flags(k + 1, j, i), 1)) &
                        topo_flags(k, j, i) = ibset(topo_flags(k, j, i), 20)
                    !--          Special flag on v grid, former nzb_v_inner + 1, required for disturb_field and
                    !--          initialization. Do not disturb directly at topography, as well as initialize v with
                    !--          zero one grid point outside of topography.
                    if (btest(topo_flags(k - 1, j, i), 2) .and. &
                        btest(topo_flags(k, j, i), 2) .and. &
                        btest(topo_flags(k + 1, j, i), 2)) &
                        topo_flags(k, j, i) = ibset(topo_flags(k, j, i), 21)
                    !--          Special flag on scalar grid, former nzb_s_inner+1. Used for lpm_sgs_tke.
                    if (btest(topo_flags(k, j, i), 0) .and. &
                        btest(topo_flags(k - 1, j, i), 0) .and. &
                        btest(topo_flags(k + 1, j, i), 0)) &
                        topo_flags(k, j, i) = ibset(topo_flags(k, j, i), 25)
                    !--          Special flag on scalar grid, nzb_diff_s_outer - 1, required in in production_e.
                    if (constant_flux_layer .or. use_surface_fluxes) then
                        if (btest(topo_flags(k, j, i), 24) .and. &
                            btest(topo_flags(k - 1, j, i), 24) .and. &
                            btest(topo_flags(k + 1, j, i), 0)) &
                            topo_flags(k, j, i) = ibset(topo_flags(k, j, i), 29)
                    else
                        if (btest(topo_flags(k, j, i), 0)) &
                            topo_flags(k, j, i) = ibset(topo_flags(k, j, i), 29)
                    end if
                    !--          Special flag on scalar grid, nzb_diff_s_outer - 1, required in
                    !--          in production_e.
                    if (constant_flux_layer .or. use_surface_fluxes) then
                        if (btest(topo_flags(k, j, i), 0) .and. &
                            btest(topo_flags(k - 1, j, i), 0) .and. &
                            btest(topo_flags(k + 1, j, i), 0)) &
                            topo_flags(k, j, i) = ibset(topo_flags(k, j, i), 30)
                    else
                        if (btest(topo_flags(k, j, i), 0)) &
                            topo_flags(k, j, i) = ibset(topo_flags(k, j, i), 30)
                    end if
                end do
                !--       Flags indicating downward facing walls.
                do k = nzb + 1, nzt + 1
                    !--          Scalar grid.
                    if (btest(topo_flags(k - 1, j, i), 0) .and. &
                        .not. btest(topo_flags(k, j, i), 0)) &
                        topo_flags(k, j, i) = ibset(topo_flags(k, j, i), 13)
                    !--          Downward facing wall on u grid.
                    if (btest(topo_flags(k - 1, j, i), 1) .and. &
                        .not. btest(topo_flags(k, j, i), 1)) &
                        topo_flags(k, j, i) = ibset(topo_flags(k, j, i), 15)
                    !--          Downward facing wall on v grid.
                    if (btest(topo_flags(k - 1, j, i), 2) .and. &
                        .not. btest(topo_flags(k, j, i), 2)) &
                        topo_flags(k, j, i) = ibset(topo_flags(k, j, i), 17)
                    !--          Downward facing wall on w grid.
                    if (btest(topo_flags(k - 1, j, i), 3) .and. &
                        .not. btest(topo_flags(k, j, i), 3)) &
                        topo_flags(k, j, i) = ibset(topo_flags(k, j, i), 19)
                end do
                !--       Flags indicating upward facing walls.
                do k = nzb, nzt
                    !--          Upward facing wall on scalar grid.
                    if (.not. btest(topo_flags(k, j, i), 0) .and. &
                        btest(topo_flags(k + 1, j, i), 0)) &
                        topo_flags(k, j, i) = ibset(topo_flags(k, j, i), 12)
                    !--          Upward facing wall on u grid.
                    if (.not. btest(topo_flags(k, j, i), 1) .and. &
                        btest(topo_flags(k + 1, j, i), 1)) &
                        topo_flags(k, j, i) = ibset(topo_flags(k, j, i), 14)

                    !--          Upward facing wall on v grid.
                    if (.not. btest(topo_flags(k, j, i), 2) .and. &
                        btest(topo_flags(k + 1, j, i), 2)) &
                        topo_flags(k, j, i) = ibset(topo_flags(k, j, i), 16)

                    !--          Upward facing wall on w grid.
                    if (.not. btest(topo_flags(k, j, i), 3) .and. &
                        btest(topo_flags(k + 1, j, i), 3)) &
                        topo_flags(k, j, i) = ibset(topo_flags(k, j, i), 18)
                    !--          Special flag on scalar grid, former nzb_s_inner.
                    if (btest(topo_flags(k, j, i), 0) .or. &
                        btest(topo_flags(k, j, i), 12) .or. &
                        btest(topo_flags(k, j, i), 13)) &
                        topo_flags(k, j, i) = ibset(topo_flags(k, j, i), 22)
                    !--          Special flag on scalar grid, nzb_diff_s_inner - 1, required for flow_statistics.
                    if (constant_flux_layer .or. use_surface_fluxes) then
                        if (btest(topo_flags(k, j, i), 0) .and. &
                            btest(topo_flags(k + 1, j, i), 0)) &
                            topo_flags(k, j, i) = ibset(topo_flags(k, j, i), 23)
                    else
                        if (btest(topo_flags(k, j, i), 22)) &
                            topo_flags(k, j, i) = ibset(topo_flags(k, j, i), 23)
                    end if

                end do
                topo_flags(nzt + 1, j, i) = ibset(topo_flags(nzt + 1, j, i), 22)
                topo_flags(nzt + 1, j, i) = ibset(topo_flags(nzt + 1, j, i), 23)
                !--       Set flags indicating that topography is close by in horizontal direction, i.e. flags that
                !--       infold the topography. These will be used to set advection flags for passive scalars,
                !--       where due to large gradients near buildings stationary numerical oscillations can produce
                !--       unrealistically high concentrations. This is only necessary if WS-scheme is applied for
                !--       scalar advection. Note, these flags will be only used for passive scalars such as chemical
                !--       species or aerosols.
                if (scalar_advec == 'ws-scheme') then
                    do k = nzb, nzt
                        if (btest(topo_flags(k, j, i), 0) .and. ( &
                            any(.not. btest(topo_flags(k, j - 3:j + 3, i - 1), 0)) .or. &
                            any(.not. btest(topo_flags(k, j - 3:j + 3, i - 2), 0)) .or. &
                            any(.not. btest(topo_flags(k, j - 3:j + 3, i - 3), 0)) .or. &
                            any(.not. btest(topo_flags(k, j - 3:j + 3, i + 1), 0)) .or. &
                            any(.not. btest(topo_flags(k, j - 3:j + 3, i + 2), 0)) .or. &
                            any(.not. btest(topo_flags(k, j - 3:j + 3, i + 3), 0)) .or. &
                            any(.not. btest(topo_flags(k, j - 1, i - 3:i + 3), 0)) .or. &
                            any(.not. btest(topo_flags(k, j - 2, i - 3:i + 3), 0)) .or. &
                            any(.not. btest(topo_flags(k, j - 3, i - 3:i + 3), 0)) .or. &
                            any(.not. btest(topo_flags(k, j + 1, i - 3:i + 3), 0)) .or. &
                            any(.not. btest(topo_flags(k, j + 2, i - 3:i + 3), 0)) .or. &
                            any(.not. btest(topo_flags(k, j + 3, i - 3:i + 3), 0)) &
                            )) &
                            then
                            topo_flags(k, j, i) = ibset(topo_flags(k, j, i), 31)
                        end if
                    end do
                end if
            end do
        end do
        !-- Set further flags to indicate topography / atmosphere at locations where advective/diffusive
        !-- fluxes are defined. Fluxes of scalars are defined on the u-, v- or w-grid points, so no
        !-- additional information is required. Same for fluxes of u in x-direction, v in y-direction and
        !-- w in z-direction, which are all defined on the scalar grid. Hence, additional information is
        !-- only required for u in y- and z-direction, for v in x- and z-direction, as well as for w in
        !-- x- and y-direction. The information has been already computed in cut_cell_topography_mod.
        !-- bit 7: u in y-direction, bit 10: u in z-direction, bit 11: v in x-direction,
        !-- bit 26: v in z-direction, bit 27: w in x-direction, bit 28: w in y-direction.
        !-- In case of non-cut-cell topography, the corresponding flags are actually not required and are
        !-- set to true everywhere.
        topo_flags = ibset(topo_flags, 7)
        topo_flags = ibset(topo_flags, 10)
        topo_flags = ibset(topo_flags, 11)
        topo_flags = ibset(topo_flags, 26)
        topo_flags = ibset(topo_flags, 27)
        topo_flags = ibset(topo_flags, 28)
        ! if (cut_cell_topography) then
        !     !--    If no advective flux is defined, indicated by BTEST( topo, nr ), clear the flag.
        !     topo_flags = merge(ibclr(topo_flags, 7), ibset(topo_flags, 7), btest(topo, 7))
        !     topo_flags = merge(ibclr(topo_flags, 10), ibset(topo_flags, 10), btest(topo, 9))
        !     topo_flags = merge(ibclr(topo_flags, 11), ibset(topo_flags, 11), btest(topo, 10))
        !     topo_flags = merge(ibclr(topo_flags, 26), ibset(topo_flags, 26), btest(topo, 26))
        !     topo_flags = merge(ibclr(topo_flags, 27), ibset(topo_flags, 27), btest(topo, 27))
        !     topo_flags = merge(ibclr(topo_flags, 28), ibset(topo_flags, 28), btest(topo, 28))
        ! end if

        !-- Finally, set identification flags indicating natural terrain or buildings.
        !-- Natural terrain grid points. Information on the type of the surface is stored in bit 1 of
        !-- 3D Integer array topo. However, this bit is only set when topography is read from file. In order
        !-- to run the land-surface model also without topography information, set bit 1 explicitely in this
        !-- case.
        !-- Natural terrain grid points.
        !-- If no topography is initialized, the land-surface is at k = nzb.
        if (trim(topography) /= 'read_from_file') then
            topo_flags(nzb, :, :) = ibset(topo_flags(nzb, :, :), 5)
        else
            do i = nxl, nxr
                do j = nys, nyn
                    do k = nzb, nzt + 1
                        !--             Natural terrain grid point.
                        if (btest(topo(k, j, i), 11)) &
                            topo_flags(k, j, i) = ibset(topo_flags(k, j, i), 5)
                    end do
                end do
            end do
        end if
        !-- Building grid points.
        do i = nxl, nxr
            do j = nys, nyn
                do k = nzb, nzt + 1
                    if (btest(topo(k, j, i), 12)) &
                        topo_flags(k, j, i) = ibset(topo_flags(k, j, i), 6)
                end do
            end do
        end do
        !-- Set flag 4, indicating new topography grid points due to filtering.
        do i = nxl, nxr
            do j = nys, nyn
                do k = nzb, nzt + 1
                    if (btest(topo(k, j, i), 15)) &
                        topo_flags(k, j, i) = ibset(topo_flags(k, j, i), 4)
                end do
            end do
        end do

        !-- Exchange ghost points for wall flags
        call exchange_horiz_int(topo_flags, nys, nyn, nxl, nxr, nzt, nbgp)
        !-- Set boundary conditions also for flags. Can be interpreted as Neumann boundary conditions for
        !-- topography.
        if (.not. bc_ns_cyc) then
            if (nys == 0) then
                do i = 1, nbgp
                    topo_flags(:, nys - i, :) = topo_flags(:, nys, :)
                end do
            end if
            if (nyn == ny) then
                do i = 1, nbgp
                    topo_flags(:, nyn + i, :) = topo_flags(:, nyn, :)
                end do
            end if
        end if
        if (.not. bc_lr_cyc) then
            if (nxl == 0) then
                do i = 1, nbgp
                    topo_flags(:, :, nxl - i) = topo_flags(:, :, nxl)
                end do
            end if
            if (nxr == nx) then
                do i = 1, nbgp
                    topo_flags(:, :, nxr + i) = topo_flags(:, :, nxr)
                end do
            end if
        end if
        !-- Pre-calculate topography top indices.
        allocate (topo_top_ind(nysg:nyng, nxlg:nxrg, 0:6))
        !-- Uppermost topography index on scalar grid.
        ibit = 12
        topo_top_ind(:, :, 0) = maxloc(merge(1, 0, btest(topo_flags(:, :, :), ibit)), DIM=1) - 1
        !-- Uppermost topography index on u grid.
        ibit = 14
        topo_top_ind(:, :, 1) = maxloc(merge(1, 0, btest(topo_flags(:, :, :), ibit)), DIM=1) - 1
        !-- Uppermost topography index on v grid.
        ibit = 16
        topo_top_ind(:, :, 2) = maxloc(merge(1, 0, btest(topo_flags(:, :, :), ibit)), DIM=1) - 1
        !-- Uppermost topography index on w grid.
        ibit = 18
        topo_top_ind(:, :, 3) = maxloc(merge(1, 0, btest(topo_flags(:, :, :), ibit)), DIM=1) - 1
        !-- Uppermost topography index on scalar outer grid.
        ibit = 24
        topo_top_ind(:, :, 4) = maxloc(merge(1, 0, btest(topo_flags(:, :, :), ibit)), DIM=1) - 1
        !-- Uppermost topography index including full-3D geometry.
        ibit = 12
        do k = nzb, nzt + 1
            where (btest(topo_flags(k, :, :), ibit)) topo_top_ind(:, :, 5) = k
        end do
        !-- Pre-calculate top index of uppermost terrain grid point. This is used for terrain-following
        !-- masked data output.
        topo_top_ind(:, :, 6) = minloc(merge(1, 0, btest(topo_flags, 5)), DIM=1) - 1

    end subroutine topography_set_flags

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! -------------------------------------------------------------------------------------------------!
    !> Set temporary topography flags and reference buildings on top of underlying orography.
    !--------------------------------------------------------------------------------------------------!
    subroutine process_topography

        integer(int32) :: dim_builds !< total number of buildings within the model domain
        integer(int32) :: i !< running index along x-direction
        integer(int32) :: j !< running index along y-direction
        integer(int32) :: k !< running index along z-direction with respect to numeric grid
        integer(int32) :: k2 !< running index along z-direction with respect to netcdf grid
        integer(int32) :: nr !< index variable indication maximum terrain height for respective building ID
        integer(int32) :: topo_top_index !< orography top index, used to map 3D buildings onto terrain
        integer(iwp) :: n

        real(real32) :: ocean_offset !< offset to consider inverse vertical coordinate at topography
        real(wp) :: dz_stretched !< stretched vertical grid spacing
        real(wp) :: dz_level_end !< distance between calculated height level for u/v-grid and user-specified end level for stretching
        !< definition

        real(real32), dimension(:), allocatable :: oro_max !< maximum terrain height occupied by an building with certain id
        real(real32), dimension(:), allocatable :: oro_max_l !< maximum terrain height occupied by an building with certain id,
        !< on local subdomain

        !-- In the following, buildings and orography are further preprocessed before they are mapped on the
        !-- LES grid.
        !-- Buildings are mapped on top of the orography by maintaining the roof shape of the building. This
        !-- can be achieved by referencing building on top of the maximum terrain height within the area
        !-- occupied by the respective building. As buildings and terrain height are defined only locally,
        !-- parallelization of this referencing is required (a building can be distributed between different
        !-- cores).
        !-- In a first step, determine the number of buildings with different building id on each PE. In a
        !-- next step, all building ids are gathered into one array which is present to all PEs. For each
        !-- building ID, the maximum terrain height occupied by the respective building is computed and
        !-- distributed to each PE.
        !-- Finally, for each building id and its respective reference orography, builidings are mapped on
        !-- top.
        !--
        !-- First, set topography flags. Bit 1 indicates orography, bit 2 buildings.
        !-- Grid point nzb is topography on the staggered grid, but
        !-- in case of vertically shifted nests, the boundary vertical coordinate does not have its standard
        !-- location 0.0, meaning that the respective boundary is "open".

        ! print *, "Break point."
        topo = ibclr(topo, 0)

        ! print *, "Break point 2."

        ! if (nest_shift_z == 0.0_real32) topo(nzb, :, :) = ibset(topo(nzb, :, :), 0)
        topo(nzb, :, :) = ibset(topo(nzb, :, :), 0)

        ! print *, "Break point 3."
        !-- In order to map topography on PALM grid also in case of ocean simulations, pre-calculate an
        !-- offset value.

        if (.not. allocated(zu)) then
            allocate (zu(nzb:nzt + 1))
            zu(0) = 0.0_wp
            zu(1) = zu(0) + dz(1) * 0.5_wp

            !--    Determine u and v height levels considering the possibility of grid stretching in several
            !--    heights.
            n = 1
            dz_stretch_level_start_index = nzt + 1
            dz_stretch_level_end_index = nzt + 1
            dz_stretched = dz(1)

            !--    The default value of dz_stretch_level_start is negative, thus the first condition is true
            !--    even if no stretching shall be applied. Hence, the second condition is also necessary.
            do k = 2, nzt + 1 - symmetry_flag
                if (dz_stretch_level_start(n) <= zu(k - 1) .and. &
                    dz_stretch_level_start(n) /= -9999999.9_wp) then
                    dz_stretched = dz_stretched * dz_stretch_factor_array(n)

                    if (dz(n) > dz(n + 1)) then
                        dz_stretched = max(dz_stretched, dz(n + 1)) !Restrict dz_stretched to the user-specified (higher) dz
                    else
                        dz_stretched = min(dz_stretched, dz(n + 1)) !Restrict dz_stretched to the user-specified (lower) dz
                    end if

                    if (dz_stretch_level_start_index(n) == nzt + 1) dz_stretch_level_start_index(n) = k - 1

                end if

                zu(k) = zu(k - 1) + dz_stretched

                !--       Make sure that the stretching ends exactly at dz_stretch_level_end
                dz_level_end = abs(zu(k) - dz_stretch_level_end(n))

                if (dz_level_end < dz(n + 1) / 3.0) then
                    zu(k) = dz_stretch_level_end(n)
                    dz_stretched = dz(n + 1)
                    dz_stretch_level_end_index(n) = k
                    n = n + 1
                end if
            end do
            print *, "zu initialized."
        end if
        if (.not. allocated(zw)) then
            allocate (zw(nzb:nzt + 1))
            zw(0) = 0.0_wp
                do k = 1, nzt
                zw(k) = (zu(k) + zu(k + 1)) * 0.5_wp
            end do
            print *, "zw initialized."
        end if
        if (.not. allocated(dzu) .or. .not. allocated(dzw)) then
            allocate (dzu(1:nzt + 1))
            allocate (dzw(1:nzt + 1))
            allocate (ddzu(1:nzt + 1))
            allocate (ddzw(1:nzt + 1))

            !-- Compute grid lengths.
            do k = 1, nzt + 1
                dzu(k) = zu(k) - zu(k - 1)
                ddzu(k) = 1.0_wp / dzu(k)
                dzw(k) = zw(k) - zw(k - 1)
                ddzw(k) = 1.0_wp / dzw(k)
            end do
            print *, "dzu and dzw initialized."
        end if

        ocean_offset = merge(zw(0), 0.0_real32, ocean_mode)

        ! #if defined( __debug )
        !  print *, "Value of 'ocean_offset' is: " // ocean_offset
        write (*, '(A,G0)') "Value of 'ocean_offset' is: ", ocean_offset
        ! ! call message("process_topography", "-1", 0, 3, myid, 6, 1)
        ! #endif
        !-- Reference buildings on top of orography. This is not necessary if topography is read from ASCII
        !-- file as no distinction between buildings and terrain height can be made. Moreover, this is also
        !-- not necessary if urban-surface and land-surface model are used at the same time.
        if (input_pids_static) then

            if (buildings_f % from_file) then
                !--       Determine maximumum terrain height occupied by the respective building and temporalily
                !--       store on oro_max. Before, check whether any buildings are defined within the domain.
                if (allocated(build_ids_unique)) then
                    dim_builds = size(build_ids_unique)
                else
                    dim_builds = 0
                end if

                allocate (oro_max_l(1:dim_builds))
                allocate (oro_max(1:dim_builds))
                oro_max_l = 0.0_real32

                do nr = 1, dim_builds
                    oro_max_l(nr) = maxval(merge(terrain_height_f % var(nys:nyn, nxl:nxr), 0.0_real32, &
                                                 building_id_f % var(nys:nyn, nxl:nxr) == &
                                                 build_ids_unique(nr)))
                end do

                ! #if defined( __parallel )
                !                 if (dim_builds >= 1) then
                !                     call MPI_ALLREDUCE(oro_max_l, oro_max, size(oro_max), MPI_REAL, MPI_MAX, comm2d, &
                !                                        ierr)
                !                 end if
                ! #else
                oro_max = oro_max_l
                ! #endif
                !--       Finally, determine discrete grid height of maximum orography occupied by a building. Use
                !--       all-or-nothing approach, i.e. if terrain exceeds the scalar level the grid box is fully
                !--       terrain and the maximum terrain is set to the zw level.
                oro_max_l = 0.0
                do nr = 1, dim_builds
                    do k = nzb, nzt
                        if (zu(k) - ocean_offset <= oro_max(nr)) oro_max_l(nr) = zw(k) - ocean_offset
                    end do
                    oro_max(nr) = oro_max_l(nr)
                end do
            end if
            !--    Allocate array for storing terrain height under buildings.
            if (buildings_f % from_file) then
                allocate (buildings_f % oro_max(nysg:nyng, nxlg:nxrg))
                buildings_f % oro_max = buildings_f % fill1
            end if
            !--    Map orography as well as buildings onto grid.
            do i = nxl, nxr
                do j = nys, nyn
                    topo_top_index = 0
                    !--          Obtain index in global building_id array.
                    if (buildings_f % from_file) then
                        if (building_id_f % var(j, i) /= building_id_f % fill) then
                            !--                Determine index where maximum terrain height occupied by the respective building
                            !--                height is stored.
                            nr = minloc(abs(build_ids_unique - building_id_f % var(j, i)), DIM=1)
                            !--                Save grid-indexed oro_max.
                            buildings_f % oro_max(j, i) = oro_max(nr)
                        end if
                    end if
                    do k = nzb, nzt
                        !--             In a first step, if grid point is below or equal the given terrain height, grid
                        !--             point is flagged to be of type natural.
                        !--             Please note, in case there is also a building which is lower than the vertical grid
                        !--             spacing, initialization of surface attributes will not be correct as given surface
                        !--             information will not be in accordance to the classified grid points.
                        !--             Hence, in this case, also a building flag.
                        if (zu(k) - ocean_offset <= terrain_height_f % var(j, i)) then
                            topo(k, j, i) = ibset(topo(k, j, i), 0)
                            topo(k, j, i) = ibset(topo(k, j, i), 11)
                            topo_top_index = k ! topo_top_index + 1
                        end if
                        !--             Set building grid points. Here, only consider 2D buildings.
                        !--             3D buildings require separate treatment.  Note, in case of surface-mounted
                        !--             buildings slanted terrain is not considered.
                        if (buildings_f % from_file .and. buildings_f % lod == 1) then
                            !--                Fill-up the terrain to the level of maximum orography within the building-covered
                            !--                area.
                            if (building_id_f % var(j, i) /= building_id_f % fill) then
                                !--                   Note, oro_max is always on zw level.
                                if (zu(k) - ocean_offset < oro_max(nr)) then
                                    topo(k, j, i) = ibset(topo(k, j, i), 0)
                                    topo(k, j, i) = ibset(topo(k, j, i), 11)
                                elseif (zu(k) - ocean_offset <= oro_max(nr) + buildings_f % var_2d(j, i)) then
                                    topo(k, j, i) = ibset(topo(k, j, i), 0)
                                    topo(k, j, i) = ibset(topo(k, j, i), 12)
                                end if
                            end if
                        end if
                    end do
                    !--          Special treatment for non grid-resolved buildings. This case, the uppermost terrain
                    !--          grid point is flagged as building as well, even though no building exists at all.
                    !--          However, the surface element will be identified as urban-surface and the input data
                    !--          provided by the drivers is consistent to the surface classification. Else, all non
                    !--          grid-resolved buildings would vanish and identified as terrain grid points, which,
                    !--          however, won't be consistent with the input data.
                    if (buildings_f % from_file .and. buildings_f % lod == 1) then
                        if (building_id_f % var(j, i) /= building_id_f % fill) then
                            do k = nzb, nzt
                                if (zw(k) - ocean_offset == oro_max(nr)) then
                                    if (buildings_f % var_2d(j, i) <= zu(k + 1) - zw(k)) then
                                        topo(k, j, i) = ibset(topo(k, j, i), 12)
                                    end if
                                end if
                            end do
                        end if
                    end if
                    !--          Map 3D buildings onto terrain height.
                    !--          In case of any slopes, map building on top of maximum terrain height covered by the
                    !--          building. In other words, extend building down to the respective local terrain-surface
                    !--          height.
                    if (buildings_f % from_file .and. buildings_f % lod == 2) then
                        if (building_id_f % var(j, i) /= building_id_f % fill) then
                            !--                Extend building down to the terrain surface, i.e. fill-up surface irregularities
                            !--                below a building. Note, oro_max is already a discrete height according to the
                            !--                all-or-nothing approach, i.e. grid box is either topography or atmosphere,
                            !--                terrain top is defined at upper bound of the grid box.
                            !--                Hence, check for zw in this case.
                            !--                Note, do this only for buildings which are surface mounted, i.e. building types
                            !--                1-6. Below bridges, which are represented exclusively by building type 7, terrain
                            !--                shape should be maintained.
                            if (building_type_f % from_file) then
                                if (building_type_f % var(j, i) /= 7) then
                                    do k = topo_top_index + 1, nzt + 1
                                        if (zu(k) - ocean_offset <= oro_max(nr)) then
                                            topo(k, j, i) = ibset(topo(k, j, i), 0)
                                            topo(k, j, i) = ibset(topo(k, j, i), 11)
                                        end if
                                    end do
                                    !--                      After surface irregularities are smoothen, determine lower start index
                                    !--                      where building starts.
                                    do k = nzb, nzt
                                        if (zu(k) - ocean_offset <= oro_max(nr)) topo_top_index = k
                                    end do
                                end if
                            end if
                            !--                Finally, map building on top.
                            k2 = 0
                            do k = topo_top_index, nzt + 1
                                if (k2 <= buildings_f % nz - 1) then
                                    if (buildings_f % var_3d(k2, j, i) == 1) then
                                        topo(k, j, i) = ibset(topo(k, j, i), 0)
                                        topo(k, j, i) = ibset(topo(k, j, i), 12)
                                    end if
                                end if
                                k2 = k2 + 1
                            end do
                        end if
                    end if
                end do
            end do
            !--    Horizontal exchange the oro_max array, which is required to for initialization of
            !--    building-surface properties.
            if (allocated(buildings_f % oro_max)) then
                call exchange_horiz_2d(buildings_f % oro_max)
                call set_lateral_neumann_bc(buildings_f % oro_max)
            end if
            !--    Deallocate temporary arrays required for processing and reading data
            if (allocated(oro_max)) deallocate (oro_max)
            if (allocated(oro_max_l)) deallocate (oro_max_l)
            if (allocated(build_ids_unique)) deallocate (build_ids_unique)
        !-- Topography input via ASCII format.
        else
            ocean_offset = merge(zw(0), 0.0_real32, ocean_mode)
            !--    Initialize topography bit 0 (indicates obstacle) everywhere to zero and clear all grid points
            !--    at nzb, where always a surface is defined.
            !--    Further, set also bit 1 (indicates terrain) at nzb, which is further used for masked data
            !--    output and further processing. Note, in the ASCII case no distinction is made between
            !--    buildings and terrain, so that setting of bit 1 and 2 at the same time has no effect.
            topo = ibclr(topo, 0)
            topo(nzb, :, :) = ibset(topo(nzb, :, :), 0)
            topo(nzb, :, :) = ibset(topo(nzb, :, :), 11)
            do i = nxl, nxr
                do j = nys, nyn
                    do k = nzb, nzt
                        !--             Flag topography for all grid points which are below the local topography height.
                        !--             Note, each topography is flagged as building (bit 12) as well as terrain (bit 11) in
                        !--             order to employ urban-surface as well as land-surface model.
                        if ((zu(k) - ocean_offset) <= buildings_f % var_2d(j, i)) then
                            topo(k, j, i) = ibset(topo(k, j, i), 0)
                            topo(k, j, i) = ibset(topo(k, j, i), 11)
                            topo(k, j, i) = ibset(topo(k, j, i), 12)
                        end if
                    end do
                end do
            end do
        end if
        !-- Exchange ghost points. Further, in case of non-cyclic boundary conditions Neumann BC
        !-- are set for the topography.
        call exchange_horiz_int(topo, nys, nyn, nxl, nxr, nzt, nbgp)
        call topography_set_non_cyc_bc(topo)

    end subroutine process_topography

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! -------------------------------------------------------------------------------------------------!
    !> Setting of further topography inidices on the grid.
    !--------------------------------------------------------------------------------------------------!
    subroutine topography_set_indices

        integer(int32) :: i !< grid index in x-direction
        integer(int32) :: j !< grid index in y-direction
        integer(int32) :: k !< grid index in z-direction
        integer(int32) :: k_top !< topography top index on local PE
        integer(int32) :: nzb_local_max !< vertical grid index of maximum topography height
        integer(int32) :: nzb_local_min !< vertical grid index of minimum topography height

        !-- Determine the maximum level of topography. It is used for steering the degradation of order of
        !-- the applied advection scheme, as well in the lpm.
        k_top = 0
        do i = nxl, nxr
            do j = nys, nyn
                do k = nzb, nzt + 1
                    k_top = max(k_top, merge(k, 0, btest(topo(k, j, i), 0)))
                end do
            end do
        end do
        ! #if defined( __parallel )
        !         call MPI_ALLREDUCE(k_top, nzb_max, 1, MPI_INTEGER, MPI_MAX, comm2d, ierr)
        ! #else
        nzb_max = k_top
        ! #endif
        !-- Increment nzb_max by 1 in order to allow for proper diverengence correction.
        !-- Further, in case topography extents up to the model top, limit to nzt.
        nzb_max = min(nzb_max + 1, nzt)
        !-- Determine minimum index of topography. Usually, this will be nzb. In case there is elevated
        !-- topography, however, the lowest topography will be higher.
        !-- This index is e.g. used to calculate mean first-grid point atmosphere temperature, surface
        !-- pressure and density, etc. .
        topo_min_level = 0
        ! #if defined( __parallel )
        !         call MPI_ALLREDUCE(minval(topo_top_ind(nys:nyn, nxl:nxr, 0)), topo_min_level, 1, MPI_INTEGER, &
        !                            MPI_MIN, comm2d, ierr)
        ! #else
        topo_min_level = minval(topo_top_ind(nys:nyn, nxl:nxr, 0))

        ! #endif
        !-- Check topography for consistency with model domain. Therefore, use maximum and minium
        !-- topography-top indices. Note, minimum topography top index is already calculated.
        if (trim(topography) /= 'flat') then
            ! #if defined( __parallel )
            !             call MPI_ALLREDUCE(maxval(topo_top_ind(nys:nyn, nxl:nxr, 0)), nzb_local_max, 1, &
            !                                MPI_INTEGER, MPI_MAX, comm2d, ierr)
            ! #else
            nzb_local_max = maxval(topo_top_ind(nys:nyn, nxl:nxr, 0))
            ! #endif
            nzb_local_min = topo_min_level
            !--    Consistency checks
            if (nzb_local_min < 0 .or. nzb_local_max > (nz + 1)) then
                write (*, *) 'nzb_local values are outside the model domain', &
                    '&MINVAL( nzb_local ) = ', nzb_local_min, &
                    '&MAXVAL( nzb_local ) = ', nzb_local_max
                ! call message('topography_mod', 'PAC0346', 1, 2, 0, 6, 0)
            end if
        end if
        !-- Allocate and set the arrays containing the topography height (for output reasons only).
        if (trim(topography) /= 'flat') then

            if (nxr == nx .and. nyn /= ny) then
                allocate (zu_s_inner(nxl:nxr + 1, nys:nyn), zw_w_inner(nxl:nxr + 1, nys:nyn))
            elseif (nxr /= nx .and. nyn == ny) then
                allocate (zu_s_inner(nxl:nxr, nys:nyn + 1), zw_w_inner(nxl:nxr, nys:nyn + 1))
            elseif (nxr == nx .and. nyn == ny) then
                allocate (zu_s_inner(nxl:nxr + 1, nys:nyn + 1), zw_w_inner(nxl:nxr + 1, nys:nyn + 1))
            else
                allocate (zu_s_inner(nxl:nxr, nys:nyn), zw_w_inner(nxl:nxr, nys:nyn))
            end if

            zu_s_inner = 0.0_real32
            zw_w_inner = 0.0_real32
            !--    Determine local topography height on scalar and w-grid. Note, setting lateral boundary values
            !--    is not necessary, realized via topo_flags array. Further, please note that loop
            !--    bounds are different from nxl to nxr and nys to nyn on south and right model boundary, hence,
            !--    use intrinsic lbound and ubound functions to infer array bounds.
            do i = lbound(zu_s_inner, 1), ubound(zu_s_inner, 1)
                do j = lbound(zu_s_inner, 2), ubound(zu_s_inner, 2)
                    !--          Topography height on scalar grid. Therefore, determine index of upward-facing surface
                    !--          element on scalar grid.
                    zu_s_inner(i, j) = zu(topo_top_ind(j, i, 0))
                    !--          Topography height on w grid. Therefore, determine index of upward-facing surface
                    !--          element on w grid.
                    zw_w_inner(i, j) = zw(topo_top_ind(j, i, 3))
                end do
            end do
        end if

    end subroutine topography_set_indices

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! -------------------------------------------------------------------------------------------------!
    !> Routine to set boundary conditions which is repeatedly necessary during topography
    !> initialization. Note, this routine is just to avoid lengthy code.
    !--------------------------------------------------------------------------------------------------!
    subroutine topography_set_non_cyc_bc(topo_array)

        integer(int32) :: i !< running index over ghost layers

        integer(int32), dimension(nzb:nzt + 1, nysg:nyng, nxlg:nxrg) :: topo_array !< input array for 3D topography

        !-- Set boundary conditions also for flags. Can be interpreted as Neumann boundary conditions
        !-- for topography.
        if (.not. bc_ns_cyc) then
            if (nys == 0) then
                do i = 1, nbgp
                    topo_array(:, nys - i, :) = topo_array(:, nys, :)
                end do
            end if
            if (nyn == ny) then
                do i = 1, nbgp
                    topo_array(:, nyn + i, :) = topo_array(:, nyn, :)
                end do
            end if
        end if
        if (.not. bc_lr_cyc) then
            if (nxl == 0) then
                do i = 1, nbgp
                    topo_array(:, :, nxl - i) = topo_array(:, :, nxl)
                end do
            end if
            if (nxr == nx) then
                do i = 1, nbgp
                    topo_array(:, :, nxr + i) = topo_array(:, :, nxr)
                end do
            end if
        end if

    end subroutine topography_set_non_cyc_bc

end module topography_mod
