module palm
    use topography_mod, only: init_topography
    use indices, only: topo_top_ind, nysg, nyng, nxlg, nxrg
    use arrays_3d, only: zu, zw
    use kinds
    use netcdf_data_input_mod, only: input_file_static

    implicit none

    integer :: i, j, point_index
    ! integer :: topo_top_ind_3dim_start = 0
    ! integer :: topo_top_ind_3dim_end = 6
    real(wp), dimension(:, :), allocatable :: terrain_height

    logical :: read_data_flag = .false.

    private

    public :: extract_scalar_grid_terrain_height, extract_w_grid_terrain_height, terrain_height, say_hello
contains

    subroutine say_hello
        print *, "Hi, palm interface is built successfully."
        
    end subroutine say_hello

    subroutine read_data
        if (read_data_flag) then
            print *, 'No need to call read_data again.'
            return
        end if

        call init_topography

        read_data_flag = .true.

        print *, ''
        print *, 'Topography data readed.'
    end subroutine read_data

    subroutine init_terrain_height_array
        if (allocated(terrain_height)) then
            deallocate(terrain_height)
        end if

        allocate(terrain_height(nysg:nyng, nxlg:nxrg))
    end subroutine init_terrain_height_array

    subroutine extract_scalar_grid_terrain_height (file_path)
        character(len=100), intent(in) :: file_path

        call set_input_file_path (file_path)

        if (.not. read_data_flag) then
            print *, 'Read data'
            call read_data
        end if

        call init_terrain_height_array

        do j = nysg, nyng
            do i = nxlg, nxrg
                point_index = topo_top_ind(j, i, 0)
                terrain_height(j, i) = zu(point_index)
            end do
        end do
    end subroutine extract_scalar_grid_terrain_height

    subroutine extract_w_grid_terrain_height (file_path)
        character(len=100), intent(in) :: file_path

        call set_input_file_path (file_path)

        if (.not. read_data_flag) then
            print *, 'Read data'
            call read_data
        end if

        call init_terrain_height_array

        do j = nysg, nyng
            do i = nxlg, nxrg
                point_index = topo_top_ind(j, i, 0)
                terrain_height(j, i) = zw(point_index)
            end do
        end do
    end subroutine extract_w_grid_terrain_height

    subroutine set_input_file_path (file_path)
        character(len=100), intent(in) :: file_path

        input_file_static = file_path

    end subroutine set_input_file_path
end module palm
