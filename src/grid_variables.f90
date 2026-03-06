module grid_variables

    use kinds

    real(wp) :: ddx !< 1/dx
    real(wp) :: ddx2 !< 1/dx2
    real(wp) :: dx = 1.0_wp !< horizontal grid size (along x-direction)
    real(wp) :: dx2 !< dx*dx
    real(wp) :: ddy !< 1/dy
    real(wp) :: ddy2 !< 1/dy2
    real(wp) :: dy = 1.0_wp !< horizontal grid size (along y-direction)
    real(wp) :: dy2 !< dy*dy

    real(wp), dimension(:, :), allocatable :: zu_s_inner !< height of topography top on scalar grid
    real(wp), dimension(:, :), allocatable :: zw_w_inner !< height of topography top on w grid

    save

end module grid_variables