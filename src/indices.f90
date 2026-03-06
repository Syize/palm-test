module indices

    use kinds

    integer(iwp) :: nbgp = 3 !< number of boundary ghost points
    integer(iwp) :: ngp_sums !< number of vertical profile grid points time number of output profiles - used for allreduce
    !< statements in MPI calls
    integer(iwp) :: ngp_sums_ls !< number of vertical profile grid points time number of large-scale forcing profiles - used for
    !< allreduce statements in MPI calls
    integer(iwp) :: nnx = 226!< number of subdomain grid points in x-direction
    integer(iwp) :: nx = 225 !< nx+1 = total number of grid points in x-direction
    integer(iwp) :: nxl = 0 !< left-most grid index of subdomain (excluding ghost points)
    integer(iwp) :: nxlg = -3 !< left-most grid index of subdomain (including ghost points)
    integer(iwp) :: nxlu !< =nxl+1 (at left domain boundary with inflow from left), else =nxl
    !< (used for u-velocity component)
    integer(iwp) :: nxr = 225 !< right-most grid index of subdomain (excluding ghost points)
    integer(iwp) :: nxrg = 228 !< right-most grid index of subdomain (including ghost points)
    integer(iwp) :: nx_on_file !< nx of previous run in job chain
    integer(iwp) :: nny = 202 !< number of subdomain grid points in y-direction
    integer(iwp) :: ny = 201 !< ny+1 = total number of grid points in y-direction
    integer(iwp) :: nyn = 201 !< north-most grid index of subdomain (excluding ghost points)
    integer(iwp) :: nyng = 204 !< north-most grid index of subdomain (including ghost points)
    integer(iwp) :: nys = 0 !< south-most grid index of subdomain (excluding ghost points)
    integer(iwp) :: nysg = -3 !< south-most grid index of subdomain (including ghost points)
    integer(iwp) :: nysv !< =nys+1 (at south domain boundary with inflow from south), else =nys
    !< (used for v-velocity component)
    integer(iwp) :: ny_on_file !< ny of previous run in job chain
    integer(iwp) :: nnz = 60 !< number of subdomain grid points in z-direction
    integer(iwp) :: nz = 60 !< total number of grid points in z-direction
    integer(iwp) :: nzb = 0 !< bottom grid index of computational domain
    integer(iwp) :: nzb_diff !< will be removed
    integer(iwp) :: nzb_max !< vertical index of topography top
    integer(iwp) :: nzt = 60 !< nzt+1 = top grid index of computational domain
    integer(iwp) :: topo_min_level !< minimum topography-top index (usually equal to nzb)

    integer(iwp), dimension(:), allocatable :: nnx_pe !< grid points along x-direction for every PE
    integer(iwp), dimension(:), allocatable :: nny_pe !< grid points along y-direction for every PE
    integer(iwp), dimension(:), allocatable :: nxl_pe !< lower index bound along x-direction for every PE
    integer(iwp), dimension(:), allocatable :: nxr_pe !< upper index bound along x-direction for every PE
    integer(iwp), dimension(:), allocatable :: nyn_pe !< lower index bound along y-direction for every PE
    integer(iwp), dimension(:), allocatable :: nys_pe !< lower index bound along y-direction for every PE

    integer(iwp), dimension(:), allocatable :: ngp_2dh !< number of grid points of a horizontal cross section through the
    !< total domain
    integer(iwp), dimension(:), allocatable :: ngp_2dh_wgrid !< number of prognostic w-grid points of a horizontal cross section
    !< through the total domain
    integer(idp), dimension(:), allocatable :: ngp_3d !< number of grid points of the total domain
    integer(idp), dimension(:), allocatable :: ngp_3d_inner !< ! need to have 64 bit for grids > 2E9

    integer(iwp), dimension(:, :), allocatable :: ngp_2dh_outer !< number of horizontal grid points which are non-topography and
    !< non-surface-bounded
    integer(iwp), dimension(:, :), allocatable :: ngp_2dh_s_inner !< number of horizontal grid points which are non-topography

    integer(iwp), dimension(:, :, :), allocatable :: advc_flags_m !< flags used to degrade order of advection scheme for
    !< momentum
    integer(iwp), dimension(:, :, :), allocatable :: advc_flags_s !< flags used to degrade order of advection scheme for
    !< scalar quantities
    integer(iwp), dimension(:, :, :), allocatable :: topo_top_ind !< precalculated topography top indices

    integer(iwp), dimension(:, :, :), allocatable :: topo_flags !< flags to mask topography and surface-bounded grid
    !< points

    save

end module indices