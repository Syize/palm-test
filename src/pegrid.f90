!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Definition of variables which define processor topology and the exchange of ghost point layers.
!> This module must be placed in all routines containing MPI-calls.
!--------------------------------------------------------------------------------------------------!
module pegrid

    use kinds
    ! use, intrinsic :: iso_fortran_env, only: int32

    ! integer, parameter :: iwp = int32

    ! #if defined( __parallel )
    !     use MPI, &
    !         only: MPI_STATUS_SIZE
    !     private MPI_STATUS_SIZE
    ! #endif
    character(LEN=7) :: myid_char = '' !< character string containing processor id number

    integer(iwp) :: comm1dx !< communicator for domain decomposition along x
    integer(iwp) :: comm1dy !< communicator for domain decomposition along y
    integer(iwp) :: comm2d !< standard 2d (xy) communicator used in PALM for the process group the PE belongs
    !< to
    integer(iwp) :: comm_palm !< internal communicator used during the MPI setup at the beginning of a run
    integer(iwp) :: id_outflow = 0 !< myidx of procs at outflow (turbulent outflow method)
    integer(iwp) :: id_outflow_source = 0 !< myidx of procs including ouflow source plane (turbulent outflow method)
    integer(iwp) :: ierr !< standard error parameter in MPI calls
    integer(iwp) :: myid = 0 !< id number of processor element
    integer(iwp) :: myidx = 0 !< id number of processor elements with same position along x-direction
    integer(iwp) :: myidy = 0 !< id number of processor elements with same position along y-direction
    integer(iwp) :: ndim = 2 !< dimension of the virtual PE grid
    integer(iwp) :: ngp_xy !< used in atmosphere/ocean coupling: number of grid points of the subdomain
    integer(iwp) :: ngp_y !< number of subdomain grid points along y including ghost points
    integer(iwp) :: npex = -1 !< number of processor elements in x-direction
    integer(iwp) :: npey = -1 !< number of processor elements in y-direction
    integer(iwp) :: numprocs = 1 !< total number of appointed processor elements
    integer(iwp) :: numprocs_previous_run = -1 !< total number of appointed processor elements in previous run (job chain)
    integer(iwp) :: num_acc_devices = 0 !< number of devices on a node visible to OpenACC
    integer(iwp) :: pleft !< MPI id of left neigbour pe
    integer(iwp) :: pnorth !< MPI id of right neigbour pe
    integer(iwp) :: pright !< MPI id of south neigbour pe
    integer(iwp) :: psouth !< MPI id of north neigbour pe
    integer(iwp) :: req_count = 0 !< MPI return variable - checks if Send-Receive operation is already finished
    integer(iwp) :: sendrecvcount_xy !< number of subdomain gridpoints to be exchanged in direct transpositions
    !< (y --> x, or x --> y) or second (2d) transposition x --> y
    integer(iwp) :: sendrecvcount_yz !< number of subdomain gridpoints to be exchanged in third (2d) transposition
    !< y --> z
    integer(iwp) :: sendrecvcount_zx !< number of subdomain gridpoints to be exchanged in first (2d) transposition
    !< z --> x
    integer(iwp) :: sendrecvcount_zyd
    !< number of subdomain gridpoints to be exchanged in direct transpositions z --> y (used for calculating spectra)
    integer(iwp) :: target_id !< in atmosphere/ocean coupling: id of the ocean/atmosphere counterpart PE with
    !< whom the atmosphere/ocean PE exchanges data
    integer(iwp) :: tasks_per_node = -9999 !< MPI tasks per compute node
    integer(iwp) :: threads_per_task = 1 !< number of OPENMP threads per MPI task
    integer(iwp) :: type_x !< derived MPI datatype for 2-D ghost-point exchange - north / south
    integer(iwp) :: type_xy !< derived MPI datatype for 2-D ghost-point exchange - north / south
    integer(iwp) :: type_y !< derived MPI datatype for 2-D exchange in atmosphere-ocean coupler

    integer(iwp) :: req(100) !< MPI return variable indicating if send-receive operation is finished

    integer(iwp), dimension(:, :), allocatable :: hor_index_bounds !< horizontal index bounds
    integer(iwp), dimension(:, :), allocatable :: hor_index_bounds_previous_run !< horizontal index bounds of previous run

    logical :: collective_wait = .false. !< switch to set an explicit MPI barrier in front of all collective MPI calls
    logical :: non_uniform_subdomain = .false. !< subdomains are non-uniform along x and/or y
    logical :: non_uniform_data_for_transpose = .false. !< data to be transposed is non-uniformly distributed among the cores

    type virtual_pe_grid
        integer(iwp) :: mpi_communicator !< MPI communicator id
        integer(iwp) :: pleft !< MPI id of left neigbour pe
        integer(iwp) :: pright !< MPI id of right neigbour pe
        integer(iwp) :: psouth !< MPI id of south neigbour pe
        integer(iwp) :: pnorth !< MPI id of north neigbour pe
    end type virtual_pe_grid

    type(virtual_pe_grid) :: communicator_configurations(4) !< stores the four possible 2d virtual grids:
    !< cyclic, cyclic along x, cyclic along y, non-cyclic







    ! #if defined( __parallel )
    !     integer(iwp) :: ibuf(12) !< internal buffer for calculating MPI settings
    !     integer(iwp) :: pcoord(2) !< PE coordinates along x and y
    !     integer(iwp) :: status(MPI_STATUS_SIZE) !< MPI status variable used in various MPI calls
    !     integer(iwp) :: type_x_byte !< derived MPI datatype for 2-D 8-bit integer ghost-point exchange - north / south
    !     integer(iwp) :: type_y_byte !< derived MPI datatype for 2-D integer ghost-point exchange - left / right
    !     integer(iwp), dimension(MPI_STATUS_SIZE, 100) :: wait_stat !< MPI status variable used in various MPI calls
    !     integer(iwp) :: type_x_int !< derived MPI datatype for 2-D integer ghost-point exchange - north/south
    !     integer(iwp) :: type_y_int !< derived MPI datatype for 2-D integer ghost-point exchange - left/right
    !     integer(iwp), dimension(:), allocatable :: ngp_xz !< number of ghost points in xz-plane on different multigrid level
    !     integer(iwp), dimension(:), allocatable :: ngp_xz_int !< number of ghost points in xz-plane on different multigrid level
    !     integer(iwp), dimension(:), allocatable :: ngp_yz !< number of ghost points in yz-plane on different multigrid level
    !     integer(iwp), dimension(:), allocatable :: ngp_yz_int !< number of ghost points in yz-plane on different multigrid level
    !     integer(iwp), dimension(:), allocatable :: type_xz !< derived MPI datatype for 3-D integer ghost-point exchange - north /
    !     !< south
    !     integer(iwp), dimension(:), allocatable :: type_xz_int !< derived MPI datatype for 3-D integer ghost-point exchange - north /
    !     !< south
    !     integer(iwp), dimension(:), allocatable :: type_yz !< derived MPI datatype for 3-D integer ghost-point exchange - left /
    !     !< right
    !     integer(iwp), dimension(:), allocatable :: type_yz_int !< derived MPI datatype for 3-D integer ghost-point exchange - left /
    !     !< right
    !     logical :: left_border_pe = .false. !< = .TRUE. if PE is on left border of computational domain
    !     logical :: north_border_pe = .false. !< = .TRUE. if PE is on north border of computational domain
    !     logical :: reorder = .false. !< flag to allow MPI the reorder of ranking (e.g. row-major or column-major)
    !     logical :: right_border_pe = .false. !< = .TRUE. if PE is on right border of computational domain
    !     logical :: south_border_pe = .false. !< = .TRUE. if PE is on south border of computational domain
    !     logical, dimension(2) :: cyclic = (/.true., .true./) !< boundary conditions of the virtual PE grid
    !     logical, dimension(2) :: remain_dims !< internal array used to determine sub-topologies for transpositions
    ! #endif
    save

end module pegrid