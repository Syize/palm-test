module control_parameters

    use kinds

    type file_status
        logical :: opened !< file is currently open
        logical :: opened_before !< file is currently closed, but has been openend before
    end type file_status

    integer(iwp), parameter :: fl_max = 500 !< maximum number of virtual-flight measurements
    integer, parameter :: mask_xyz_dimension = 100 !< limit of mask dimensions (100 points in each direction)
    integer, parameter :: max_masks = 300 !< maximum number of masks
    integer(iwp), parameter :: var_fl_max = 20 !< maximum number of different sampling variables in virtual flight
    !< measurements
    integer(iwp), parameter :: varnamelength = 30 !< length of output variable names

    real(wp), parameter :: output_fill_value = -999999.0_wp

    type(file_status), dimension(200 + 2*max_masks) :: & !< indicates if file is open or if it has been opened before
        openfile = file_status(.false., .false.)

    character(LEN=1) :: timestep_reason = ' ' !< 'A'dvection or 'D'iffusion criterion, written to
    !< RUN_CONTROL file
    character(LEN=5) :: run_zone = ' ' !< time zone of simulation run
    character(LEN=7) :: netcdf_compression_method = 'none' !< netCDF compression mode
    character(LEN=8) :: coupling_char = '' !< appended to filenames in coupled or nested runs
    !< ('_O': ocean PE,
    !< '_NV': vertically nested atmosphere PE, '_N##': PE of
    !< nested domain ##
    character(LEN=8) :: run_time = ' ' !< time of simulation run
    character(LEN=10) :: simulated_time_chr !< simulated time, printed to RUN_CONTROL file
    character(LEN=10) :: run_date = ' ' !< date of simulation run
    character(LEN=11) :: topography_grid_convention = ' ' !< namelist parameter
    character(LEN=26) :: version_string = ' ' !< PALM version
    character(LEN=12) :: user_interface_required_revision = '2.0'
    !< required user-interface revision number (only major number is checked)
    character(LEN=16) :: conserve_volume_flow_mode = 'initial_profiles' !< namelist parameter
    character(LEN=16) :: loop_optimization = 'cache' !< namelist parameter
    character(LEN=16) :: momentum_advec = 'ws-scheme' !< namelist parameter
    character(LEN=16) :: psolver = 'poisfft' !< namelist parameter
    character(LEN=16) :: scalar_advec = 'ws-scheme' !< namelist parameter
    character(LEN=20) :: approximation = 'boussinesq' !< namelist parameter
    character(LEN=20) :: bc_e_b = 'neumann' !< namelist parameter
    character(LEN=20) :: bc_lr = 'cyclic' !< namelist parameter
    character(LEN=20) :: bc_ns = 'cyclic' !< namelist parameter
    character(LEN=20) :: bc_p_b = 'neumann' !< namelist parameter
    character(LEN=20) :: bc_p_t = 'neumann' !< namelist parameter
    character(LEN=20) :: bc_pt_b = 'dirichlet' !< namelist parameter
    character(LEN=20) :: bc_pt_t = 'initial_gradient' !< namelist parameter
    character(LEN=20) :: bc_q_b = 'dirichlet' !< namelist parameter
    character(LEN=20) :: bc_q_t = 'neumann' !< namelist parameter
    character(LEN=20) :: bc_s_b = 'dirichlet' !< namelist parameter
    character(LEN=20) :: bc_s_t = 'initial_gradient' !< namelist parameter
    character(LEN=20) :: bc_uv_b = 'dirichlet' !< namelist parameter
    character(LEN=20) :: bc_uv_t = 'dirichlet' !< namelist parameter
    character(LEN=20) :: dissipation_1d = 'detering' !< namelist parameter
    character(LEN=20) :: fft_method = 'temperton-algorithm' !< namelist parameter
    character(LEN=20) :: mixing_length_1d = 'blackadar' !< namelist parameter
    character(LEN=20) :: random_generator = 'random-parallel' !< namelist parameter
    character(LEN=20) :: reference_state = 'initial_profile' !< namelist parameter
    character(LEN=20) :: restart_data_format = 'mpi_shared_memory' !< namelist parameter
    character(LEN=20) :: restart_data_format_input = 'undefined' !< namelist parameter
    character(LEN=20) :: restart_data_format_output = 'undefined' !< namelist parameter
    character(LEN=20) :: timestep_scheme = 'runge-kutta-3' !< namelist parameter
    character(LEN=20) :: turbulence_closure = '1.5-order' !< namelist parameter
    character(LEN=23) :: origin_date_time = '2019-06-21 12:00:00 +00' !< date and time to be simulated
    character(LEN=40) :: flux_input_mode = 'application-specific' !< type of flux input: dynamic or kinematic
    character(LEN=40) :: flux_output_mode = 'application-specific' !< type of flux output: dynamic or kinematic
    character(LEN=40) :: topography = 'flat' !< namelist parameter
    character(LEN=64) :: host = '????' !< configuration identifier as given by palmrun option -c,
    !< ENVPAR namelist parameter provided by palmrun
    character(LEN=80) :: log_message !< user-defined message for debugging (sse data_log.f90)
    character(LEN=80) :: run_identifier !< run identifier as given by palmrun option -r, ENVPAR
    !< namelist parameter provided by palmrun
    character(LEN=100) :: initializing_actions = ' ' !< namelist parameter
    character(LEN=100) :: restart_string = ' ' !< for storing strings in case of writing/reading restart
    !< data
    character(LEN=210) :: run_description_header !< string containing diverse run informations as run
    !< identifier, coupling mode, host, ensemble number, run
    !< date and time
    character(LEN=1000) :: debug_string = ' ' !<.....
    character(LEN=1000) :: message_string = ' ' !< dynamic string for error message output

    character(LEN=varnamelength), dimension(300) :: data_output_pr_user = ' ' !< namelist parameter
    character(LEN=varnamelength), dimension(500) :: data_output_pr = ' ' !< namelist parameter
    character(LEN=varnamelength), dimension(500) :: data_output = ' ' !< namelist parameter
    character(LEN=varnamelength), dimension(500) :: data_output_user = ' ' !< namelist parameter
    character(LEN=varnamelength), dimension(500) :: doav = ' ' !< label array for multi-dimensional,
    !< averaged output quantities

    character(LEN=varnamelength), dimension(max_masks, 100) :: data_output_masks = ' ' !< namelist parameter
    character(LEN=varnamelength), dimension(max_masks, 100) :: data_output_masks_user = ' ' !< namelist parameter
    character(LEN=varnamelength), dimension(0:1, 500) :: do2d = ' ' !< label array for 2d output
    !< quantities
    character(LEN=varnamelength), dimension(0:1, 500) :: do3d = ' ' !< label array for 3d output
    !< quantities

    character(LEN=varnamelength), dimension(max_masks, 0:1, 100) :: domask = ' ' !< label array for multi-dimensional,
    !< masked output quantities

    integer(iwp) :: abort_mode = 1 !< abort condition (nested runs)
    integer(iwp) :: agt_time_count = 0 !< number of output intervals for agent data output
    integer(iwp) :: average_count_pr = 0 !< number of samples in vertical-profile output
    integer(iwp) :: average_count_3d = 0 !< number of samples in 3d output
    integer(iwp) :: current_timestep_number = 0 !< current timestep number, printed to RUN_CONTROL file
    integer(iwp) :: dist_range = 0 !< switch for steering the horizontal disturbance range, 1: inflow
    !< disturbances in case of non-cyclic horizontal BC, 0: otherwise
    integer(iwp) :: disturbance_level_ind_b !< lowest grid index where flow disturbance is applied
    integer(iwp) :: disturbance_level_ind_t !< highest grid index where flow disturbance is applied
    integer(iwp) :: doav_n = 0 !< number of 2d/3d output quantities subject to time averaging
    integer(iwp) :: dopr_n = 0 !< number of profile output quantities subject to time averaging
    integer(iwp) :: dopr_time_count = 0 !< number of output intervals for profile output
    integer(iwp) :: dots_time_count = 0 !< number of output intervals for timeseries output
    integer(iwp) :: dp_level_ind_b = 0 !< lowest grid index for external pressure gradient forcing
    integer(iwp) :: ensemble_member_nr = 0 !< namelist parameter
    integer(iwp) :: ibc_e_b !< integer flag for bc_e_b
    integer(iwp) :: ibc_p_b !< integer flag for bc_p_b
    integer(iwp) :: ibc_p_t !< integer flag for bc_p_t
    integer(iwp) :: ibc_pt_b !< integer flag for bc_pt_b
    integer(iwp) :: ibc_pt_t !< integer flag for bc_pt_t
    integer(iwp) :: ibc_q_b !< integer flag for bc_q_b
    integer(iwp) :: ibc_q_t !< integer flag for bc_q_t
    integer(iwp) :: ibc_s_b !< integer flag for bc_s_b
    integer(iwp) :: ibc_s_t !< integer flag for bc_s_t
    integer(iwp) :: ibc_uv_b !< integer flag for bc_uv_b
    integer(iwp) :: ibc_uv_t !< integer flag for bc_uv_t
    integer(iwp) :: inflow_disturbance_begin = -1 !< namelist parameter
    integer(iwp) :: inflow_disturbance_end = -1 !< namelist parameter
    integer(iwp) :: intermediate_timestep_count !< number of current Runge-Kutta substep
    integer(iwp) :: intermediate_timestep_count_max !< maximum number of Runge-Kutta substeps
    integer(iwp) :: io_group = 0 !< I/O group to which the PE belongs (= #PE / io_blocks)
    integer(iwp) :: io_blocks = 1 !< number of blocks for which I/O is done in sequence (total number of PEs /
    !< maximum_parallel_io_streams)
    integer(iwp) :: iran = -1234567 !< integer random number used for flow disturbances
    integer(iwp) :: length = 0 !< integer that specifies the length of a string in case of writing/reading
    !< restart data
    integer(iwp) :: masks = 0 !< counter for number of masked output quantities
    integer(iwp) :: maximum_parallel_io_streams = -1 !< maximum number of parallel io streams that the underlying parallel file
    !< system allows, set with palmrun option -w, ENVPAR namelist parameter, provided by palmrun
    integer(iwp) :: max_pr_cs = 0 !< number of chemistry profiles in output
    integer(iwp) :: max_pr_det = 0 !< number of det profiles in output
    integer(iwp) :: max_pr_salsa = 0 !< number of salsa profiles (must not change within a job chain)
    integer(iwp) :: max_pr_user = 0 !< number of user-defined profiles (must not change within a job chain)
    integer(iwp) :: max_pr_user_tmp = 0 !< number of user-defined profiles that is temporary stored to check it
    !< against max_pr_user in case of restarts
    integer(iwp) :: nr_timesteps_this_run = 0 !< number of timesteps (cpu time measurements)
    integer(iwp) :: nsor = 20 !< namelist parameter
    integer(iwp) :: nsor_ini = 100 !< namelist parameter
    integer(iwp) :: n_sor !< number of iterations to be used in SOR-scheme
    integer(iwp) :: normalizing_region = 0 !< namelist parameter
    integer(iwp) :: num_mean_inflow_profiles = 7 !< number of mean inflow profiles in case of turbulent inflow
    integer(iwp) :: num_leg = 0 !< number of different legs in virtual flight measurements
    integer(iwp) :: num_var_fl !< number of sampling/output variables in virtual flight measurements
    integer(iwp) :: num_var_fl_user = 0 !< number of user-defined sampling/output variables in virtual flight
    !< measurements
    integer(iwp) :: number_dz !< number of user-specified dz values
    integer(iwp) :: number_stretch_level_start !< number of user-specified start levels for stretching
    integer(iwp) :: number_stretch_level_end !< number of user-specified end levels for stretching
    integer(iwp) :: nz_do3d = -9999 !< namelist parameter
    integer(iwp) :: prt_time_count = 0 !< number of output intervals for particle data output
    integer(iwp) :: runnr = 0 !< number of run in job chain
    integer(iwp) :: symmetry_flag = 0 !< flag for sterring the symmetric behavior of the bottom and top boundary
    integer(iwp) :: terminate_coupled = 0 !< switch for steering termination in case of coupled runs
    integer(iwp) :: timestep_count = 0 !< number of timesteps carried out since the beginning of the initial run
    integer(iwp) :: y_shift = 0 !< namelist parameter

    integer(iwp) :: dist_nxl(0:1) !< left boundary of disturbance region
    integer(iwp) :: dist_nxr(0:1) !< right boundary of disturbance region
    integer(iwp) :: dist_nyn(0:1) !< north boundary of disturbance region
    integer(iwp) :: dist_nys(0:1) !< south boundary of disturbance region
    integer(iwp) :: do2d_no(0:1) = 0 !< number of 2d output quantities
    integer(iwp) :: do2d_xy_time_count(0:1) = 0 !< number of output intervals for 2d data (xy)
    integer(iwp) :: do2d_xz_time_count(0:1) = 0 !< number of output intervals for 2d data (xz)
    integer(iwp) :: do2d_yz_time_count(0:1) = 0 !< number of output intervals for 2d data (yz)
    integer(iwp) :: do3d_no(0:1) = 0 !< number of 3d output quantities
    integer(iwp) :: do3d_time_count(0:1) = 0 !< number of output intervals for 3d data
    integer(iwp) :: domask_no(max_masks, 0:1) = 0 !< number of masked output quantities
    integer(iwp) :: domask_time_count(max_masks, 0:1) !< number of output intervals for masked data
    integer(iwp) :: dz_stretch_level_end_index(9) !< vertical grid level index until which the vertical grid spacing
    !< is stretched
    integer(iwp) :: dz_stretch_level_start_index(9) !< vertical grid level index above which the vertical grid spacing
    !< is stretched
    integer(iwp) :: mask_size(max_masks, 3) = -1 !< size of mask array per mask and dimension (for netcdf output)
    integer(iwp) :: mask_size_l(max_masks, 3) = -1 !< subdomain size of mask array per mask and dimension
    !< (for netcdf output)
    integer(iwp) :: mask_start_l(max_masks, 3) = -1 !< subdomain start index of mask array (for netcdf output)
    integer(iwp) :: pt_vertical_gradient_level_ind(12) = -9999 !< grid index values of pt_vertical_gradient_level(s)
    integer(iwp) :: q_vertical_gradient_level_ind(12) = -9999 !< grid index values of q_vertical_gradient_level(s)
    integer(iwp) :: s_vertical_gradient_level_ind(12) = -9999 !< grid index values of s_vertical_gradient_level(s)
    integer(iwp) :: section(100, 3) !< collective array for section_xy/xz/yz
    integer(iwp) :: section_xy(100) = -9999 !< namelist parameter
    integer(iwp) :: section_xz(100) = -9999 !< namelist parameter
    integer(iwp) :: section_yz(100) = -9999 !< namelist parameter
    integer(iwp) :: ug_vertical_gradient_level_ind(12) = -9999 !< grid index values of ug_vertical_gradient_level(s)
    integer(iwp) :: vg_vertical_gradient_level_ind(12) = -9999 !< grid index values of vg_vertical_gradient_level(s)
    integer(iwp) :: subs_vertical_gradient_level_i(12) = -9999 !< grid index values of subs_vertical_gradient_level(s)
    integer(iwp), dimension(0:1) :: ntdim_2d_xy !< number of output intervals for 2d data (xy)
    integer(iwp), dimension(0:1) :: ntdim_2d_xz !< number of output intervals for 2d data (xz)
    integer(iwp), dimension(0:1) :: ntdim_2d_yz !< number of output intervals for 2d data (yz)
    integer(iwp), dimension(0:1) :: ntdim_3d !< number of output intervals for 3d data

    integer(iwp), dimension(max_masks, 0:1) :: ntdim_masks !< number of output intervals for masked data

    integer(iwp), dimension(max_masks, mask_xyz_dimension) :: mask_k_over_surface = -1 !< namelist parameter, k index of height
    !<over surface

    integer(iwp), dimension(:, :), allocatable :: mask_i !< subdomain grid index of masked output point on x-dimension
    integer(iwp), dimension(:, :), allocatable :: mask_j !< subdomain grid index of masked output point on y-dimension
    integer(iwp), dimension(:, :), allocatable :: mask_k !< subdomain grid index of masked output point on z-dimension
    integer(iwp), dimension(:, :), allocatable :: mask_i_global !< global grid index of masked output point on x-dimension
    integer(iwp), dimension(:, :), allocatable :: mask_j_global !< global grid index of masked output point on y-dimension
    integer(iwp), dimension(:, :), allocatable :: mask_k_global !< global grid index of masked output point on z-dimension

    logical :: advanced_div_correction = .false. !< namelist parameter
    logical :: agent_time_unlimited = .false. !< namelist parameter
    logical :: air_chemistry = .false. !< chemistry model switch
    logical :: allow_negative_scalar_values = .true. !< namelist parameter
    logical :: allow_roughness_limitation = .false. !< namelist parameter
    logical :: atmosphere_run_coupled_to_ocean = .false. !< atmosphere part of a coupled atmosphere-ocean run
    logical :: bc_dirichlet_l = .false. !< flag indicating dirichlet boundary condition on left model
    !< boundary
    logical :: bc_dirichlet_n = .false. !< flag indicating dirichlet boundary condition on north model
    !< boundary
    logical :: bc_dirichlet_r = .false. !< flag indicating dirichlet boundary condition on right model
    !< boundary
    logical :: bc_dirichlet_s = .false. !< flag indicating dirichlet boundary condition on south model
    !< boundary
    logical :: bc_lr_cyc = .true. !< left-right boundary condition cyclic?
    logical :: bc_lr_dirrad = .false. !< left-right boundary condition dirichlet/radiation?
    logical :: bc_lr_raddir = .false. !< left-right boundary condition radiation/dirichlet?
    logical :: bc_ns_cyc = .true. !< north-south boundary condition cyclic?
    logical :: bc_ns_dirrad = .false. !< north-south boundary condition dirichlet/radiation?
    logical :: bc_ns_raddir = .false. !< north-south boundary condition radiation/dirichlet?
    logical :: bc_radiation_l = .false. !< radiation boundary condition for outflow at left domain boundary
    logical :: bc_radiation_n = .false. !< radiation boundary condition for outflow at north domain
    !< boundary
    logical :: bc_radiation_r = .false. !< radiation boundary condition for outflow at right domain
    !< boundary
    logical :: bc_radiation_s = .false. !< radiation boundary condition for outflow at south domain
    !< boundary
    logical :: biometeorology = .false. !< biometeorology module switch
    logical :: calc_soil_moisture_during_spinup = .false. !< namelist parameter
    logical :: call_psolver_at_all_substeps = .true. !< namelist parameter
    logical :: check_realistic_q = .true. !< namelist parameter
    logical :: child_domain = .false. !< flag indicating that model is nested in a parent domain
    logical :: cloud_droplets = .false. !< namelist parameter
    logical :: conserve_volume_flow = .false. !< namelist parameter
    logical :: constant_diffusion = .false. !< diffusion coefficient constant?
    logical :: constant_flux_layer = .true. !< namelist parameter
    logical :: constant_heatflux = .true. !< heat flux at all surfaces constant?
    logical :: constant_top_heatflux = .true. !< heat flux at domain top constant?
    logical :: constant_top_momentumflux = .false. !< momentum flux at domain topconstant?
    logical :: constant_top_salinityflux = .true. !< constant salinity flux at ocean surface
    logical :: constant_top_scalarflux = .true. !< passive-scalar flux at domain top constant?
    logical :: constant_scalarflux = .true. !< passive-scalar flux at surfaces constant?
    logical :: constant_waterflux = .true. !< water flux at all surfaces constant?
    logical :: create_disturbances = .true. !< namelist parameter
    logical :: cut_cell_topography = .false. !< namelist parameter
    logical :: cyclic_fill_initialization = .false. !< switch for steering cyclic fill actions
    logical :: dmp_enabled = .false. !< namelist parameter
    logical :: data_output_during_spinup = .false. !< namelist parameter
    logical :: data_output_raw = .false. !< namelist parameter
    logical :: data_output_2d_on_each_pe = .true. !< namelist parameter
    logical :: dcep = .false. !< switch for activiating the dcep model
    logical :: debug_output = .false. !< namelist parameter
    logical :: debug_output_timestep = .false. !< namelist parameter
    logical :: disturbance_created = .false. !< flow disturbance imposed?
    logical :: do2d_at_begin = .false. !< namelist parameter
    logical :: do3d_at_begin = .false. !< namelist parameter
    logical :: do_sum = .false. !< contribute to time average of profile data?
    logical :: dp_external = .false. !< namelist parameter
    logical :: dp_smooth = .false. !< namelist parameter
    logical :: dt_fixed = .false. !< fixed timestep (namelist parameter dt set)?
    logical :: dt_3d_reached !< internal timestep for particle advection
    logical :: dt_3d_reached_l !< internal timestep for particle advection
    logical :: det_enabled = .false. !< switch for the dust transport and emission module (DET)
    logical :: enable_openacc = .true.
    !< namelist parameter, enables OpenACC offloading, only meaningful if built with OpenACC support
    logical :: fct_enabled = .false. !< switch for activating the flow control module
    logical :: first_call_mas = .true. !< call mas only once per timestep ??
    logical :: force_print_header = .false. !< namelist parameter
    logical :: galilei_transformation = .false. !< namelist parameter
    logical :: homogenize_surface_temperature = .false. !< namelist parameter
    logical :: humidity = .false. !< namelist parameter
    logical :: humidity_remote = .false. !< switch for receiving near-surface humidity flux
    !< (atmosphere-ocean coupling)
    logical :: implicit_diffusion_1d = .false. !< Crank-Nicolson scheme for diffusion term in 1d-model
    logical :: indoor_model = .false. !< switch for indoor-climate and energy-demand model
    logical :: interpolate_to_grid_center = .false. !< namelist parameter
    logical :: kolmogorov_length_scale = .false. !< switch to activate calculations in flow_statistics for the
    !< kolmogorov length scale
    logical :: large_scale_forcing = .false. !< namelist parameter
    logical :: large_scale_subsidence = .false. !< namelist parameter
    logical :: land_surface = .false. !< use land surface model?
    logical :: les_dai = .false. !< use Dai et al. turbulence closure (modified 1.5-order closure)
    !< for LES mode. Shall replace the default 1.5-order closure
    logical :: les_dynamic = .false. !< use dynamic subgrid model as turbulence closure for LES mode
    logical :: les_default = .false. !< use 1.5-order default turbulence closure for LES mode
    logical :: lsf_exception = .false. !< use of lsf with buildings (temporary)?
    logical :: lsf_surf = .true. !< use surface forcing (large scale forcing)?
    logical :: lsf_vert = .true. !< use atmospheric forcing (large scale forcing)?
    logical :: masking_method = .true. !< namelist parameter
    logical :: monotonic_limiter_z = .false. !< use monotonic flux limiter for vertical scalar advection
    logical :: nested_run = .false. !< general switch to indicate a nested run
    logical :: nesting_offline = .false. !< flag controlling offline nesting in COSMO model
    logical :: neutral = .false. !< namelist parameter
    logical :: nudging = .false. !< namelist parameter
    logical :: ocean_mode = .false. !< namelist parameter
    logical :: ocean_run_coupled_to_atmosphere = .false. !< ocean part of coupled atmosphere-ocean run
    logical :: open_debug_files = .true. !< switch for opening the debug files independent of switch debug_output
    logical :: passive_scalar = .false. !< namelist parameter
    logical :: pe_grid_prescribed = .false. !< switch to indicate if PE grid is prescribed by user
    logical :: plant_canopy = .false. !< switch for use of plant canopy model
    logical :: random_heatflux = .false. !< namelist parameter
    logical :: rans_mode = .false. !< switch between RANS and LES mode
    logical :: rans_tke_e = .false. !< use TKE-e turbulence closure for RANS mode
    logical :: rans_tke_l = .false. !< use TKE-l turbulence closure for RANS mode
    logical :: read_spinup_data = .false. !< flag to control the input of surface spinup data
    logical :: read_svf = .false. !< ENVPAR namelist parameter to steer input of svf
    !< (ENVPAR is provided by palmrun)
    logical :: run_control_header = .false. !< onetime output of RUN_CONTROL header
    logical :: salinity = .true. !< switch for using salinity
    logical :: salsa = .false. !< switch for the sectional aerosol module salsa
    logical :: scalar_rayleigh_damping = .true. !< namelist parameter
    logical :: serial_run = .false. !< switch to indicate serial or parallel run
    logical :: sloping_surface = .false. !< use sloped surface? (namelist parameter alpha_surface)
    logical :: slurb = .false. !< switch for activiating the SLUrb model
    logical :: spinup = .false. !< perform model spinup without atmosphere code?
    logical :: spinup_phase = .false. !< indicates, if model is in spinup phase
    logical :: surface_output = .false. !< output of surface data
    logical :: stop_dt = .false. !< internal switch to stop the time stepping
    logical :: synchronous_exchange = .false. !< namelist parameter
    logical :: syn_turb_gen = .false. !< flag for synthetic turbulence generator module
    logical :: temperton_fft_vec = .false. !< flag for using vectorized version of Temperton FFT
    logical :: terminate_run = .false. !< terminate run (cpu-time limit, restarts)?
    logical :: terrain_following_mapping = .false. !< namelist parameter
    logical :: topo_no_distinct = .false. !< flag controlling classification of topography surfaces
    logical :: turbulent_inflow = .false. !< namelist parameter
    logical :: turbulent_outflow = .false. !< namelist parameter
    logical :: urban_surface = .false. !< use urban surface model?
    logical :: use_contiguous_buffer = .false. !< namelist parameter
    logical :: use_fixed_date = .false. !< date of simulation does not change (namelist parameter)
    logical :: use_fixed_time = .false. !< time of simulation does not change (namelist parameter)
    logical :: use_free_convection_scaling = .false. !< namelist parameter to switch on free convection velocity scale
    !< in calculation of horizontal wind speed (surface_layer_fluxes)
    logical :: use_initial_profile_as_reference = .false. !< use of initial profiles as reference state?
    logical :: use_prescribed_profile_data = .false. !< use of prescribed wind profiles?
    !< (namelist parameters u_profile, v_profile)
    logical :: use_single_reference_value = .false. !< use of single value as reference state?
    logical :: use_sm_for_poisfft = .false.
    !< use shared-memory on node for solving the Poisson equation with FFT-method (1d-decomposition)
    logical :: use_subsidence_tendencies = .false. !< namelist parameter
    logical :: use_surface_fluxes = .false. !< namelist parameter
    logical :: use_top_fluxes = .false. !< namelist parameter
    logical :: use_ug_for_galilei_tr = .true. !< namelist parameter
    logical :: use_upstream_for_tke = .false. !< namelist parameter
    logical :: uv_radiation = .false. !< use UV-radiation module
    logical :: vdi_checks = .false. !< do internal controls after VDI 3783 Part 9
    logical :: traffic = .false. !< use traffic module
    logical :: virtual_flight = .false. !< use virtual flight model
    logical :: virtual_measurement = .false. !< control parameter to switch-on virtual measurements
    logical :: wall_adjustment = .true. !< namelist parameter
    logical :: wind_turbine = .false. !< flag for use of wind turbine model
    logical :: write_binary = .false. !< ENVPAR namelist parameter to steer restart I/O
    !< (ENVPAR is provided by palmrun)
    logical :: write_spinup_data = .false. !< ENVPAR namelist parameter to steer restart I/O
    !< (ENVPAR is provided by palmrun)
    logical :: write_svf = .false. !< ENVPAR namelist parameter to steer output of svf
    !< (ENVPAR is provided by palmrun)
    logical :: ws_scheme_sca = .false. !< use Wicker-Skamarock scheme (scalar advection)?
    logical :: ws_scheme_mom = .false. !< use Wicker-Skamarock scheme (momentum advection)?

    logical :: data_output_xy(0:1) = .false. !< output of xy cross-section data?
    logical :: data_output_xz(0:1) = .false. !< output of xz cross-section data?
    logical :: data_output_yz(0:1) = .false. !< output of yz cross-section data?

    logical, dimension(max_masks) :: mask_surface = .false. !< flag for surface-following masked output

    real(wp) :: advected_distance_x = 0.0_wp !< advected distance of model domain along x
    !< (galilei transformation)
    real(wp) :: advected_distance_y = 0.0_wp !< advected distance of model domain along y
    !< (galilei transformation)
    real(wp) :: alpha_surface = 0.0_wp !< namelist parameter
    real(wp) :: atmos_ocean_sign = 1.0_wp !< vertical-grid conversion factor
    !< (=1.0 in atmosphere, =-1.0 in ocean)
    real(wp) :: averaging_interval = 0.0_wp !< namelist parameter
    real(wp) :: averaging_interval_pr = 9999999.9_wp !< namelist parameter
    real(wp) :: bc_pt_t_val !< vertical gradient of pt near domain top
    real(wp) :: bc_q_t_val !< vertical gradient of humidity near domain top
    real(wp) :: bc_s_t_val !< vertical gradient of passive scalar near domain top
    real(wp) :: bottom_salinityflux = 0.0_wp !< namelist parameter
    real(wp) :: building_height = 50.0_wp !< namelist parameter
    real(wp) :: building_length_x = 50.0_wp !< namelist parameter
    real(wp) :: building_length_y = 50.0_wp !< namelist parameter
    real(wp) :: building_wall_left = 9999999.9_wp !< namelist parameter
    real(wp) :: building_wall_south = 9999999.9_wp !< namelist parameter
    real(wp) :: canyon_height = 50.0_wp !< namelist parameter
    real(wp) :: canyon_width_x = 9999999.9_wp !< namelist parameter
    real(wp) :: canyon_width_y = 9999999.9_wp !< namelist parameter
    real(wp) :: canyon_wall_left = 9999999.9_wp !< namelist parameter
    real(wp) :: canyon_wall_south = 9999999.9_wp !< namelist parameter
    real(wp) :: cfl_factor = -1.0_wp !< namelist parameter
    real(wp) :: cos_alpha_surface !< cosine of alpha_surface
    real(wp) :: coupling_start_time = 0.0_wp !< namelist parameter
    real(wp) :: days_since_reference_point = 0.0_wp !< days after atmosphere-ocean coupling has been activated,
    !< or after spinup phase of LSM has been finished
    real(wp) :: disturbance_amplitude = 0.25_wp !< namelist parameter
    real(wp) :: disturbance_energy_limit = 0.01_wp !< namelist parameter
    real(wp) :: disturbance_level_b = -9999999.9_wp !< namelist parameter
    real(wp) :: disturbance_level_t = -9999999.9_wp !< namelist parameter
    real(wp) :: dp_level_b = 0.0_wp !< namelist parameter
    real(wp) :: dt = -1.0_wp !< namelist parameter
    real(wp) :: dt_averaging_input = 0.0_wp !< namelist parameter
    real(wp) :: dt_averaging_input_pr = 9999999.9_wp !< namelist parameter
    real(wp) :: dt_coupling = 9999999.9_wp !< namelist parameter
    real(wp) :: dt_data_output = 9999999.9_wp !< namelist parameter
    real(wp) :: dt_data_output_av = 9999999.9_wp !< namelist parameter
    real(wp) :: dt_disturb = 9999999.9_wp !< namelist parameter
    real(wp) :: dt_dopr = 9999999.9_wp !< namelist parameter
    real(wp) :: dt_dopr_listing = 9999999.9_wp !< namelist parameter
    real(wp) :: dt_dots = 9999999.9_wp !< namelist parameter
    real(wp) :: dt_do2d_xy = 9999999.9_wp !< namelist parameter
    real(wp) :: dt_do2d_xz = 9999999.9_wp !< namelist parameter
    real(wp) :: dt_do2d_yz = 9999999.9_wp !< namelist parameter
    real(wp) :: dt_do3d = 9999999.9_wp !< namelist parameter
    real(wp) :: dt_max = 20.0_wp !< namelist parameter
    real(wp) :: dt_restart = 9999999.9_wp !< namelist parameter
    real(wp) :: dt_run_control = 60.0_wp !< namelist parameter
    real(wp) :: dt_run_control_spinup = 3600.0_wp !< namelist parameter
    real(wp) :: dt_write_agent_data = 9999999.9_wp !< namelist parameter
    real(wp) :: dt_3d = 0.01_wp !< time step
    real(wp) :: dz_max = 999.0_wp !< namelist parameter
    real(wp) :: dz_stretch_factor = 1.08_wp !< namelist parameter
    real(wp) :: dz_stretch_level = -9999999.9_wp !< namelist parameter
    real(wp) :: e_init = 0.0_wp !< namelist parameter
    real(wp) :: e_min = 0.0_wp !< namelist parameter
    real(wp) :: end_time = 0.0_wp !< namelist parameter
    real(wp) :: f = 0.0_wp !< Coriolis parameter
    real(wp) :: fs = 0.0_wp !< Coriolis parameter
    real(wp) :: implicit_timestep_factor = 5.0_wp !< factor by which the timestep is enlarged when
    !< using Crank-Nicolson scheme
    real(wp) :: km_constant = -1.0_wp !< namelist parameter
    real(wp) :: latitude = 55.0_wp !< namelist parameter
    real(wp) :: longitude = 0.0_wp !< namelist parameter
    real(wp) :: mask_scale_x = 1.0_wp !< namelist parameter
    real(wp) :: mask_scale_y = 1.0_wp !< namelist parameter
    real(wp) :: mask_scale_z = 1.0_wp !< namelist parameter
    real(wp) :: maximum_cpu_time_allowed = 0.0_wp !< given wall time for run
    real(wp) :: molecular_viscosity = 1.461e-5_wp !< molecular viscosity (used in lsm and lpm)
    real(wp) :: multi_agent_system_end = 9999999.9_wp !< namelist parameter (see documentation)
    real(wp) :: multi_agent_system_start = 0.0_wp !< namelist parameter (see documentation)
    real(wp) :: omega = 7.29212e-5_wp !< namelist parameter
    real(wp) :: omega_sor = 1.8_wp !< namelist parameter
    real(wp) :: outflow_damping_factor = 0.0_wp !< namelist parameter
    real(wp) :: outflow_damping_width = 0.0_wp !< namelist parameter
    real(wp) :: outflow_source_plane = -9999999.9_wp !< namelist parameter
    real(wp) :: output_3d_file_size = 0.0_wp !< file size of 3d NetCDF output file
    real(wp) :: particle_maximum_age = 9999999.9_wp !< namelist parameter
    real(wp) :: prandtl_number = 1.0_wp !< namelist parameter
    real(wp) :: pt_damping_factor = 0.0_wp !< namelist parameter
    real(wp) :: pt_damping_width = 0.0_wp !< namelist parameter
    real(wp) :: pt_reference = 9999999.9_wp !< namelist parameter
    real(wp) :: pt_slope_offset = 0.0_wp !< temperature difference between left and right
    !< boundary of total domain
    real(wp) :: pt_surface = 300.0_wp !< namelist parameter
    real(wp) :: pt_surface_heating_rate = 0.0_wp !< namelist parameter
    real(wp) :: pt_surface_initial_change = 0.0_wp !< namelist parameter
    real(wp) :: q_surface = 0.0_wp !< namelist parameter
    real(wp) :: q_surface_initial_change = 0.0_wp !< namelist parameter
    real(wp) :: rayleigh_damping_factor = 0.0_wp !< namelist parameter
    real(wp) :: rayleigh_damping_height = -1.0_wp !< namelist parameter
    real(wp) :: restart_file_size !< size of restart file in mbyte
    real(wp) :: restart_time = 9999999.9_wp !< namelist parameter
    real(wp) :: rho_cp !< rho at surface * c_p
    real(wp) :: rho_reference !< reference state of density
    real(wp) :: rho_surface !< surface value of density
    real(wp) :: rotation_angle = 0.0_wp !< clockwise rotation of model North relative to real North [deg]
    real(wp) :: roughness_length = 0.1_wp !< namelist parameter
    real(wp) :: simulated_time = 0.0_wp !< elapsed simulated time
    real(wp) :: simulated_time_at_begin !< elapsed simulated time of previous run (job chain)
    real(wp) :: sin_alpha_surface !< sine of alpha_surface (sloped surface)
    real(wp) :: skip_time_data_output = 0.0_wp !< namelist parameter
    real(wp) :: skip_time_data_output_av = 9999999.9_wp !< namelist parameter
    real(wp) :: skip_time_dopr = 9999999.9_wp !< namelist parameter
    real(wp) :: skip_time_do2d_xy = 9999999.9_wp !< namelist parameter
    real(wp) :: skip_time_do2d_xz = 9999999.9_wp !< namelist parameter
    real(wp) :: skip_time_do2d_yz = 9999999.9_wp !< namelist parameter
    real(wp) :: skip_time_do3d = 9999999.9_wp !< namelist parameter
    real(wp) :: spinup_pt_amplitude = 0.0_wp !< namelist parameter
    real(wp) :: spinup_pt_mean = 9999999.9_wp !< namelist parameter
    real(wp) :: spinup_time = 0.0_wp !< namelist parameter
    real(wp) :: surface_heatflux = 9999999.9_wp !< namelist parameter
    real(wp) :: surface_pressure = 1013.25_wp !< namelist parameter
    real(wp) :: surface_scalarflux = 9999999.9_wp !< namelist parameter
    real(wp) :: surface_waterflux = 9999999.9_wp !< namelist parameter
    real(wp) :: s_surface = 0.0_wp !< namelist parameter
    real(wp) :: s_surface_initial_change = 0.0_wp !< namelist parameter
    real(wp) :: termination_time_needed = 35.0_wp !< namelist parameter
    real(wp) :: time_coupling = 0.0_wp !< time since last coupling (surface_coupler)
    real(wp) :: time_disturb = 0.0_wp !< time since last flow disturbance
    real(wp) :: time_dopr = 0.0_wp !< time since last profile output
    real(wp) :: time_dopr_av = 0.0_wp !< time since last averaged profile output
    real(wp) :: time_dopr_listing = 0.0_wp !< time since last profile output (ASCII) on file
    real(wp) :: time_dosp = 0.0_wp !< time since last spectra output
    real(wp) :: time_dosp_av = 0.0_wp !< time since last averaged spectra output
    real(wp) :: time_dots = 0.0_wp !< time since last timeseries output
    real(wp) :: time_do2d_xy = 0.0_wp !< time since last xy cross-section output
    real(wp) :: time_do2d_xz = 0.0_wp !< time since last xz cross-section output
    real(wp) :: time_do2d_yz = 0.0_wp !< time since last yz cross-section output
    real(wp) :: time_do3d = 0.0_wp !< time since last 3d output
    real(wp) :: time_do_av = 0.0_wp !< time since last averaged-data output
    real(wp) :: time_do_sla = 0.0_wp !< time since last
    real(wp) :: time_restart = 9999999.9_wp !< time at which run shall be terminated and restarted
    real(wp) :: time_run_control = 0.0_wp !< time since last RUN_CONTROL output
    real(wp) :: time_since_reference_point = 0.0_wp !< time after atmosphere-ocean coupling has been activated, or time
    !< after spinup phase of LSM has been finished
    real(wp) :: top_heatflux = 9999999.9_wp !< namelist parameter
    real(wp) :: top_momentumflux_u = 9999999.9_wp !< namelist parameter
    real(wp) :: top_momentumflux_v = 9999999.9_wp !< namelist parameter
    real(wp) :: top_salinityflux = 9999999.9_wp !< namelist parameter
    real(wp) :: top_scalarflux = 9999999.9_wp !< namelist parameter
    real(wp) :: tunnel_height = 9999999.9_wp !< namelist parameter
    real(wp) :: tunnel_length = 9999999.9_wp !< namelist parameter
    real(wp) :: tunnel_width_x = 9999999.9_wp !< namelist parameter
    real(wp) :: tunnel_width_y = 9999999.9_wp !< namelist parameter
    real(wp) :: tunnel_wall_depth = 9999999.9_wp !< namelist parameter
    real(wp) :: ug_surface = 0.0_wp !< namelist parameter
    real(wp) :: u_bulk = 0.0_wp !< namelist parameter
    real(wp) :: u_gtrans = 0.0_wp !< transformed wind component (galilei transformation)
    real(wp) :: vg_surface = 0.0_wp !< namelist parameter
    real(wp) :: vpt_reference = 9999999.9_wp !< reference state of virtual potential temperature
    real(wp) :: v_bulk = 0.0_wp !< namelist parameter
    real(wp) :: v_gtrans = 0.0_wp !< transformed wind component (galilei transformation)
    real(wp) :: zeta_max = 20.0_wp !< namelist parameter
    real(wp) :: zeta_min = -20.0_wp !< namelist parameter
    real(wp) :: z0h_factor = 1.0_wp !< namelist parameter

    real(wp) :: do2d_xy_last_time(0:1) = -1.0_wp !< time of previous xy output
    real(wp) :: do2d_xz_last_time(0:1) = -1.0_wp !< time of previous xz output
    real(wp) :: do2d_yz_last_time(0:1) = -1.0_wp !< time of previous yz output
    real(wp) :: dpdxy(1:2) = 0.0_wp !< namelist parameter
    real(wp) :: dt_domask(max_masks) = 9999999.9_wp !< namelist parameter
    real(wp) :: dz(10) = -1.0_wp !< namelist parameter
    real(wp) :: dz_stretch_level_start(9) = -9999999.9_wp !< namelist parameter
    real(wp) :: dz_stretch_level_end(9) = 9999999.9_wp !< namelist parameter
    real(wp) :: dz_stretch_factor_array(9) = 1.08_wp !< namelist parameter
    real(wp) :: mask_scale(3) !< collective array for mask_scale_x/y/z
    real(wp) :: pt_vertical_gradient(12) = 0.0_wp
    !< namelist parameter (dimensioned with 12 to store two extra points at bottom/top)
    real(wp) :: pt_vertical_gradient_level(12) = huge(1.0_wp)
    !< namelist parameter (dimensioned with 12 to store two extra points at bottom/top)
    real(wp) :: q_vertical_gradient(12) = 0.0_wp !< namelist parameter (dimensioned with 12 to store two extra points at bottom/top)
    real(wp) :: q_vertical_gradient_level(12) = huge(1.0_wp)
    !< namelist parameter (dimensioned with 12 to store two extra points at bottom/top)
    real(wp) :: section_xy_m(100) = -9999999.9_wp !< namelist parameter
    real(wp) :: section_xz_m(100) = -9999999.9_wp !< namelist parameter
    real(wp) :: section_yz_m(100) = -9999999.9_wp !< namelist parameter
    real(wp) :: s_vertical_gradient(12) = 0.0_wp !< namelist parameter (dimensioned with 12 to store two extra points at bottom/top)
    real(wp) :: s_vertical_gradient_level(12) = huge(1.0_wp)
    !< namelist parameter (dimensioned with 12 to store two extra points at bottom/top)
    real(wp) :: skip_time_domask(max_masks) = 9999999.9_wp !< namelist parameter
    real(wp) :: threshold(20) = 0.0_wp !< namelist parameter
    real(wp) :: time_domask(max_masks) = 0.0_wp !< namelist parameter
    real(wp) :: tsc(10) = (/1.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, & !< array used for controlling time-integration at different substeps
                            0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp/)
    real(wp) :: u_profile(200) = 9999999.9_wp !< namelist parameter
    real(wp) :: uv_heights(200) = 9999999.9_wp !< namelist parameter
    real(wp) :: v_profile(200) = 9999999.9_wp !< namelist parameter
    real(wp) :: ug_vertical_gradient(12) = 0.0_wp
    !< namelist parameter (dimensioned with 12 to store two extra points at bottom/top)
    real(wp) :: ug_vertical_gradient_level(12) = huge(1.0_wp)
    !< namelist parameter (dimensioned with 12 to store two extra points at bottom/top)
    real(wp) :: vg_vertical_gradient(12) = 0.0_wp
    !< namelist parameter (dimensioned with 12 to store two extra points at bottom/top)
    real(wp) :: vg_vertical_gradient_level(12) = huge(1.0_wp)
    !< namelist parameter (dimensioned with 12 to store two extra points at bottom/top)
    real(wp) :: volume_flow(1:3) = 0.0_wp !< volume flow through 1:yz-plane, 2: xz-plane, 3: xy-plane
    !< (nest childs only)
    real(wp) :: volume_flow_area(1:3) = 0.0_wp !< area of the respective volume flow planes
    real(wp) :: volume_flow_initial(1:3) = 0.0_wp !< initial volume flow (t=0) through the respective volume flow
    !< planes
    real(wp) :: wall_heatflux(0:5) = 0.0_wp !< namelist parameter
    real(wp) :: wall_humidityflux(0:5) = 0.0_wp !< namelist parameter
    real(wp) :: wall_scalarflux(0:5) = 0.0_wp !< namelist parameter
    real(wp) :: subs_vertical_gradient(12) = 0.0_wp
    !< namelist parameter (dimensioned with 12 to store two extra points at bottom/top)
    real(wp) :: subs_vertical_gradient_level(12) = huge(1.0_wp)
    !< namelist parameter (dimensioned with 12 to store two extra points at bottom/top)

    real(wp), dimension(:), allocatable :: dp_smooth_factor !< smoothing factor for external pressure gradient forcing

    real(wp), dimension(max_masks, mask_xyz_dimension) :: mask_x = -1.0_wp !< namelist parameter
    real(wp), dimension(max_masks, mask_xyz_dimension) :: mask_y = -1.0_wp !< namelist parameter
    real(wp), dimension(max_masks, mask_xyz_dimension) :: mask_z = -1.0_wp !< namelist parameter

    real(wp), dimension(max_masks, 3) :: mask_x_loop = -1.0_wp !< namelist parameter
    real(wp), dimension(max_masks, 3) :: mask_y_loop = -1.0_wp !< namelist parameter
    real(wp), dimension(max_masks, 3) :: mask_z_loop = -1.0_wp !< namelist parameter

!
!-- Internal mask arrays ("mask,dimension,selection")
    real(wp), dimension(:, :, :), allocatable :: mask !< collective array for mask_x/y/z
    real(wp), dimension(:, :, :), allocatable :: mask_loop !< collective array for mask_x/y/z_loop

    save

end module control_parameters