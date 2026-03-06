!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Definition of all arrays defined on the computational grid.
!--------------------------------------------------------------------------------------------------!
module arrays_3d

    use kinds

    real(wp), dimension(:), allocatable :: d_exner !< ratio of potential and actual temperature
    real(wp), dimension(:), allocatable :: ddzu !< 1/dzu
    real(wp), dimension(:), allocatable :: ddzu_pres !< modified ddzu for pressure solver
    real(wp), dimension(:), allocatable :: dd2zu !< 1/(dzu(k)+dzu(k+1))
    real(wp), dimension(:), allocatable :: drho_air !< inverse air density profile on the uv grid
    real(wp), dimension(:), allocatable :: drho_air_zw !< inverse air density profile on the w grid
    real(wp), dimension(:), allocatable :: dzu !< vertical grid size (u-grid)
    real(wp), dimension(:), allocatable :: ddzw !< 1/dzw
    real(wp), dimension(:), allocatable :: dzw !< vertical grid size (w-grid)
    real(wp), dimension(:), allocatable :: exner !< ratio of actual and potential temperature
    real(wp), dimension(:), allocatable :: heatflux_input_conversion !< conversion factor array for heatflux input
    real(wp), dimension(:), allocatable :: heatflux_output_conversion !< conversion factor array for heatflux output
    real(wp), dimension(:), allocatable :: hyp !< hydrostatic pressure
    real(wp), dimension(:), allocatable :: hyrho !< density of air calculated with hydrostatic pressure
    real(wp), dimension(:), allocatable :: momentumflux_input_conversion !< conversion factor array for momentumflux input
    real(wp), dimension(:), allocatable :: momentumflux_output_conversion !< conversion factor array for momentumflux output
    real(wp), dimension(:), allocatable :: odf_x !< outflow damping factor in x-direction
    real(wp), dimension(:), allocatable :: odf_y !< outflow damping factor in y-direction
    real(wp), dimension(:), allocatable :: ptdf_x !< damping factor for potential temperature in
    !< x-direction
    real(wp), dimension(:), allocatable :: ptdf_y !< damping factor for potential temperature in
    !< y-direction
    real(wp), dimension(:), allocatable :: pt_init !< initial profile of potential temperature
    real(wp), dimension(:), allocatable :: q_init !< initial profile of total water mixing ratio
    !< (or total water content with active cloud physics)
    real(wp), dimension(:), allocatable :: rdf !< rayleigh damping factor for velocity components
    real(wp), dimension(:), allocatable :: rdf_sc !< rayleigh damping factor for scalar quantities
    real(wp), dimension(:), allocatable :: ref_state !< reference state of potential temperature
    !< (and density in case of ocean simulation)
    real(wp), dimension(:), allocatable :: rho_air !< air density profile on the uv grid
    real(wp), dimension(:), allocatable :: rho_air_zw !< air density profile on the w grid
    real(wp), dimension(:), allocatable :: s_init !< initial profile of passive scalar concentration
    real(wp), dimension(:), allocatable :: sa_init !< initial profile of salinity (ocean)
    real(wp), dimension(:), allocatable :: scalarflux_input_conversion !< conversion factor array for scalarflux input
    real(wp), dimension(:), allocatable :: scalarflux_output_conversion !< conversion factor array for scalarflux output
    real(wp), dimension(:), allocatable :: ug !< geostrophic wind component in x-direction
    real(wp), dimension(:), allocatable :: u_init !< initial profile of horizontal velocity component u
    real(wp), dimension(:), allocatable :: u_stokes_zu !< u-component of Stokes drift velocity at zu levels
    real(wp), dimension(:), allocatable :: u_stokes_zw !< u-component of Stokes drift velocity at zw levels
    real(wp), dimension(:), allocatable :: vg !< geostrophic wind component in y-direction
    real(wp), dimension(:), allocatable :: v_init !< initial profile of horizontal velocity component v
    real(wp), dimension(:), allocatable :: v_stokes_zu !< v-component of Stokes drift velocity at zu levels
    real(wp), dimension(:), allocatable :: v_stokes_zw !< v-component of Stokes drift velocity at zw levels
    real(wp), dimension(:), allocatable :: waterflux_input_conversion !< conversion factor array for waterflux input
    real(wp), dimension(:), allocatable :: waterflux_output_conversion !< conversion factor array for waterflux output
    real(wp), dimension(:), allocatable :: w_subs !< subsidence/ascent velocity
    real(wp), dimension(:), allocatable :: x !< horizontal grid coordinate of v- and s-grid (in m)
    real(wp), dimension(:), allocatable :: xu !< horizontal grid coordinate of u-grid (in m)
    real(wp), dimension(:), allocatable :: y !< horizontal grid coordinate of u- and s-grid (in m)
    real(wp), dimension(:), allocatable :: yv !< horizontal grid coordinate of v-grid (in m)
    real(wp), dimension(:), allocatable :: zu !< vertical grid coordinate of u-, v-, and s-grid (in m)
    real(wp), dimension(:), allocatable :: zw !< vertical grid coordinate of w-grid (in m)

    real(wp), dimension(:, :), allocatable :: c_u !< phase speed of u-velocity component
    real(wp), dimension(:, :), allocatable :: c_v !< phase speed of v-velocity component
    real(wp), dimension(:, :), allocatable :: c_w !< phase speed of w-velocity component
    real(wp), dimension(:, :), allocatable :: diss_s_diss !< artificial numerical dissipation flux at south face of grid
    !< box - TKE dissipation
    real(wp), dimension(:, :), allocatable :: diss_s_e !< artificial numerical dissipation flux at south face of grid
    !< box - subgrid-scale TKE
    real(wp), dimension(:, :), allocatable :: diss_s_nc !< artificial numerical dissipation flux at south face of grid
    !< box - clouddrop-number concentration
    real(wp), dimension(:, :), allocatable :: diss_s_ng !< artificial numerical dissipation flux at south face of grid
    !< box - graupel number concentration
    real(wp), dimension(:, :), allocatable :: diss_s_ni !< artificial numerical dissipation flux at south face of grid
    !< box - ice crystal-number concentration
    real(wp), dimension(:, :), allocatable :: diss_s_nr !< artificial numerical dissipation flux at south face of grid
    !< box - raindrop-number concentration
    real(wp), dimension(:, :), allocatable :: diss_s_ns !< artificial numerical dissipation flux at south face of grid
    !< box - snow-number concentration
    real(wp), dimension(:, :), allocatable :: diss_s_pt !< artificial numerical dissipation flux at south face of grid
    !< box - potential temperature
    real(wp), dimension(:, :), allocatable :: diss_s_q !< artificial numerical dissipation flux at south face of grid
    !< box - total water mixing ratio
    real(wp), dimension(:, :), allocatable :: diss_s_qc !< artificial numerical dissipation flux at south face of grid
    !< box - cloudwater mixing ratio
    real(wp), dimension(:, :), allocatable :: diss_s_qg !< artificial numerical dissipation flux at south face of grid
    !< box - graupel mixing ratio
    real(wp), dimension(:, :), allocatable :: diss_s_qi !< artificial numerical dissipation flux at south face of grid
    !< box - ice crystal mixing ratio
    real(wp), dimension(:, :), allocatable :: diss_s_qr !< artificial numerical dissipation flux at south face of grid
    !< box - rainwater mixing ratio
    real(wp), dimension(:, :), allocatable :: diss_s_qs !< artificial numerical dissipation flux at south face of grid
    !< box - snow mixing ratio
    real(wp), dimension(:, :), allocatable :: diss_s_s !< artificial numerical dissipation flux at south face of grid
    !< box - passive scalar
    real(wp), dimension(:, :), allocatable :: diss_s_sa !< artificial numerical dissipation flux at south face of grid
    !< box - salinity
    real(wp), dimension(:, :), allocatable :: diss_s_u !< artificial numerical dissipation flux at south face of grid
    !< box - u-component
    real(wp), dimension(:, :), allocatable :: diss_s_v !< artificial numerical dissipation flux at south face of grid
    !< box - v-component
    real(wp), dimension(:, :), allocatable :: diss_s_w !< artificial numerical dissipation flux at south face of grid
    !< box - w-component
    real(wp), dimension(:, :), allocatable :: flux_s_diss !< 6th-order advective flux at south face of grid box -
    !< TKE dissipation
    real(wp), dimension(:, :), allocatable :: flux_s_e !< 6th-order advective flux at south face of grid box -
    !< subgrid-scale TKE
    real(wp), dimension(:, :), allocatable :: flux_s_nc !< 6th-order advective flux at south face of grid box -
    !< clouddrop-number concentration
    real(wp), dimension(:, :), allocatable :: flux_s_ng !< 6th-order advective flux at south face of grid box -
    !< graupel-number concentration
    real(wp), dimension(:, :), allocatable :: flux_s_ni !< 6th-order advective flux at south face of grid box -
    !< icecrystal-number concentration
    real(wp), dimension(:, :), allocatable :: flux_s_nr !< 6th-order advective flux at south face of grid box -
    !< raindrop-number concentration
    real(wp), dimension(:, :), allocatable :: flux_s_ns !< 6th-order advective flux at south face of grid box -
    !< graupel-number concentration
    real(wp), dimension(:, :), allocatable :: flux_s_pt !< 6th-order advective flux at south face of grid box -
    !< potential temperature
    real(wp), dimension(:, :), allocatable :: flux_s_q !< 6th-order advective flux at south face of grid box -
    !< total water mixing ratio
    real(wp), dimension(:, :), allocatable :: flux_s_qc !< 6th-order advective flux at south face of grid box -
    !< cloudwater mixing ratio
    real(wp), dimension(:, :), allocatable :: flux_s_qg !< 6th-order advective flux at south face of grid box -
    !< graupel mixing ratio
    real(wp), dimension(:, :), allocatable :: flux_s_qi !< 6th-order advective flux at south face of grid box -
    !< ice crystal mixing ratio
    real(wp), dimension(:, :), allocatable :: flux_s_qr !< 6th-order advective flux at south face of grid box -
    !< rainwater mixing ratio
    real(wp), dimension(:, :), allocatable :: flux_s_qs !< 6th-order advective flux at south face of grid box -
    !< snow mixing ratio
    real(wp), dimension(:, :), allocatable :: flux_s_s !< 6th-order advective flux at south face of grid box -
    !< passive scalar
    real(wp), dimension(:, :), allocatable :: flux_s_sa !< 6th-order advective flux at south face of grid box -
    !< salinity
    real(wp), dimension(:, :), allocatable :: flux_s_u !< 6th-order advective flux at south face of grid box -
    !< u-component
    real(wp), dimension(:, :), allocatable :: flux_s_v !< 6th-order advective flux at south face of grid box -
    !< v-component
    real(wp), dimension(:, :), allocatable :: flux_s_w !< 6th-order advective flux at south face of grid box -
    !< w-component
    real(wp), dimension(:, :), allocatable :: mean_inflow_profiles !< used for turbulent inflow (non-cyclic boundary conditions)
    real(wp), dimension(:, :), allocatable :: precipitation_amount !< precipitation amount due to gravitational settling
    !< (bulk microphysics)
    real(wp), dimension(:, :), allocatable :: pt_slope_ref !< potential temperature in rotated coordinate system
    !< (in case of sloped surface)
    real(wp), dimension(:, :), allocatable :: total_2d_a !< horizontal array to store the total domain data, used for
    !< atmosphere-ocean coupling (atmosphere data)
    real(wp), dimension(:, :), allocatable :: total_2d_o !< horizontal array to store the total domain data, used for
    !< atmosphere-ocean coupling (ocean data)

    real(wp), dimension(:, :, :), allocatable :: d !< divergence
    real(wp), dimension(:, :, :), allocatable :: de_dx !< gradient of sgs tke in x-direction (lpm)
    real(wp), dimension(:, :, :), allocatable :: de_dy !< gradient of sgs tke in y-direction (lpm)
    real(wp), dimension(:, :, :), allocatable :: de_dz !< gradient of sgs tke in z-direction (lpm)
    real(wp), dimension(:, :, :), allocatable :: diss_l_diss !< artificial numerical dissipation flux at left face of grid box -
    !< TKE dissipation
    real(wp), dimension(:, :, :), allocatable :: diss_l_e !< artificial numerical dissipation flux at left face of grid box -
    !< subgrid-scale TKE
    real(wp), dimension(:, :, :), allocatable :: diss_l_nc !< artificial numerical dissipation flux at left face of grid box -
    !< clouddrop-number concentration
    real(wp), dimension(:, :, :), allocatable :: diss_l_ng !< artificial numerical dissipation flux at left face of grid box -
    !< graupel-number concentration
    real(wp), dimension(:, :, :), allocatable :: diss_l_ni !< artificial numerical dissipation flux at left face of grid box -
    !< ice crystal-number concentration
    real(wp), dimension(:, :, :), allocatable :: diss_l_nr !< artificial numerical dissipation flux at left face of grid box -
    !< raindrop-number concentration
    real(wp), dimension(:, :, :), allocatable :: diss_l_ns !< artificial numerical dissipation flux at left face of grid box -
    !< snow-number concentration
    real(wp), dimension(:, :, :), allocatable :: diss_l_pt !< artificial numerical dissipation flux at left face of grid box -
    !< potential temperature
    real(wp), dimension(:, :, :), allocatable :: diss_l_q !< artificial numerical dissipation flux at left face of grid box -
    !< total water mixing ratio
    real(wp), dimension(:, :, :), allocatable :: diss_l_qc !< artificial numerical dissipation flux at left face of grid box -
    !< cloudwater
    real(wp), dimension(:, :, :), allocatable :: diss_l_qg !< artificial numerical dissipation flux at left face of grid box -
    !< graupel
    real(wp), dimension(:, :, :), allocatable :: diss_l_qi !< artificial numerical dissipation flux at left face of grid box -
    !< ice crystal
    real(wp), dimension(:, :, :), allocatable :: diss_l_qr !< artificial numerical dissipation flux at left face of grid box -
    !< rainwater
    real(wp), dimension(:, :, :), allocatable :: diss_l_qs !< artificial numerical dissipation flux at left face of grid box -
    !< snow
    real(wp), dimension(:, :, :), allocatable :: diss_l_s !< artificial numerical dissipation flux at left face of grid box -
    !< passive scalar
    real(wp), dimension(:, :, :), allocatable :: diss_l_sa !< artificial numerical dissipation flux at left face of grid box -
    !< salinity
    real(wp), dimension(:, :, :), allocatable :: diss_l_u !< artificial numerical dissipation flux at left face of grid box -
    !< u-component
    real(wp), dimension(:, :, :), allocatable :: diss_l_v !< artificial numerical dissipation flux at left face of grid box -
    !< v-component
    real(wp), dimension(:, :, :), allocatable :: diss_l_w !< artificial numerical dissipation flux at left face of grid box -
    !< w-component
    real(wp), dimension(:, :, :), allocatable :: flux_l_diss !< 6th-order advective flux at south face of grid box - TKE dissipation
    real(wp), dimension(:, :, :), allocatable :: flux_l_e !< 6th-order advective flux at south face of grid box - subgrid-scale
    !< TKE
    real(wp), dimension(:, :, :), allocatable :: flux_l_nc !< 6th-order advective flux at south face of grid box - clouddrop-number
    !< concentration
    real(wp), dimension(:, :, :), allocatable :: flux_l_ng !< 6th-order advective flux at south face of grid box -
    !< graupel-number concentration
    real(wp), dimension(:, :, :), allocatable :: flux_l_ni !< 6th-order advective flux at south face of grid box -
    !< ice crystal-number concentration
    real(wp), dimension(:, :, :), allocatable :: flux_l_nr !< 6th-order advective flux at south face of grid box - raindrop-number
    !< concentration
    real(wp), dimension(:, :, :), allocatable :: flux_l_ns
    !< 6th-order advective flux at south face of grid box - snow-number concentration
    real(wp), dimension(:, :, :), allocatable :: flux_l_pt !< 6th-order advective flux at south face of grid box - potential
    !< temperature
    real(wp), dimension(:, :, :), allocatable :: flux_l_q !< 6th-order advective flux at south face of grid box - mixing ratio
    real(wp), dimension(:, :, :), allocatable :: flux_l_qc !< 6th-order advective flux at south face of grid box - cloudwater
    real(wp), dimension(:, :, :), allocatable :: flux_l_qg !< 6th-order advective flux at south face of grid box - graupel
    real(wp), dimension(:, :, :), allocatable :: flux_l_qi !< 6th-order advective flux at south face of grid box - ice crystal
    real(wp), dimension(:, :, :), allocatable :: flux_l_qr !< 6th-order advective flux at south face of grid box - rainwater
    real(wp), dimension(:, :, :), allocatable :: flux_l_qs !< 6th-order advective flux at south face of grid box - snow
    real(wp), dimension(:, :, :), allocatable :: flux_l_s !< 6th-order advective flux at south face of grid box - passive scalar
    real(wp), dimension(:, :, :), allocatable :: flux_l_sa !< 6th-order advective flux at south face of grid box - salinity
    real(wp), dimension(:, :, :), allocatable :: flux_l_u !< 6th-order advective flux at south face of grid box - u-component
    real(wp), dimension(:, :, :), allocatable :: flux_l_v !< 6th-order advective flux at south face of grid box - v-component
    real(wp), dimension(:, :, :), allocatable :: flux_l_w !< 6th-order advective flux at south face of grid box - w-component
    real(wp), dimension(:, :, :), allocatable, target :: kh !< eddy diffusivity for heat
    real(wp), dimension(:, :, :), allocatable, target :: km !< eddy diffusivity for momentum
    real(wp), dimension(:, :, :), allocatable :: prr !< total precipitation rate (all phases)
    real(wp), dimension(:, :, :), allocatable :: prr_cloud !< cloud sedimentation rate
    real(wp), dimension(:, :, :), allocatable :: prr_graupel !< graupel precipitation rate
    real(wp), dimension(:, :, :), allocatable :: prr_ice !< ice precipitation rate
    real(wp), dimension(:, :, :), allocatable :: prr_rain !< rain precipitation rate
    real(wp), dimension(:, :, :), allocatable :: prr_snow !< snow precipitation rate
    real(wp), dimension(:, :, :), allocatable :: p_loc !< local array in multigrid/sor solver containing the pressure which is
    !< iteratively advanced in each iteration step
    real(wp), dimension(:, :, :), allocatable :: tend !< tendency field (time integration)

    real(wp), dimension(:, :, :), allocatable, target :: diss_1 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: diss_2 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: diss_3 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: e_1 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: e_2 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: e_3 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: p !< pointer: perturbation pressure
    real(wp), dimension(:, :, :), allocatable, target :: prho_1 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: nc_1 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: nc_2 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: nc_3 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: ng_1 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: ng_2 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: ng_3 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: ni_1 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: ni_2 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: ni_3 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: nr_1 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: nr_2 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: nr_3 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: ns_1 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: ns_2 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: ns_3 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: pt_1 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: pt_2 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: pt_3 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: q_1 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: q_2 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: q_3 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: qc_1 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: qc_2 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: qc_3 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: qf_1 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: qg_1 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: qg_2 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: qg_3 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: qi_1 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: qi_2 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: qi_3 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: ql_v !< pointer: volume of liquid water
    real(wp), dimension(:, :, :), allocatable, target :: ql_vp !< pointer: liquid water weighting factor
    real(wp), dimension(:, :, :), allocatable, target :: ql_1 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: ql_2 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: qr_1 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: qr_2 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: qr_3 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: qs_1 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: qs_2 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: qs_3 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: rho_1 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: s_1 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: s_2 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: s_3 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: sa_1 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: sa_2 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: sa_3 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: u_1 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: u_2 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: u_3 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: v_1 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: v_2 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: v_3 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: vpt_1 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: w_1 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: w_2 !< pointer for swapping of timelevels for respective quantity
    real(wp), dimension(:, :, :), allocatable, target :: w_3 !< pointer for swapping of timelevels for respective quantity

    real(wp), dimension(:, :, :), pointer, contiguous :: diss !< pointer: TKE dissipation
    real(wp), dimension(:, :, :), pointer, contiguous :: diss_p !< pointer: prognostic value of TKE dissipation
    real(wp), dimension(:, :, :), pointer, contiguous :: e !< pointer: subgrid-scale turbulence kinetic energy (sgs tke)
    real(wp), dimension(:, :, :), pointer, contiguous :: e_p !< pointer: prognostic value of sgs tke
    real(wp), dimension(:, :, :), pointer, contiguous :: nc !< pointer: cloud drop number density
    real(wp), dimension(:, :, :), pointer, contiguous :: nc_p !< pointer: prognostic value of cloud drop number density
    real(wp), dimension(:, :, :), pointer, contiguous :: ng !< pointer: graupel number density
    real(wp), dimension(:, :, :), pointer, contiguous :: ng_p !< pointer: prognostic value of graupel number density
    real(wp), dimension(:, :, :), pointer, contiguous :: ni !< pointer: ice crystal number density
    real(wp), dimension(:, :, :), pointer, contiguous :: ni_p !< pointer: prognostic value of ice crystal number density
    real(wp), dimension(:, :, :), pointer, contiguous :: nr !< pointer: rain drop number density
    real(wp), dimension(:, :, :), pointer, contiguous :: nr_p !< pointer: prognostic value of rain drop number density
    real(wp), dimension(:, :, :), pointer, contiguous :: ns !< pointer: snow number density
    real(wp), dimension(:, :, :), pointer, contiguous :: ns_p !< pointer: prognostic value of snow number density
    real(wp), dimension(:, :, :), pointer, contiguous :: prho !< pointer: potential density
    real(wp), dimension(:, :, :), pointer, contiguous :: pt !< pointer: potential temperature
    real(wp), dimension(:, :, :), pointer, contiguous :: pt_p !< pointer: prognostic value of potential temperature
    real(wp), dimension(:, :, :), pointer, contiguous :: q !< pointer: mixing ratio
    real(wp), dimension(:, :, :), pointer, contiguous :: q_p !< pointer: prognostic value of mixing ratio
    real(wp), dimension(:, :, :), pointer, contiguous :: qc !< pointer: cloud water content
    real(wp), dimension(:, :, :), pointer, contiguous :: qc_p !< pointer: prognostic value cloud water content
    real(wp), dimension(:, :, :), pointer, contiguous :: qf !< pointer: frozen water content
    real(wp), dimension(:, :, :), pointer, contiguous :: qg !< pointer: graupel water content
    real(wp), dimension(:, :, :), pointer, contiguous :: qg_p !< pointer: prognostic value graupel water content
    real(wp), dimension(:, :, :), pointer, contiguous :: qi !< pointer: ice crystal content
    real(wp), dimension(:, :, :), pointer, contiguous :: qi_p !< pointer: prognostic value ice crystal content
    real(wp), dimension(:, :, :), pointer, contiguous :: ql !< pointer: liquid water content
    real(wp), dimension(:, :, :), pointer, contiguous :: ql_c !< pointer: change in liquid water content due to
    !< condensation/evaporation during last time step
    real(wp), dimension(:, :, :), pointer, contiguous :: qr !< pointer: rain water content
    real(wp), dimension(:, :, :), pointer, contiguous :: qr_p !< pointer: prognostic value of rain water content
    real(wp), dimension(:, :, :), pointer, contiguous :: qs !< pointer: rain water content
    real(wp), dimension(:, :, :), pointer, contiguous :: qs_p !< pointer: prognostic value of rain water content
    real(wp), dimension(:, :, :), pointer, contiguous :: rho_ocean !< pointer: density of ocean
    real(wp), dimension(:, :, :), pointer, contiguous :: s !< pointer: passive scalar
    real(wp), dimension(:, :, :), pointer, contiguous :: s_p !< pointer: prognostic value of passive scalar
    real(wp), dimension(:, :, :), pointer, contiguous :: sa !< pointer: ocean salinity
    real(wp), dimension(:, :, :), pointer, contiguous :: sa_p !< pointer: prognostic value of ocean salinity
    real(wp), dimension(:, :, :), pointer, contiguous :: tdiss_m !< pointer: weighted tendency of diss for previous sub-timestep
    !< (Runge-Kutta)
    real(wp), dimension(:, :, :), pointer, contiguous :: te_m !< pointer: weighted tendency of e for previous sub-timestep
    !< (Runge-Kutta)
    real(wp), dimension(:, :, :), pointer, contiguous :: tnc_m !< pointer: weighted tendency of nc for previous sub-timestep
    !< (Runge-Kutta)
    real(wp), dimension(:, :, :), pointer, contiguous :: tng_m
    !< pointer: weighted tendency of ng for previous sub-timestep (Runge-Kutta)
    real(wp), dimension(:, :, :), pointer, contiguous :: tni_m !< pointer: weighted tendency of ni for previous sub-timestep
    !< (Runge-Kutta)
    real(wp), dimension(:, :, :), pointer, contiguous :: tnr_m !< pointer: weighted tendency of nr for previous sub-timestep
    !< (Runge-Kutta)
    real(wp), dimension(:, :, :), pointer, contiguous :: tns_m
    !< pointer: weighted tendency of ns for previous sub-timestep (Runge-Kutta)
    real(wp), dimension(:, :, :), pointer, contiguous :: tpt_m !< pointer: weighted tendency of pt for previous sub-timestep
    !< (Runge-Kutta)
    real(wp), dimension(:, :, :), pointer, contiguous :: tq_m !< pointer: weighted tendency of q for previous sub-timestep
    !< (Runge-Kutta)
    real(wp), dimension(:, :, :), pointer, contiguous :: tqc_m !< pointer: weighted tendency of qc for previous sub-timestep
    !< (Runge-Kutta)
    real(wp), dimension(:, :, :), pointer, contiguous :: tqg_m
    !< pointer: weighted tendency of qg for previous sub-timestep (Runge-Kutta)
    real(wp), dimension(:, :, :), pointer, contiguous :: tqi_m !< pointer: weighted tendency of qi for previous sub-timestep
    !< (Runge-Kutta)
    real(wp), dimension(:, :, :), pointer, contiguous :: tqr_m !< pointer: weighted tendency of qr for previous sub-timestep
    !< (Runge-Kutta)
    real(wp), dimension(:, :, :), pointer, contiguous :: tqs_m
    !< pointer: weighted tendency of qs for previous sub-timestep (Runge-Kutta)
    real(wp), dimension(:, :, :), pointer, contiguous :: ts_m !< pointer: weighted tendency of s for previous sub-timestep
    !< (Runge-Kutta)
    real(wp), dimension(:, :, :), pointer, contiguous :: tsa_m !< pointer: weighted tendency of sa for previous sub-timestep
    !< (Runge-Kutta)
    real(wp), dimension(:, :, :), pointer, contiguous :: tu_m !< pointer: weighted tendency of u for previous sub-timestep
    !< (Runge-Kutta)
    real(wp), dimension(:, :, :), pointer, contiguous :: tv_m !< pointer: weighted tendency of v for previous sub-timestep
    !< (Runge-Kutta)
    real(wp), dimension(:, :, :), pointer, contiguous :: tw_m !< pointer: weighted tendency of w for previous sub-timestep
    !< (Runge-Kutta)
    real(wp), dimension(:, :, :), pointer, contiguous :: u !< pointer: horizontal velocity component u (x-direction)
    real(wp), dimension(:, :, :), pointer, contiguous :: u_p !< pointer: prognostic value of u
    real(wp), dimension(:, :, :), pointer, contiguous :: v !< pointer: horizontal velocity component v (y-direction)
    real(wp), dimension(:, :, :), pointer, contiguous :: v_p !< pointer: prognostic value of v
    real(wp), dimension(:, :, :), pointer, contiguous :: vpt !< pointer: virtual potential temperature
    real(wp), dimension(:, :, :), pointer, contiguous :: w !< pointer: vertical velocity component w (z-direction)
    real(wp), dimension(:, :, :), pointer, contiguous :: w_p !< pointer: prognostic value of w

    save

end module arrays_3d