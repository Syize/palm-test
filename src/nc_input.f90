!> @file netcdf_data_input_mod.f90
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
! Copyright 2022-2022 pecanode GmbH
!--------------------------------------------------------------------------------------------------!
! Authors:
! --------
! @author Matthias Suehring
! @author Edward C. Chan
! @author Emanuele Russo
! Description:
! ------------
!> Modulue contains routines to input data according to Palm input data
!> standart using dynamic and static input files.
!> @todo - Chemistry: revise reading of netcdf file and ajdust formatting according to standard!!!
!> @todo - Order input alphabetically
!> @todo - Revise error messages and error numbers
!> @todo - Input of missing quantities (chemical species, emission rates)
!> @todo - Definition and input of still missing variable attributes
!> @todo - Input of initial geostrophic wind profiles with cyclic conditions.
!> @todo - remove z dimension from default_emission_data nad preproc_emission_data and correpsonding
!>         subroutines get_var_5d_real and get_var_5d_dynamic
!> @todo - decpreciate chem_emis_att_type@nspec
!> @todo - depreciate subroutines get_variable_4d_to_3d_real and get_variable_5d_to_4d_real
!> @todo - introduce useful debug_message(s)
!--------------------------------------------------------------------------------------------------!
module netcdf_data_input_mod


    ! #if defined( __parallel )
    !     use MPI
    ! #endif
    use boundary_settings_mod, &
        only: set_lateral_neumann_bc

    use control_parameters, &
        only: air_chemistry, &
              bc_dirichlet_l, &
              bc_dirichlet_n, &
              bc_dirichlet_r, &
              bc_dirichlet_s, &
              bc_lr_cyc, &
              bc_ns_cyc, &
              bc_radiation_l, &
              bc_radiation_n, &
              bc_radiation_r, &
              bc_radiation_s, &
              coupling_char, &
              cut_cell_topography, &
              dz, &
              humidity, &
              initializing_actions, &
              io_blocks, &
              io_group, &
              land_surface, &
              message_string, &
              neutral, &
              number_dz, &
              pt_surface, &
              q_surface, &
              topo_no_distinct, &
              urban_surface

    ! use cpulog, &
    !     only: cpu_log, &
    !           log_point_s

    use exchange_horiz_mod, &
        only: exchange_horiz, &
              exchange_horiz_2d, &
              exchange_horiz_2d_byte

    use, intrinsic :: ieee_arithmetic

    use kinds

    use indices, &
        only: nbgp, &
              nx, &
              nxl, &
              nxlg, &
              nxlu, &
              nxr, &
              nxrg, &
              ny, &
              nyn, &
              nyng, &
              nys, &
              nysg, &
              nysv, &
              nz, &
              nzb, &
              nzt, &
              topo_flags

    ! use kinds

    ! #if defined ( __netcdf )
    use netcdf

    ! #endif
    use pegrid

    implicit none

    ! integer, parameter :: iwp = int32
    ! real, parameter :: wp = real32
    ! integer, parameter :: ibp = int8

    !-- Define type for dimensions.
    type dims_xy
        integer(iwp) :: nx !< dimension length in x
        integer(iwp) :: ny !< dimension length in y
        integer(iwp) :: nz !< dimension length in z
        real(wp), dimension(:), allocatable :: x !< dimension array in x
        real(wp), dimension(:), allocatable :: y !< dimension array in y
        real(wp), dimension(:), allocatable :: z !< dimension array in z
    end type dims_xy

    type init_type

        character(LEN=16) :: init_char = 'init_atmosphere_' !< leading substring for init variables
        character(LEN=23) :: origin_time = '2000-01-01 00:00:00 +00' !< reference time of input data

        character(LEN=100), dimension(:), allocatable :: var_names_chem !< list of chemistry variable names that can potentially be
        !< on file

        integer(iwp) :: lod_msoil !< level of detail - soil moisture
        integer(iwp) :: lod_pt !< level of detail - pt
        integer(iwp) :: lod_q !< level of detail - q
        integer(iwp) :: lod_tsoil !< level of detail - soil temperature
        integer(iwp) :: lod_u !< level of detail - u-component
        integer(iwp) :: lod_v !< level of detail - v-component
        integer(iwp) :: lod_w !< level of detail - w-component
        integer(iwp) :: nx !< number of scalar grid points along x in dynamic input file
        integer(iwp) :: nxu !< number of u grid points along x in dynamic input file
        integer(iwp) :: ny !< number of scalar grid points along y in dynamic input file
        integer(iwp) :: nyv !< number of v grid points along y in dynamic input file
        integer(iwp) :: nzs !< number of vertical soil levels in dynamic input file
        integer(iwp) :: nzu !< number of vertical levels on scalar grid in dynamic input file
        integer(iwp) :: nzw !< number of vertical levels on w grid in dynamic input file

        integer(iwp), dimension(:), allocatable :: lod_chem !< level of detail - chemistry variables

        logical :: from_file_msoil = .false. !< flag indicating whether soil moisture is already initialized from file
        logical :: from_file_pt = .false. !< flag indicating whether pt is already initialized from file
        logical :: from_file_q = .false. !< flag indicating whether q is already initialized from file
        logical :: from_file_tsoil = .false. !< flag indicating whether soil temperature is already initialized from file
        logical :: from_file_u = .false. !< flag indicating whether u is already initialized from file
        logical :: from_file_ug = .false. !< flag indicating whether ug is already initialized from file
        logical :: from_file_v = .false. !< flag indicating whether v is already initialized from file
        logical :: from_file_vg = .false. !< flag indicating whether ug is already initialized from file
        logical :: from_file_w = .false. !< flag indicating whether w is already initialized from file

        logical, dimension(:), allocatable :: from_file_chem !< flag indicating whether chemistry variable is read from file

        real(wp) :: fill_msoil !< fill value for soil moisture
        real(wp) :: fill_pt !< fill value for pt
        real(wp) :: fill_q !< fill value for q
        real(wp) :: fill_tsoil !< fill value for soil temperature
        real(wp) :: fill_u !< fill value for u
        real(wp) :: fill_v !< fill value for v
        real(wp) :: fill_w !< fill value for w
        real(wp) :: latitude = 0.0_wp !< latitude of the lower left corner
        real(wp) :: longitude = 0.0_wp !< longitude of the lower left corner
        real(wp) :: origin_x = 500000.0_wp !< UTM easting of the lower left corner
        real(wp) :: origin_y = 0.0_wp !< UTM northing of the lower left corner
        real(wp) :: origin_z = 0.0_wp !< reference height of input data
        real(wp) :: rotation_angle = 0.0_wp !< rotation angle of input data

        real(wp), dimension(:), allocatable :: fill_chem !< fill value - chemistry variables
        real(wp), dimension(:), allocatable :: msoil_1d !< initial vertical profile of soil moisture
        real(wp), dimension(:), allocatable :: pt_init !< initial vertical profile of pt
        real(wp), dimension(:), allocatable :: q_init !< initial vertical profile of q
        real(wp), dimension(:), allocatable :: tsoil_1d !< initial vertical profile of soil temperature
        real(wp), dimension(:), allocatable :: u_init !< initial vertical profile of u
        real(wp), dimension(:), allocatable :: ug_init !< initial vertical profile of ug
        real(wp), dimension(:), allocatable :: v_init !< initial vertical profile of v
        real(wp), dimension(:), allocatable :: vg_init !< initial vertical profile of ug
        real(wp), dimension(:), allocatable :: w_init !< initial vertical profile of w
        real(wp), dimension(:), allocatable :: z_soil !< vertical levels in soil in dynamic input file, used for interpolation
        real(wp), dimension(:), allocatable :: zu_atmos !< vertical levels at scalar grid in dynamic input file, used for
        !< interpolation
        real(wp), dimension(:), allocatable :: zw_atmos !< vertical levels at w grid in dynamic input file, used for
        !< interpolation

        real(wp), dimension(:, :), allocatable :: chem_init !< initial vertical profiles of chemistry variables

        real(wp), dimension(:, :, :), allocatable :: msoil_3d !< initial 3d soil moisture provided by dynamic file and
        !< interpolated onto soil grid
        real(wp), dimension(:, :, :), allocatable :: tsoil_3d !< initial 3d soil temperature provided by dynamic file and
        !< interpolated onto soil grid

    end type init_type
    !-- Data type for the general information of chemistry emissions, do not dependent on the particular chemical species
    type chem_emis_att_type
        !--    DIMENSIONS
        integer(iwp) :: nspec = 0 !< no of chem species provided in emission_values
        integer(iwp) :: n_emiss_species = 0 !< no of chem species provided in emission_values
        !< same function as nspec, which will be depreciated
        integer(iwp) :: ncat = 0 !< number of emission categories
        integer(iwp) :: nvoc = 0 !< number of VOC components
        integer(iwp) :: npm = 0 !< number of PM components
        integer(iwp) :: nnox = 2 !< number of NOx components: NO and NO2
        integer(iwp) :: nsox = 2 !< number of SOX components: SO and SO4
        integer(iwp) :: nhoursyear !< number of hours of a specific year in the HOURLY mode
        !< of the default mode
        integer(iwp) :: nmonthdayhour !< number of month days and hours in the MDH mode
        !< of the default mode
        integer(iwp) :: dt_emission !< Number of emissions timesteps for one year
        !< in the pre-processed emissions case
        !--    1d emission input variables
        character(LEN=25), allocatable, dimension(:) :: pm_name !< Names of PM components
        character(LEN=25), allocatable, dimension(:) :: cat_name !< Emission category names
        character(LEN=25), allocatable, dimension(:) :: species_name !< Names of emission chemical species
        character(LEN=25), allocatable, dimension(:) :: voc_name !< Names of VOCs components
        character(LEN=25) :: units !< Units

        integer(iwp) :: i_hour !< indices for assigning emission values at different
        !< timesteps
        integer(iwp), allocatable, dimension(:) :: cat_index !< Indices for emission categories
        integer(iwp), allocatable, dimension(:) :: species_index !< Indices for emission chem species

        real(wp), allocatable, dimension(:) :: xm !< Molecular masses of emission chem species

        !--    2d emission input variables
        real(wp), allocatable, dimension(:, :) :: hourly_emis_time_factor !< Time factors for HOURLY emissions (DEFAULT mode)
        real(wp), allocatable, dimension(:, :) :: mdh_emis_time_factor !< Time factors for MDH emissions (DEFAULT mode)
        real(wp), allocatable, dimension(:, :) :: nox_comp !< Composition of NO and NO2
        real(wp), allocatable, dimension(:, :) :: sox_comp !< Composition of SO2 and SO4
        real(wp), allocatable, dimension(:, :) :: voc_comp !< Composition of VOC components (not fixed)

        !--    3d emission input variables
        real(wp), allocatable, dimension(:, :, :) :: pm_comp !< Composition of PM components (not fixed)

    end type chem_emis_att_type

    !-- Data type for the values of chemistry emissions
    type chem_emis_val_type

        !REAL(wp),ALLOCATABLE, DIMENSION(:,:)     :: stack_height           !< stack height
        real(wp), allocatable, dimension(:, :, :) :: default_emission_data !< Emission input values for LOD1 (DEFAULT mode)
        real(wp), allocatable, dimension(:, :, :, :) :: preproc_emission_data !< Emission input values for LOD2 (PRE-PROCESSED mode)

    end type chem_emis_val_type

    !-- Define data structures for different input data types.
    !-- 8-bit Integer 2D
    type int_2d_8bit
        integer(ibp) :: fill = -127 !< fill value
        integer(ibp), dimension(:, :), allocatable :: var !< respective variable

        logical :: from_file = .false. !< flag indicating whether an input variable is available and read from file or default
        !< values are used
    end type int_2d_8bit
    !-- 8-bit Integer 3D
    type int_3d_8bit
        integer(ibp) :: fill = -127 !< fill value
        integer(ibp), dimension(:, :, :), allocatable :: var_3d !< respective variable

        logical :: from_file = .false. !< flag indicating whether an input variable is available and read from file or default
        !< values are used
    end type int_3d_8bit
    !-- 32-bit Integer 2D
    type int_2d_32bit
        integer(iwp) :: fill = -9999 !< fill value
        integer(iwp), dimension(:, :), allocatable :: var !< respective variable

        logical :: from_file = .false. !< flag indicating whether an input variable is available and read from file or default
        !< values are used
    end type int_2d_32bit
    !-- Define data type to read 1D or 3D real variables.
    type real_1d_3d
        logical :: from_file = .false. !< flag indicating whether an input variable is available and read from file or default
        !< values are used

        integer(iwp) :: lod = -1 !< level-of-detail

        real(wp) :: fill = -9999.9_wp !< fill value

        real(wp), dimension(:), allocatable :: var1d !< respective 1D variable
        real(wp), dimension(:, :, :), allocatable :: var3d !< respective 3D variable
    end type real_1d_3d
    !-- Define data type to read 2D real variables
    type real_2d
        logical :: from_file = .false. !< flag indicating whether an input variable is available and read from file or default
        !< values are used

        integer(iwp) :: lod !< level-of-detail

        real(wp) :: fill = -9999.9_wp !< fill value
        real(wp), dimension(:, :), allocatable :: var !< respective variable
    end type real_2d

    !-- Define data type to read 3D real variables
    type real_3d
        logical :: from_file = .false. !< flag indicating whether an input variable is available and read from file or default
        !< values are used

        integer(iwp) :: nz !< number of grid points along vertical dimension

        real(wp) :: fill = -9999.9_wp !< fill value
        real(wp), dimension(:, :, :), allocatable :: var !< respective variable
    end type real_3d
    !-- Define data structure where the dimension and type of the input depends on the given level of
    !-- detail.
    !-- For buildings, the input is either 2D float, or 3d byte.
    type build_in
        integer(iwp) :: lod = 1 !< level of detail
        integer(ibp) :: fill2 = -127 !< fill value for lod = 2
        integer(iwp) :: nz !< number of vertical layers in file
        integer(ibp), dimension(:, :, :), allocatable :: var_3d !< 3d variable (lod = 2)

        real(wp), dimension(:), allocatable :: z !< vertical coordinate for 3D building, used for consistency check

        logical :: from_file = .false. !< flag indicating whether an input variable is available and read from file or default
        !< values are used

        real(wp) :: fill1 = -9999.9_wp !< fill values for lod = 1
        real(wp), dimension(:, :), allocatable :: var_2d !< 2d variable (lod = 1)
        real(wp), dimension(:, :), allocatable :: oro_max !< terraing height under particular buildings
    end type build_in

    !-- For soil_type, the input is either 2D or 3D one-byte integer.
    type soil_in
        integer(iwp) :: lod = 1 !< level of detail
        integer(ibp) :: fill = -127 !< fill value for lod = 2
        integer(ibp), dimension(:, :), allocatable :: var_2d !< 2d variable (lod = 1)
        integer(ibp), dimension(:, :, :), allocatable :: var_3d !< 3d variable (lod = 2)

        logical :: from_file = .false. !< flag indicating whether an input variable is available and read from file or default
        !< values are used
    end type soil_in

    !-- Define data type for fractions between surface types
    type fracs
        integer(iwp) :: nf !< total number of fractions
        integer(iwp), dimension(:), allocatable :: nfracs !< dimension array for fraction

        logical :: from_file = .false. !< flag indicating whether an input variable is available and read from file or default
        !< values are used

        real(wp) :: fill = -9999.9_wp !< fill value
        real(wp), dimension(:, :, :), allocatable :: frac !< respective fraction between different surface types
    end type fracs
    !-- Data type for parameter lists, Depending on the given level of detail, the input is 3D or 4D
    type pars
        integer(iwp) :: lod = 1 !< level of detail
        integer(iwp) :: np !< total number of parameters
        integer(iwp) :: nz !< vertical dimension - number of soil layers

        logical :: from_file = .false. !< flag indicating whether an input variable is available and read from file or default
        !< values are used

        real(wp) :: fill = -9999.9_wp !< fill value
        real(wp), dimension(:, :, :), allocatable :: pars_xy !< respective parameters, level of detail = 1
        real(wp), dimension(:, :, :, :), allocatable :: pars_xyz !< respective parameters, level of detail = 2
    end type pars
    !-- Data type for surface parameter lists
    type pars_surf
        integer(iwp) :: np !< total number of parameters
        integer(iwp) :: nsurf !< number of local surfaces
        integer(iwp), dimension(:, :, :), allocatable :: index_ji !< index for beginning and end of surfaces at (j,i)
        integer(iwp), dimension(:, :), allocatable :: coords !< (k,j,i,norm_z,norm_y,norm_x)
        !< k,j,i:                surface position
        !< norm_z,norm_y,norm_x: surface normal vector

        logical :: from_file = .false. !< flag indicating whether an input variable is available and read from file or default
        !< values are used

        real(wp) :: fill = -9999.9_wp !< fill value
        real(wp), dimension(:, :), allocatable :: pars !< respective parameters per surface
    end type pars_surf
    !-- Define type for global file attributes
    !-- Please refer to the PALM data standard for a detailed description of each attribute.
    type global_atts_type
        character(LEN=200) :: acronym = ' ' !< acronym of institution
        character(LEN=7) :: acronym_char = 'acronym' !< name of attribute
        character(LEN=200) :: author = ' ' !< first name, last name, email adress
        character(LEN=6) :: author_char = 'author' !< name of attribute
        character(LEN=200) :: campaign = 'PALM-4U' !< name of campaign
        character(LEN=8) :: campaign_char = 'campaign' !< name of attribute
        character(LEN=200) :: comment = ' ' !< comment to data
        character(LEN=7) :: comment_char = 'comment' !< name of attribute
        character(LEN=200) :: contact_person = ' ' !< first name, last name, email adress
        character(LEN=14) :: contact_person_char = 'contact_person' !< name of attribute
        character(LEN=200) :: conventions = 'CF-1.7' !< netCDF convention
        character(LEN=11) :: conventions_char = 'Conventions' !< name of attribute
        character(LEN=23) :: creation_time = ' ' !< creation time of data set
        character(LEN=13) :: creation_time_char = 'creation_time' !< name of attribute
        character(LEN=200) :: data_content = 'airmeteo' !< content of data set
        character(LEN=12) :: data_content_char = 'data_content' !< name of attribute
        character(LEN=200) :: dependencies = ' ' !< dependencies of data set
        character(LEN=12) :: dependencies_char = 'dependencies' !< name of attribute
        character(LEN=200) :: history = ' ' !< information about data processing
        character(LEN=7) :: history_char = 'history' !< name of attribute
        character(LEN=200) :: institution = ' ' !< name of responsible institution
        character(LEN=11) :: institution_char = 'institution' !< name of attribute
        character(LEN=200) :: keywords = ' ' !< keywords of data set
        character(LEN=8) :: keywords_char = 'keywords' !< name of attribute
        character(LEN=200) :: licence = ' ' !< licence of data set
        character(LEN=7) :: licence_char = 'licence' !< name of attribute
        character(LEN=200) :: location = ' ' !< place which refers to data set
        character(LEN=8) :: location_char = 'location' !< name of attribute
        character(LEN=10) :: origin_lat_char = 'origin_lat' !< name of attribute
        character(LEN=10) :: origin_lon_char = 'origin_lon' !< name of attribute
        character(LEN=23) :: origin_time = '2000-01-01 00:00:00 +00' !< reference time
        character(LEN=11) :: origin_time_char = 'origin_time' !< name of attribute
        character(LEN=8) :: origin_x_char = 'origin_x' !< name of attribute
        character(LEN=8) :: origin_y_char = 'origin_y' !< name of attribute
        character(LEN=8) :: origin_z_char = 'origin_z' !< name of attribute
        character(LEN=12) :: palm_version_char = 'palm_version' !< name of attribute
        character(LEN=200) :: references = ' ' !< literature referring to data set
        character(LEN=10) :: references_char = 'references' !< name of attribute
        character(LEN=14) :: rotation_angle_char = 'rotation_angle' !< name of attribute
        character(LEN=200) :: site = ' ' !< name of model domain
        character(LEN=4) :: site_char = 'site' !< name of attribute
        character(LEN=200) :: source = ' ' !< source of data set
        character(LEN=6) :: source_char = 'source' !< name of attribute
        character(LEN=200) :: title = ' ' !< title of data set
        character(LEN=5) :: title_char = 'title' !< name of attribute
        character(LEN=7) :: version_char = 'version' !< name of attribute

        integer(iwp) :: version !< version of data set

        real(wp) :: fillvalue = -9999.0 !< default fill value
        real(wp) :: origin_lat !< latitude of lower left corner
        real(wp) :: origin_lon !< longitude of lower left corner
        real(wp) :: origin_x !< easting (UTM coordinate) of lower left corner
        real(wp) :: origin_y !< northing (UTM coordinate) of lower left corner
        real(wp) :: origin_z !< reference height
        real(wp) :: palm_version !< PALM version of data set
        real(wp) :: rotation_angle !< rotation angle of coordinate system of data set
    end type global_atts_type
    !-- Define type for coordinate reference system (crs)
    type crs_type
        character(LEN=200) :: epsg_code = 'EPSG:25831' !< EPSG code
        character(LEN=200) :: grid_mapping_name = 'transverse_mercator' !< name of grid mapping
        character(LEN=200) :: long_name = 'coordinate reference system' !< name of variable crs
        character(LEN=200) :: units = 'm' !< unit of crs

        real(wp) :: false_easting = 500000.0_wp !< false easting
        real(wp) :: false_northing = 0.0_wp !< false northing
        real(wp) :: inverse_flattening = 298.257223563_wp !< 1/f (default for WGS84)
        real(wp) :: latitude_of_projection_origin = 0.0_wp !< latitude of projection origin
        real(wp) :: longitude_of_central_meridian = 3.0_wp !< longitude of central meridian of UTM zone (default: zone 31)
        real(wp) :: longitude_of_prime_meridian = 0.0_wp !< longitude of prime meridian
        real(wp) :: scale_factor_at_central_meridian = 0.9996_wp !< scale factor of UTM coordinates
        real(wp) :: semi_major_axis = 6378137.0_wp !< length of semi major axis (default for WGS84)
    end type crs_type

    !-- Define variables
    type(crs_type) :: coord_ref_sys !< coordinate reference system

    type(init_type) :: init_3d !< data structure for the initialization of the 3D flow and soil fields
    type(init_type) :: init_model !< data structure for the initialization of the model

    !-- Define 2D variables of type NC_BYTE
    type(int_2d_8bit) :: albedo_type_f !< input variable for albedo type
    type(int_2d_8bit) :: building_type_f !< input variable for building type
    type(int_2d_8bit) :: pavement_type_f !< input variable for pavenment type
    type(int_2d_8bit) :: street_crossing_f !< input variable for water type
    type(int_2d_8bit) :: street_type_f !< input variable for water type
    type(int_2d_8bit) :: vegetation_type_f !< input variable for vegetation type
    type(int_2d_8bit) :: water_type_f !< input variable for water type
    !-- Define 2D variables of type NC_INT
    type(int_2d_32bit) :: building_id_f !< input variable for building ID
    !-- Define 2D variables of type NC_FLOAT
    type(real_2d) :: terrain_height_f !< input variable for terrain height
    !-- Define 3D variables of type NC_FLOAT
    type(real_3d) :: root_area_density_lsm_f !< input variable for root area density - parametrized vegetation
    !-- Define input variable for buildings
    type(build_in) :: buildings_f !< input variable for buildings
    !-- Define input variables for soil_type
    type(soil_in) :: soil_type_f !< input variable for soil type

    type(fracs) :: surface_fraction_f !< input variable for surface fraction

    type(pars) :: albedo_pars_f !< input variable for albedo parameters
    type(pars) :: pavement_pars_f !< input variable for pavement parameters
    type(pars) :: pavement_subsurface_pars_f !< input variable for pavement parameters
    type(pars) :: soil_pars_f !< input variable for soil parameters
    type(pars) :: vegetation_pars_f !< input variable for vegetation parameters
    type(pars) :: water_pars_f !< input variable for water parameters

    type(pars_surf) :: building_surface_pars_f !< input variable for building surface parameters

    type(chem_emis_att_type) :: chem_emis_att !< Input Information of Chemistry Emission Data from
    !< netcdf
    type(chem_emis_val_type), allocatable, dimension(:) :: chem_emis !< Input Chemistry Emission Data from netcdf

    character(LEN=3) :: char_lod = 'lod' !< name of level-of-detail attribute in NetCDF file

    character(LEN=10) :: char_fill = '_FillValue' !< name of fill value attribute in NetCDF file

    character(LEN=100) :: input_file_static = 'PIDS_STATIC' !< Name of file which comprises static input data
    character(LEN=100) :: input_file_dynamic = 'PIDS_DYNAMIC' !< Name of file which comprises dynamic input data
    character(LEN=100) :: input_file_chem = 'PIDS_CHEM' !< Name of file which comprises chemistry input data
    character(LEN=100) :: input_file_vm = 'PIDS_VM' !< Name of file which comprises virtual measurement data
    character(LEN=100) :: input_file_wtm = 'PIDS_WTM' !< Name of file which comprises wind turbine model input data

    character(LEN=25), dimension(:), allocatable :: string_values !< output of string variables read from netcdf input files
    character(LEN=50), dimension(:), allocatable :: vars_pids !< variable in input file

    integer(iwp) :: id_emis !< NetCDF id of input file for chemistry emissions: TBD: It  has to be removed

    integer(iwp) :: nc_stat !< return value of nf90 function call
    integer(iwp) :: num_var_pids !< number of variables in file
    integer(iwp) :: pids_id !< file id

    logical :: input_pids_static = .true. !< Flag indicating whether Palm-input-data-standard file containing static
    !< information exists
    logical :: input_pids_dynamic = .false. !< Flag indicating whether Palm-input-data-standard file containing dynamic
    !< information exists
    logical :: input_pids_chem = .false. !< Flag indicating whether Palm-input-data-standard file containing chemistry
    !< information exists
    logical :: input_pids_vm = .false. !< Flag indicating whether input file for virtual measurements exist
    logical :: input_pids_wtm = .false. !< Flag indicating whether input file for wind turbine model exists

    logical :: collective_read = .false. !< Enable NetCDF collective read

    real(wp), dimension(8) :: crs_list !< list of coord_ref_sys values

    type(global_atts_type) :: input_file_atts !< global attributes of input file

    save

    private

    interface netcdf_data_input_1d
        module procedure netcdf_data_input_1d_real
    end interface netcdf_data_input_1d

    interface netcdf_data_input_2d_xy_real
        module procedure netcdf_data_input_2d_xy_real
    end interface netcdf_data_input_2d_xy_real

    interface netcdf_data_input_check_dynamic
        module procedure netcdf_data_input_check_dynamic
    end interface netcdf_data_input_check_dynamic

    interface netcdf_data_input_check_static
        module procedure netcdf_data_input_check_static
    end interface netcdf_data_input_check_static

    ! interface netcdf_data_input_chemistry_data
    !     module procedure netcdf_data_input_chemistry_data
    ! end interface netcdf_data_input_chemistry_data

    interface netcdf_data_input_parameter_lists
        module procedure netcdf_data_input_parameter_lists
    end interface netcdf_data_input_parameter_lists

    interface get_dimension_length
        module procedure get_dimension_length
    end interface get_dimension_length

    interface inquire_fill_value
        module procedure inquire_fill_value_int
        module procedure inquire_fill_value_real
    end interface inquire_fill_value

    interface netcdf_data_input_inquire_file
        module procedure netcdf_data_input_inquire_file
    end interface netcdf_data_input_inquire_file

    interface netcdf_data_input_init
        module procedure netcdf_data_input_init
    end interface netcdf_data_input_init

    interface netcdf_data_input_init_3d
        module procedure netcdf_data_input_init_3d
    end interface netcdf_data_input_init_3d

    interface netcdf_data_input_surface_data
        module procedure netcdf_data_input_surface_data
    end interface netcdf_data_input_surface_data

    interface get_variable
        module procedure get_variable_1d_char
        module procedure get_variable_1d_int
        module procedure get_variable_1d_real
        module procedure get_variable_2d_int8
        module procedure get_variable_2d_int32
        module procedure get_variable_2d_real
        module procedure get_variable_2d_list_int
        module procedure get_variable_2d_list_real
        module procedure get_variable_2d_real_dynamic
        module procedure get_variable_2d_real_time_slice
        module procedure get_variable_3d_int8
        module procedure get_variable_3d_real
        module procedure get_variable_3d_real_dynamic
        module procedure get_variable_4d_real_dynamic
        module procedure get_variable_4d_to_3d_real
        module procedure get_variable_4d_real
        module procedure get_variable_5d_to_4d_real
        module procedure get_variable_5d_real ! temp subroutine 4 reading 5D NC arrays
        module procedure get_variable_5d_real_dynamic ! 2B removed as z is out of emission_values
        module procedure get_variable_string
        module procedure get_variable_string_generic ! generic string function

    end interface get_variable

    interface get_variable_pr
        module procedure get_variable_pr
    end interface get_variable_pr

    interface get_attribute
        module procedure get_attribute_real
        module procedure get_attribute_int8
        module procedure get_attribute_int32
        module procedure get_attribute_string
    end interface get_attribute

    interface list_building_ids
        module procedure list_building_ids
    end interface list_building_ids

    !-- Public data structures
    public dims_xy, &
        pars, &
        int_2d_8bit, &
        real_1d_3d, &
        real_2d, &
        real_3d
    !-- Public variables
    public albedo_pars_f, &
        albedo_type_f, &
        buildings_f, &
        building_id_f, &
        building_surface_pars_f, &
        building_type_f, &
        char_fill, &
        char_lod, &
        chem_emis, &
        chem_emis_att, &
        chem_emis_att_type, &
        chem_emis_val_type, &
        coord_ref_sys, &
        crs_list, &
        init_3d, &
        init_model, &
        input_file_atts, &
        input_file_dynamic, &
        input_file_static, &
        input_file_vm, &
        input_file_wtm, &
        input_pids_static, &
        input_pids_dynamic, &
        input_pids_vm, &
        input_pids_wtm, &
        num_var_pids, &
        pavement_pars_f, &
        pavement_subsurface_pars_f, &
        pavement_type_f, &
        pids_id, &
        root_area_density_lsm_f, &
        soil_pars_f, &
        soil_type_f, &
        street_crossing_f, &
        street_type_f, &
        surface_fraction_f, &
        terrain_height_f, &
        vars_pids, &
        vegetation_pars_f, &
        vegetation_type_f, &
        water_pars_f, &
        water_type_f
    !-- Public subroutines
    public check_existence, &
        close_input_file, &
        get_attribute, &
        get_dimension_length, &
        get_variable, &
        get_variable_pr, &
        inquire_fill_value, &
        inquire_num_variables, &
        inquire_variable_names, &
        netcdf_data_input_1d, &
        netcdf_data_input_2d_xy_real, &
        list_building_ids, &
        netcdf_data_input_check_dynamic, &
        netcdf_data_input_check_static, &
        ! netcdf_data_input_chemistry_data, &
        netcdf_data_input_inquire_file, &
        netcdf_data_input_init, &
        netcdf_data_input_init_3d, &
        netcdf_data_input_parameter_lists, &
        netcdf_data_input_surface_data, &
        open_read_file

contains

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Reads single 1D variables from an input file.
    !--------------------------------------------------------------------------------------------------!
    subroutine netcdf_data_input_1d_real(pids_filename, var_in, dimname, val_1d)

        character(LEN=*) :: dimname !< corresponding dimension name
        character(LEN=*) :: pids_filename !< name of input file
        character(LEN=*) :: var_in !< input variable

        integer(iwp) :: nd1 !< dimension size

        real(wp), dimension(:), allocatable :: val_1d !< data array

        ! #if defined( __netcdf )
        call open_read_file(trim(pids_filename), pids_id)

        call inquire_num_variables(pids_id, num_var_pids)
        !-- Allocate memory to store variable names and read them.
        allocate (vars_pids(1:num_var_pids))
        call inquire_variable_names(pids_id, vars_pids)

        if (check_existence(vars_pids, trim(var_in))) then
            !--    Inquire corresponding dimension length.
            call get_dimension_length(pids_id, nd1, trim(dimname))
            !--    Allocate space for the variable to read.
            allocate (val_1d(0:nd1 - 1))
            !--    Read data.
            call get_variable(pids_id, trim(var_in), val_1d)

        end if

        deallocate (vars_pids)

        call close_input_file(pids_id)

    ! #endif
    end subroutine netcdf_data_input_1d_real

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Reads single 2D x-y variables from an input file.
    !--------------------------------------------------------------------------------------------------!
    subroutine netcdf_data_input_2d_xy_real(pids_filename, var_in, tmp, var_available, &
                                            includes_ghost_layers)

        character(LEN=*) :: pids_filename !< name of input file
        character(LEN=*) :: var_in !< input variable

        logical, optional :: includes_ghost_layers !< flag to steer resizing of the array including ghost points
        logical :: var_available !< flag to indicate whether a variable is present in the static input file

        type(real_2d), intent(INOUT) :: tmp !< data array

        ! #if defined( __netcdf )
        call open_read_file(trim(pids_filename), pids_id)

        call inquire_num_variables(pids_id, num_var_pids)
        !-- Allocate memory to store variable names and read them.
        allocate (vars_pids(1:num_var_pids))
        call inquire_variable_names(pids_id, vars_pids)

        var_available = check_existence(vars_pids, trim(var_in))
        if (var_available) then
            !--    Check if 2D array has been already allocated. If this is the case, re-allocate array.
            if (allocated(tmp % var)) deallocate (tmp % var)
            if (present(includes_ghost_layers)) then
                if (includes_ghost_layers) then
                    allocate (tmp % var(nysg:nyng, nxlg:nxrg))
                else
                    allocate (tmp % var(nys:nyn, nxl:nxr))
                end if
            else
                allocate (tmp % var(nys:nyn, nxl:nxr))
            end if
            !--    Read _FillValue attribute.
            call get_attribute(pids_id, char_fill, tmp % fill, .false., trim(var_in))
            !--    Read data.
            if (present(includes_ghost_layers)) then
                if (includes_ghost_layers) then
                    call get_variable(pids_id, trim(var_in), tmp % var, nxl, nxr, nys, nyn, nbgp=nbgp)
                else
                    call get_variable(pids_id, trim(var_in), tmp % var, nxl, nxr, nys, nyn, nbgp=0)
                end if
            else
                call get_variable(pids_id, trim(var_in), tmp % var, nxl, nxr, nys, nyn, nbgp=0)
            end if

            if (present(includes_ghost_layers)) then
                if (includes_ghost_layers) then
                    call exchange_horiz_2d(tmp % var)
                    call set_lateral_neumann_bc(tmp % var)
                end if
            end if
        end if

        deallocate (vars_pids)

        call close_input_file(pids_id)

    ! #endif
    end subroutine netcdf_data_input_2d_xy_real

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Inquires whether NetCDF input files according to Palm-input-data standard exist. Moreover, basic
    !> checks are performed.
    !--------------------------------------------------------------------------------------------------!
    subroutine netcdf_data_input_inquire_file

        implicit none

        ! #if defined ( __netcdf )
        inquire (FILE=trim(input_file_static)//trim(coupling_char), &
                 EXIST=input_pids_static)
        inquire (FILE=trim(input_file_dynamic)//trim(coupling_char), &
                 EXIST=input_pids_dynamic)
        inquire (FILE=trim(input_file_chem)//trim(coupling_char), &
                 EXIST=input_pids_chem)
        inquire (FILE=trim(input_file_vm)//trim(coupling_char), &
                 EXIST=input_pids_vm)
        inquire (FILE=trim(input_file_wtm)//trim(coupling_char), &
                 EXIST=input_pids_wtm)

        ! #endif
        !-- As long as topography can be input via ASCII format, no distinction between building and terrain
        !-- can be done. This case, classify all surfaces as default type. Same in case land-surface and
        !-- urban-surface model are not applied.
        if (.not. input_pids_static) then
            topo_no_distinct = .true.
        end if

    end subroutine netcdf_data_input_inquire_file

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Reads global attributes and coordinate reference system required for initialization of the model.
    !--------------------------------------------------------------------------------------------------!
    subroutine netcdf_data_input_init

        implicit none

        character(LEN=500) :: faulty_attribute !< string to gather attributes with incorrect settings

        integer(iwp) :: id_mod !< NetCDF id of input file
        integer(iwp) :: var_id_crs !< NetCDF id of variable crs

        logical :: input_pids_static_for_all_domains !< parameter to check if all domains have a static input file
        logical :: trigger_error_message !< control flag to check whether all attributes are set appropriately

        !-- Define default list of crs attributes. This is required for coordinate transformation.
        crs_list = (/coord_ref_sys % semi_major_axis, &
                     coord_ref_sys % inverse_flattening, &
                     coord_ref_sys % longitude_of_prime_meridian, &
                     coord_ref_sys % longitude_of_central_meridian, &
                     coord_ref_sys % scale_factor_at_central_meridian, &
                     coord_ref_sys % latitude_of_projection_origin, &
                     coord_ref_sys % false_easting, &
                     coord_ref_sys % false_northing/)

        if (input_pids_static) then

            ! #if defined ( __netcdf )
            !--    Open file in read-only mode
            call open_read_file(trim(input_file_static)//trim(coupling_char), id_mod)
            !--    Read global attributes
            call get_attribute(id_mod, input_file_atts % origin_lat_char, &
                               input_file_atts % origin_lat, .true.)

            call get_attribute(id_mod, input_file_atts % origin_lon_char, &
                               input_file_atts % origin_lon, .true.)

            call get_attribute(id_mod, input_file_atts % origin_time_char, &
                               input_file_atts % origin_time, .true., ignore_error=.true.)

            call get_attribute(id_mod, input_file_atts % origin_x_char, &
                               input_file_atts % origin_x, .true.)

            call get_attribute(id_mod, input_file_atts % origin_y_char, &
                               input_file_atts % origin_y, .true.)

            call get_attribute(id_mod, input_file_atts % origin_z_char, &
                               input_file_atts % origin_z, .true.)

            call get_attribute(id_mod, input_file_atts % rotation_angle_char, &
                               input_file_atts % rotation_angle, .true.)

            call get_attribute(id_mod, input_file_atts % author_char, &
                               input_file_atts % author, .true., ignore_error=.true.)

            call get_attribute(id_mod, input_file_atts % contact_person_char, &
                               input_file_atts % contact_person, .true., ignore_error=.true.)
            call get_attribute(id_mod, input_file_atts % institution_char, &
                               input_file_atts % institution, .true., ignore_error=.true.)
            call get_attribute(id_mod, input_file_atts % acronym_char, &
                               input_file_atts % acronym, .true., ignore_error=.true.)

            call get_attribute(id_mod, input_file_atts % campaign_char, &
                               input_file_atts % campaign, .true., ignore_error=.true.)
            call get_attribute(id_mod, input_file_atts % location_char, &
                               input_file_atts % location, .true., ignore_error=.true.)
            call get_attribute(id_mod, input_file_atts % site_char, &
                               input_file_atts % site, .true., ignore_error=.true.)

            call get_attribute(id_mod, input_file_atts % source_char, &
                               input_file_atts % source, .true., ignore_error=.true.)
            call get_attribute(id_mod, input_file_atts % references_char, &
                               input_file_atts % references, .true., ignore_error=.true.)
            call get_attribute(id_mod, input_file_atts % keywords_char, &
                               input_file_atts % keywords, .true., ignore_error=.true.)
            call get_attribute(id_mod, input_file_atts % licence_char, &
                               input_file_atts % licence, .true., ignore_error=.true.)
            call get_attribute(id_mod, input_file_atts % comment_char, &
                               input_file_atts % comment, .true., ignore_error=.true.)
            !--    Check for proper setting of global attributes, i.e. they must not be set to NaN values.
            trigger_error_message = .false.
            faulty_attribute = ''
            if (IEEE_IS_NAN(input_file_atts % origin_lat)) then
                trigger_error_message = .true.
                faulty_attribute = trim(faulty_attribute) // input_file_atts % origin_lat_char // ', '
            end if
            if (IEEE_IS_NAN(input_file_atts % origin_lon)) then
                trigger_error_message = .true.
                faulty_attribute = trim(faulty_attribute) // input_file_atts % origin_lon_char // ', '
            end if
            if (IEEE_IS_NAN(input_file_atts % origin_x)) then
                trigger_error_message = .true.
                faulty_attribute = trim(faulty_attribute) // input_file_atts % origin_x_char // ', '
            end if
            if (IEEE_IS_NAN(input_file_atts % origin_y)) then
                trigger_error_message = .true.
                faulty_attribute = trim(faulty_attribute) // input_file_atts % origin_y_char // ', '
            end if
            if (IEEE_IS_NAN(input_file_atts % origin_z)) then
                trigger_error_message = .true.
                faulty_attribute = trim(faulty_attribute) // input_file_atts % origin_z_char // ', '
            end if
            if (IEEE_IS_NAN(input_file_atts % rotation_angle)) then
                trigger_error_message = .true.
                faulty_attribute = trim(faulty_attribute) // input_file_atts % rotation_angle_char // ', '
            end if

            if (trigger_error_message) then
                print *, 'static driver global attribute(s): &' // trim(faulty_attribute) // &
                                 '&is (are) "NaN"'
                ! ! call message('netcdf_data_input_init', 'DRV0042', 1, 2, 0, 6, 0)
            end if
            !--    Read coordinate reference system if available
            nc_stat = NF90_INQ_VARID(id_mod, 'crs', var_id_crs)
            if (nc_stat == NF90_NOERR) then
                call get_attribute(id_mod, 'epsg_code', coord_ref_sys % epsg_code, .false., 'crs')
                call get_attribute(id_mod, 'false_easting', coord_ref_sys % false_easting, .false., 'crs')
                call get_attribute(id_mod, 'false_northing', coord_ref_sys % false_northing, .false., 'crs')
                call get_attribute(id_mod, 'grid_mapping_name', coord_ref_sys % grid_mapping_name, &
                                   .false., 'crs')
                call get_attribute(id_mod, 'inverse_flattening', coord_ref_sys % inverse_flattening, &
                                   .false., 'crs')
                call get_attribute(id_mod, 'latitude_of_projection_origin', &
                                   coord_ref_sys % latitude_of_projection_origin, .false., 'crs')
                call get_attribute(id_mod, 'long_name', coord_ref_sys % long_name, .false., 'crs')
                call get_attribute(id_mod, 'longitude_of_central_meridian', &
                                   coord_ref_sys % longitude_of_central_meridian, .false., 'crs')
                call get_attribute(id_mod, 'longitude_of_prime_meridian', &
                                   coord_ref_sys % longitude_of_prime_meridian, .false., 'crs')
                call get_attribute(id_mod, 'scale_factor_at_central_meridian', &
                                   coord_ref_sys % scale_factor_at_central_meridian, .false., 'crs')
                call get_attribute(id_mod, 'semi_major_axis', coord_ref_sys % semi_major_axis, .false., &
                                   'crs')
                call get_attribute(id_mod, 'units', coord_ref_sys % units, .false., 'crs')
            else
                !--       Calculate central meridian from origin_lon
                coord_ref_sys % longitude_of_central_meridian = &
                    ceiling(input_file_atts % origin_lon / 6.0_wp) * 6.0_wp - 3.0_wp
            end if
            !--    Finally, close input file
            call close_input_file(id_mod)

            ! #endif
            !--    Copy latitude, longitude, origin_z, rotation angle on init type.
            !--    Note: A shifting height might have already been saved to orgin_z in init_grid; therefore,
            !--          do not override but add the reference height from the input file.
            init_model % latitude = input_file_atts % origin_lat
            init_model % longitude = input_file_atts % origin_lon
            init_model % origin_x = input_file_atts % origin_x
            init_model % origin_y = input_file_atts % origin_y
            init_model % origin_z = init_model % origin_z + input_file_atts % origin_z
            init_model % rotation_angle = input_file_atts % rotation_angle

            !--    Update list of crs attributes. This is required for coordinate transformation.
            crs_list = (/coord_ref_sys % semi_major_axis, &
                         coord_ref_sys % inverse_flattening, &
                         coord_ref_sys % longitude_of_prime_meridian, &
                         coord_ref_sys % longitude_of_central_meridian, &
                         coord_ref_sys % scale_factor_at_central_meridian, &
                         coord_ref_sys % latitude_of_projection_origin, &
                         coord_ref_sys % false_easting, &
                         coord_ref_sys % false_northing/)
        end if
        !-- In case of nested runs, each model domain might have different longitude and latitude, which
        !-- would result in different Coriolis parameters and sun-zenith angles. To avoid this, longitude
        !-- and latitude in each model domain will be set to the values of the root model. Please note,
        !-- this synchronization is required already here. However, synchronization is only done if
        !-- these information is provided by a static input file in all domains. Therefore, first
        !-- check if a static input file is available for all domains.
        input_pids_static_for_all_domains = input_pids_static


        ! #if defined( __parallel )
        !         call MPI_ALLREDUCE(MPI_IN_PLACE, input_pids_static_for_all_domains, 1, MPI_LOGICAL, MPI_LAND, &
        !                            MPI_COMM_WORLD, ierr)
        ! #endif
        if (input_pids_static_for_all_domains) then
        ! #if defined( __parallel )
        !             call MPI_BCAST(init_model % latitude, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        !             call MPI_BCAST(init_model % longitude, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        ! #endif
        end if

    end subroutine netcdf_data_input_init

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Reads Chemistry NETCDF Input data, such as emission values, emission species, etc.
    !--------------------------------------------------------------------------------------------------!
    ! subroutine netcdf_data_input_chemistry_data(emt_att, emt)

    !     use chem_modules, &
    !         only: emiss_lod, &
    !               surface_csflux_name, &
    !               time_fac_type

    !     implicit none

    !     type(chem_emis_att_type), intent(INOUT) :: emt_att
    !     type(chem_emis_val_type), allocatable, dimension(:), intent(INOUT) :: emt

    !     integer(iwp) :: i, j, k !< generic counters
    !     integer(iwp) :: ispec !< index for number of emission species in input
    !     integer(iwp) :: len_dims !< Length of dimension
    !     integer(iwp) :: num_vars !< number of variables in netcdf input file


    !     !-- dum_var_4d are designed to read in emission_values from the chemistry netCDF file.
    !     !-- Currently the vestigial "z" dimension in emission_values makes it a 5D array, hence the
    !     !-- corresponding dum_var_5d array. When the "z" dimension is removed completely, dum_var_4d will be
    !     !-- used instead
    !     !    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:)    ::  dum_var_4d  !< temp array 4 4D chem emission data
    !     real(wp), allocatable, dimension(:, :, :, :, :) :: dum_var_5d !< temp array 4 5D chem emission data

    !     !-- Start processing data
    !     !-- Emission LOD 0 (Parameterized mode)
    !     if (emiss_lod == 0) then


    !         !--    For reference
    !         !       IF (TRIM(mode_emis) == "PARAMETERIZED" .OR. TRIM(mode_emis) == "parameterized") THEN
    !         ispec = 1
    !         emt_att % n_emiss_species = 0

    !         !--    Number of species
    !         do while (trim(surface_csflux_name(ispec)) /= 'novalue')

    !             emt_att % n_emiss_species = emt_att % n_emiss_species + 1
    !             ispec = ispec + 1
    !             !--       Followling line retained for compatibility with salsa_mod which still uses emt_att%nspec
    !             !--       heavily
    !             emt_att % nspec = emt_att % nspec + 1

    !         end do

    !         !--    Allocate emission values data type arrays
    !         allocate (emt(emt_att % n_emiss_species))


    !         !--    Read emission species names
    !         !--    Allocate space for strings
    !         allocate (emt_att % species_name(emt_att % n_emiss_species))

    !         do ispec = 1, emt_att % n_emiss_species
    !             emt_att % species_name(ispec) = trim(surface_csflux_name(ispec))
    !         end do

    !     !-- LOD 1 (default mode) and LOD 2 (pre-processed mode)
    !     else


    !         ! #if defined ( __netcdf )
    !         if (.not. input_pids_chem) return


    !         !--    First we allocate memory space for the emission species and then we differentiate between
    !         !--    LOD 1 (default mode) and LOD 2 (pre-processed mode)
    !         !--    Open emission data file ( {palmcase}_chemistry )
    !         call open_read_file(trim(input_file_chem)//trim(coupling_char), id_emis)

    !         !--    Inquire number of variables
    !         call inquire_num_variables(id_emis, num_vars)

    !         !--    Get general dimension lengths: only # species and # categories.
    !         !--    Other dimensions depend on the emission mode or specific components.
    !         call get_dimension_length(id_emis, emt_att % n_emiss_species, 'nspecies')

    !         !--    Backward compatibility for salsa_mod
    !         emt_att % nspec = emt_att % n_emiss_species

    !         !--    Allocate emission values data type arrays
    !         allocate (emt(emt_att % n_emiss_species))


    !         !--    Reading in species names
    !         !--    Allocate memory for species names
    !         allocate (emt_att % species_name(emt_att % n_emiss_species))

    !         !--    Retrieve variable name (again, should use n_emiss_strlen)
    !         call get_variable(id_emis, 'emission_name', string_values, emt_att % n_emiss_species)
    !         emt_att % species_name = string_values

    !         !--    Dealocate string_values previously allocated in get_variable call
    !         if (allocated(string_values)) deallocate (string_values)


    !         !--    Reading in species indices
    !         !--    Allocate memory for species indices
    !         allocate (emt_att % species_index(emt_att % n_emiss_species))

    !         !--    Retrieve variable data
    !         call get_variable(id_emis, 'emission_index', emt_att % species_index)

    !         !--    Now the routine has to distinguish between chemistry emission LOD 1 (default mode) and LOD 2
    !         !--    (pre-processed mode)
    !         !--    Start of emission LOD 1 (default mode)
    !         if (emiss_lod == 1) then

    !             !--       For reference
    !             !          IF (TRIM(mode_emis) == "DEFAULT" .OR. TRIM(mode_emis) == "default") THEN
    !             !--       Get number of emission categories
    !             call get_dimension_length(id_emis, emt_att % ncat, 'ncat')

    !             !--       Reading in emission categories indices
    !             allocate (emt_att % cat_index(emt_att % ncat))

    !             !--       Retrieve variable data
    !             call get_variable(id_emis, 'emission_cat_index', emt_att % cat_index)



    !             !-- Loop through individual species to get basic information on VOC/PM/NOX/SOX
    !             !---------------------------------------------------------------------------------------------------
    !             !-- Note - Check array indices for reading in names and species in LOD1 (default mode) for the
    !             !--        various mode splits as all id_emis conditionals have been removed from get_var_functions.
    !             !--        In theory this would mean all arrays should be read from 0 to n-1 (C convention) as
    !             !--        opposed to 1 to n (FORTRAN convention). Keep this in mind !!
    !             !--------------------------------------------------------------------------------------------------
    !             do ispec = 1, emt_att % n_emiss_species

    !                 !--          VOC data (name and composition)
    !                 if (trim(emt_att % species_name(ispec)) == "VOC" .or. &
    !                     trim(emt_att % species_name(ispec)) == "voc") then

    !                     !--             VOC name
    !                     call get_dimension_length(id_emis, emt_att % nvoc, 'nvoc')
    !                     allocate (emt_att % voc_name(emt_att % nvoc))
    !                     call get_variable(id_emis, "emission_voc_name", string_values, emt_att % nvoc)
    !                     emt_att % voc_name = string_values
    !                     if (allocated(string_values)) deallocate (string_values)

    !                     !--             VOC composition
    !                     allocate (emt_att % voc_comp(emt_att % ncat, emt_att % nvoc))
    !                     call get_variable(id_emis, "composition_voc", emt_att % voc_comp, 1, emt_att % ncat, &
    !                                       1, emt_att % nvoc, nbgp=0)

    !                 end if ! VOC

    !                 !--          PM data (name and composition)
    !                 if (trim(emt_att % species_name(ispec)) == "PM" .or. &
    !                     trim(emt_att % species_name(ispec)) == "pm") then

    !                     !--             PM name
    !                     call get_dimension_length(id_emis, emt_att % npm, 'npm')
    !                     allocate (emt_att % pm_name(emt_att % npm))
    !                     call get_variable(id_emis, "pm_name", string_values, emt_att % npm)
    !                     emt_att % pm_name = string_values
    !                     if (allocated(string_values)) deallocate (string_values)

    !                     !--             PM composition (PM1, PM2.5 and PM10)
    !                     len_dims = 3 ! PM1, PM2.5, PM10
    !                     allocate (emt_att % pm_comp(emt_att % ncat, emt_att % npm, len_dims))
    !                     call get_variable(id_emis, "composition_pm", emt_att % pm_comp, 1, emt_att % ncat, 1, &
    !                                       emt_att % npm, 1, len_dims, nbgp=0)

    !                 end if ! PM

    !                 !--          NOX (NO and NO2)
    !                 if (trim(emt_att % species_name(ispec)) == "NOX" .or. &
    !                     trim(emt_att % species_name(ispec)) == "nox") then

    !                     allocate (emt_att % nox_comp(emt_att % ncat, emt_att % nnox))
    !                     call get_variable(id_emis, "composition_nox", emt_att % nox_comp, 1, emt_att % ncat, &
    !                                       1, emt_att % nnox, nbgp=0)

    !                 end if ! NOX

    !                 !--          SOX (SO2 and SO4)
    !                 if (trim(emt_att % species_name(ispec)) == "SOX" .or. &
    !                     trim(emt_att % species_name(ispec)) == "sox") then

    !                     allocate (emt_att % sox_comp(emt_att % ncat, emt_att % nsox))
    !                     call get_variable(id_emis, "composition_sox", emt_att % sox_comp, 1, emt_att % ncat, &
    !                                       1, emt_att % nsox, nbgp=0)

    !                 end if ! SOX

    !             end do ! do ispec


    !             !--       Emission time scaling factors (hourly and MDH data)
    !             !--       Hour
    !             if (trim(time_fac_type) == "HOUR" .or. trim(time_fac_type) == "hour") then

    !                 call get_dimension_length(id_emis, emt_att % nhoursyear, 'nhoursyear')
    !                 allocate (emt_att % hourly_emis_time_factor(emt_att % ncat, emt_att % nhoursyear))
    !                 call get_variable(id_emis, "emission_time_factors", emt_att % hourly_emis_time_factor, &
    !                                   1, emt_att % ncat, 1, emt_att % nhoursyear, nbgp=0)

    !             !--       MDH
    !             elseif (trim(time_fac_type) == "MDH" .or. trim(time_fac_type) == "mdh") then

    !                 call get_dimension_length(id_emis, emt_att % nmonthdayhour, 'nmonthdayhour')
    !                 allocate (emt_att % mdh_emis_time_factor(emt_att % ncat, emt_att % nmonthdayhour))
    !                 call get_variable(id_emis, "emission_time_factors", emt_att % mdh_emis_time_factor, &
    !                                   1, emt_att % ncat, 1, emt_att % nmonthdayhour, nbgp=0)

    !             !--       Error (time factor undefined)
    !             else

    !                 print *, 'We are in the DEFAULT chemistry emissions mode: ' // &
    !                                  '     !no time-factor type specified!' // &
    !                                  'Please specify the value of time_fac_type:' // &
    !                                  '         either "MDH" or "HOUR"'
    !                 ! call message('netcdf_data_input_chemistry_data', 'DRV0002', 2, 2, 0, 6, 0)

    !             end if ! time_fac_type


    !             !--       Read in default (LOD1) emissions from chemisty netCDF file per species
    !             !--       Note - At the moment the data is read in per species, but in the future it would be much
    !             !--              more sensible to read in per species per time step to reduce memory consumption
    !             !--              and, to a smaller degree, dimensionality of data exchange (I expect this will be
    !             !--              necessary when the problem size is large).
    !             do ispec = 1, emt_att % n_emiss_species

    !                 !--          Allocate space for species specific emission values.
    !                 !--          Note - This array is extended by 1 cell in each horizontal direction to compensate for
    !                 !--                 an apparent linear offset. The reason of this offset is not known but it has
    !                 !--                 been determined to take place beyond the scope of this module, and has little to
    !                 !--                 do with index conventions.
    !                 !--                 That is, setting the array horizontal limit from nx0:nx1 to 1:(nx1-nx0+1) or
    !                 !--                 nx0+1:nx1+1 did not result in correct or definite behavior.
    !                 !--                 This must be looked at at some point by the Hannover team but for now this
    !                 !--                 workaround is deemed reasonable.
    !                 if (.not. allocated(emt(ispec) % default_emission_data)) then
    !                     allocate (emt(ispec) % default_emission_data(emt_att % ncat, nys:nyn + 1, nxl:nxr + 1))
    !                 end if
    !                 !--          Allocate dummy variable w/ index order identical to that shown in the netCDF header
    !                 allocate (dum_var_5d(1, nys:nyn, nxl:nxr, 1, emt_att % ncat))
    !                 !--          Get variable.  Be very careful
    !                 !--          I am using get_variable_5d_real_dynamic (note logical argument at the end)
    !                 !--          1) use Fortran index convention (i.e., 1 to N)
    !                 !--          2) index order must be in reverse order from above allocation order
    !                 call get_variable(id_emis, "emission_values", dum_var_5d, 1, ispec, &
    !                                   nxl + 1, nys + 1, 1, emt_att % ncat, 1, &
    !                                   nxr - nxl + 1, nyn - nys + 1, emt_att % dt_emission, .false.)
    !                 !--          Assign temp array to data structure then deallocate temp array
    !                 !--          Note - Indices are shifted from nx0:nx1 to nx0+1:nx1+1 to offset the emission data
    !                 !--                 array to counter said domain offset.
    !                 do k = 1, emt_att % ncat
    !                     do j = nys + 1, nyn + 1
    !                         do i = nxl + 1, nxr + 1
    !                             emt(ispec) % default_emission_data(k, j, i) = dum_var_5d(1, j - 1, i - 1, 1, k)
    !                         end do
    !                     end do
    !                 end do

    !                 deallocate (dum_var_5d)

    !             end do ! ispec
    !             !--       Units
    !             call get_attribute(id_emis, "units", emt_att % units, .false., "emission_values")


    !         !--       End default mode
    !         !--    Start LOD 2 (pre-processed mode)
    !         elseif (emiss_lod == 2) then


    !             !      For reference
    !             !       ELSEIF( TRIM( mode_emis ) == "PRE-PROCESSED" .OR. TRIM( mode_emis ) == "pre-processed")  THEN
    !             !--       For LOD 2 only VOC and emission data need to be read
    !             !--       Note - Check array indices for reading in names and species in LOD2 (pre-processed mode)
    !             !--       for the various mode splits as all id_emis conditionals have been removed from
    !             !--       get_var_functions. In theory this would mean all arrays should be read from 0 to n-1
    !             !--       (C convention) as opposed to 1 to n (FORTRAN convention). Keep this in mind!!
    !             do ispec = 1, emt_att % n_emiss_species

    !                 !--          VOC data (name and composition)
    !                 if (trim(emt_att % species_name(ispec)) == "VOC" .or. &
    !                     trim(emt_att % species_name(ispec)) == "voc") then

    !                     !--             VOC name
    !                     call get_dimension_length(id_emis, emt_att % nvoc, 'nvoc')
    !                     allocate (emt_att % voc_name(emt_att % nvoc))
    !                     call get_variable(id_emis, "emission_voc_name", string_values, emt_att % nvoc)
    !                     emt_att % voc_name = string_values
    !                     if (allocated(string_values)) deallocate (string_values)

    !                     !--             VOC composition
    !                     allocate (emt_att % voc_comp(emt_att % ncat, emt_att % nvoc))
    !                     call get_variable(id_emis, "composition_voc", emt_att % voc_comp, 1, emt_att % ncat, &
    !                                       1, emt_att % nvoc, nbgp=0)
    !                 end if ! VOC

    !             end do ! ispec

    !             !--       Emission data
    !             call get_dimension_length(id_emis, emt_att % dt_emission, 'time')


    !             !--       Read in pre-processed (LOD2) emissions from chemisty netCDF file per species
    !             !--       Note - At the moment the data is read in per species, but in the future it would be much
    !             !--              more sensible to read in per species per time step to reduce memory consumption
    !             !--              and, to a smaller degree, dimensionality of data exchange (I expect this will be
    !             !--              necessary when the problem size is large).
    !             do ispec = 1, emt_att % n_emiss_species


    !                 !--          Allocate space for species specific emission values.
    !                 !--          Note - This array is extended by 1 cell in each horizontal direction to compensate for
    !                 !--                 an apparent linear offset. The reason of this offset is not known but it has
    !                 !--                 been determined to take place beyond the scope of this module, and has little to
    !                 !--                 do with index conventions.
    !                 !--                 That is, setting the array horizontal limit from nx0:nx1 to 1:(nx1-nx0+1) or
    !                 !--                 nx0+1:nx1+1 did not result in correct or definite behavior.
    !                 !--                 This must be looked at at some point by the Hannover team but for now this
    !                 !--                 workaround is deemed reasonable.
    !                 if (.not. allocated(emt(ispec) % preproc_emission_data)) then
    !                     allocate (emt(ispec) % preproc_emission_data( &
    !                               emt_att % dt_emission, 1, nys:nyn + 1, nxl:nxr + 1))
    !                 end if
    !                 !--          Allocate dummy variable w/ index order identical to that shown in the netCDF header
    !                 allocate (dum_var_5d(emt_att % dt_emission, 1, nys:nyn, nxl:nxr, 1))
    !                 !--          Get variable.  Be very careful.
    !                 !--          I am using get_variable_5d_real_dynamic (note logical argument at the end)
    !                 !--          1) use Fortran index convention (i.e., 1 to N)
    !                 !--          2) index order must be in reverse order from above allocation order
    !                 call get_variable(id_emis, "emission_values", dum_var_5d, ispec, &
    !                                   nxl + 1, nys + 1, 1, 1, 1, &
    !                                   nxr - nxl + 1, nyn - nys + 1, 1, emt_att % dt_emission, .false.)

    !                 !--          Assign temp array to data structure then deallocate temp array.
    !                 !--          Note - Indices are shifted from nx0:nx1 to nx0+1:nx1+1 to offset the emission data
    !                 !--                 array to counter mentioned unkonwn offset.
    !                 do k = 1, emt_att % dt_emission
    !                     do j = nys + 1, nyn + 1
    !                         do i = nxl + 1, nxr + 1
    !                             emt(ispec) % preproc_emission_data(k, 1, j, i) = dum_var_5d(k, 1, j - 1, i - 1, 1)
    !                         end do
    !                     end do
    !                 end do

    !                 deallocate (dum_var_5d)

    !             end do ! ispec
    !             !--       Units
    !             call get_attribute(id_emis, "units", emt_att % units, .false., "emission_values")

    !         end if ! LOD1 & LOD2 (default and pre-processed mode)

    !         call close_input_file(id_emis)


    !     ! #endif
    !     end if ! LOD0 (parameterized mode)

    ! end subroutine netcdf_data_input_chemistry_data

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Reads the parameter lists for each grid point (3D data).
    !--------------------------------------------------------------------------------------------------!
    subroutine netcdf_data_input_parameter_lists(pids_filename, var_in, dimname, dimname_layer, &
                                                 parlist)

        character(LEN=*), intent(IN) :: dimname !< main dimension name
        character(LEN=*), intent(IN), optional :: dimname_layer !< layer dimension name
        character(LEN=*), intent(IN) :: pids_filename !< name of input file
        character(LEN=*), intent(IN) :: var_in !< input variable

        type(pars), intent(INOUT) :: parlist !< input data type

        integer(iwp) :: k, kz !< running index along parameters and layers

        real(wp), dimension(:, :), allocatable :: tmp_2d !< temporary array to hold horizontal slices of 3d-data

        ! #if defined( __netcdf )
        call open_read_file(trim(pids_filename), pids_id)

        call inquire_num_variables(pids_id, num_var_pids)
        !-- Allocate memory to store variable names and read them.
        allocate (vars_pids(1:num_var_pids))
        call inquire_variable_names(pids_id, vars_pids)

        if (check_existence(vars_pids, trim(var_in))) then
            parlist % from_file = .true.
            !--    Inquire FillValue attribute.
            call get_attribute(pids_id, char_fill, parlist % fill, .false., trim(var_in))
            !--    Inquire corresponding dimension length.
            call get_dimension_length(pids_id, parlist % np, trim(dimname))

            if (present(dimname_layer)) then
                !--       Inquire layer dimension length.
                call get_dimension_length(pids_id, parlist % nz, trim(dimname_layer))

                !--       Read parameter list provided at x,y positions.
                allocate (parlist % pars_xyz(0:parlist % np - 1, 0:parlist % nz - 1, nysg:nyng, nxlg:nxrg))
                call get_variable(pids_id, trim(var_in), parlist % pars_xyz, nxl, nxr, nys, nyn, &
                                  0, parlist % nz - 1, 0, parlist % np - 1, nbgp=nbgp)

            else
                !--       Read parameter list provided at x,y positions.
                allocate (parlist % pars_xy(0:parlist % np - 1, nysg:nyng, nxlg:nxrg))
                call get_variable(pids_id, trim(var_in), parlist % pars_xy, nxl, nxr, nys, nyn, &
                                  0, parlist % np - 1, nbgp=nbgp)
            end if

        else
            parlist % from_file = .false.
        end if

        deallocate (vars_pids)

        call close_input_file(pids_id)

        ! #endif
        !-- Parameter lists that are read from file need to be extended by the ghost layers. This is
        !-- required because wall surfaces access the +-1 surrounding grid points.
        if (parlist % from_file) then
            allocate (tmp_2d(nysg:nyng, nxlg:nxrg))

            if (present(dimname_layer)) then
                do k = 0, parlist % np - 1
                    do kz = 0, parlist % nz - 1
                        tmp_2d(:, :) = parlist % pars_xyz(k, kz, :, :)
                        call exchange_horiz_2d(tmp_2d)
                        call set_lateral_neumann_bc(tmp_2d)
                        parlist % pars_xyz(k, kz, :, :) = tmp_2d(:, :)
                    end do
                end do
            else
                do k = 0, parlist % np - 1
                    tmp_2d(:, :) = parlist % pars_xy(k, :, :)
                    call exchange_horiz_2d(tmp_2d)
                    call set_lateral_neumann_bc(tmp_2d)
                    parlist % pars_xy(k, :, :) = tmp_2d(:, :)
                end do
            end if

            deallocate (tmp_2d)
        end if

    end subroutine netcdf_data_input_parameter_lists

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Reads surface classification data, such as vegetation and soil type, etc. .
    !--------------------------------------------------------------------------------------------------!
    subroutine netcdf_data_input_surface_data

        character(LEN=100), dimension(:), allocatable :: var_names !< variable names in static input file

        integer(iwp) :: id_surf !< NetCDF id of input file
        integer(iwp) :: k !< running index along z-direction
        integer(iwp) :: k2 !< running index
        integer(iwp) :: num_vars !< number of variables in input file
        integer(iwp) :: nz_soil !< number of soil layers in file

        integer(ibp), dimension(:, :), allocatable :: tmp_2d_byte !< temporary array to hold horizontal slices of 3d-data

        real(wp), dimension(:, :), allocatable :: tmp_2d !< temporary array to hold horizontal slices of 3d-data

        !-- If not static input file is available, skip this routine
        if (.not. input_pids_static) return
        !-- Measure CPU time
        ! call cpu_log(log_point_s(82), 'NetCDF input', 'start')
        !-- Skip the following if no land-surface or urban-surface module are applied. This case, no one of
        !-- the following variables is used anyway.
        if (.not. land_surface .and. .not. urban_surface) return

        ! #if defined ( __netcdf )
        !-- Open file in read-only mode
        call open_read_file(trim(input_file_static)//trim(coupling_char), id_surf)
        !-- Inquire all variable names.
        !-- This will be used to check whether an optional input variable exists or not.
        call inquire_num_variables(id_surf, num_vars)

        allocate (var_names(1:num_vars))
        call inquire_variable_names(id_surf, var_names)
        !-- Read vegetation type and required attributes
        if (check_existence(var_names, 'vegetation_type')) then
            vegetation_type_f % from_file = .true.
            call get_attribute(id_surf, char_fill, vegetation_type_f % fill, .false., 'vegetation_type')

            allocate (vegetation_type_f % var(nysg:nyng, nxlg:nxrg))

            call get_variable(id_surf, 'vegetation_type', vegetation_type_f % var, nxl, nxr, nys, nyn, &
                              nbgp=nbgp)
        else
            vegetation_type_f % from_file = .false.
        end if

        !-- Read soil type and required attributes
        if (check_existence(var_names, 'soil_type')) then
            soil_type_f % from_file = .true.
            !--    Note, lod is currently not on file; skip for the moment
            !       CALL get_attribute( id_surf, char_lod,                                                      &
            !                           soil_type_f%lod,                                                        &
            !                           .FALSE., 'soil_type' )
            call get_attribute(id_surf, char_fill, soil_type_f % fill, .false., 'soil_type')

            if (soil_type_f % lod == 1) then

                allocate (soil_type_f % var_2d(nysg:nyng, nxlg:nxrg))

                call get_variable(id_surf, 'soil_type', soil_type_f % var_2d, nxl, nxr, nys, nyn, &
                                  nbgp=nbgp)

            elseif (soil_type_f % lod == 2) then
                !--       Obtain number of soil layers from file.
                call get_dimension_length(id_surf, nz_soil, 'zsoil')

                allocate (soil_type_f % var_3d(0:nz_soil, nysg:nyng, nxlg:nxrg))

                call get_variable(id_surf, 'soil_type', soil_type_f % var_3d, nxl, nxr, nys, nyn, 0, &
                                  nz_soil, nbgp=nbgp)

            end if
        else
            soil_type_f % from_file = .false.
        end if

        !-- Read pavement type and required attributes
        if (check_existence(var_names, 'pavement_type')) then
            pavement_type_f % from_file = .true.
            call get_attribute(id_surf, char_fill, pavement_type_f % fill, .false., 'pavement_type')

            allocate (pavement_type_f % var(nysg:nyng, nxlg:nxrg))

            call get_variable(id_surf, 'pavement_type', pavement_type_f % var, nxl, nxr, nys, nyn, &
                              nbgp=nbgp)
        else
            pavement_type_f % from_file = .false.
        end if

        !-- Read water type and required attributes
        if (check_existence(var_names, 'water_type')) then
            water_type_f % from_file = .true.
            call get_attribute(id_surf, char_fill, water_type_f % fill, .false., 'water_type')

            allocate (water_type_f % var(nysg:nyng, nxlg:nxrg))

            call get_variable(id_surf, 'water_type', water_type_f % var, nxl, nxr, nys, nyn, nbgp=nbgp)

        else
            water_type_f % from_file = .false.
        end if
        !-- Read relative surface fractions of vegetation, pavement and water.
        if (check_existence(var_names, 'surface_fraction')) then
            surface_fraction_f % from_file = .true.
            call get_attribute(id_surf, char_fill, surface_fraction_f % fill, .false., 'surface_fraction')
            !--    Inquire number of surface fractions
            call get_dimension_length(id_surf, surface_fraction_f % nf, 'nsurface_fraction')
            !--    Allocate dimension array and input array for surface fractions
            allocate (surface_fraction_f % nfracs(0:surface_fraction_f % nf - 1))
            allocate (surface_fraction_f % frac(0:surface_fraction_f % nf - 1, nysg:nyng, nxlg:nxrg))
            !--    Get dimension of surface fractions
            call get_variable(id_surf, 'nsurface_fraction', surface_fraction_f % nfracs)
            !--    Read surface fractions
            call get_variable(id_surf, 'surface_fraction', surface_fraction_f % frac, &
                              nxl, nxr, nys, nyn, 0, surface_fraction_f % nf - 1, nbgp=nbgp)
        else
            surface_fraction_f % from_file = .false.
        end if
        !-- Read building surface parameters
        if (check_existence(var_names, 'building_surface_pars')) then
            building_surface_pars_f % from_file = .true.
            call get_attribute(id_surf, char_fill, building_surface_pars_f % fill, .false., &
                               'building_surface_pars')
            !--    Read building_surface_pars
            call get_variable_surf(id_surf, 'building_surface_pars', building_surface_pars_f)
        else
            building_surface_pars_f % from_file = .false.
        end if

        !-- Read albedo type and required attributes
        if (check_existence(var_names, 'albedo_type')) then
            albedo_type_f % from_file = .true.
            call get_attribute(id_surf, char_fill, albedo_type_f % fill, .false., 'albedo_type')

            allocate (albedo_type_f % var(nysg:nyng, nxlg:nxrg))

            call get_variable(id_surf, 'albedo_type', albedo_type_f % var, nxl, nxr, nys, nyn, nbgp=nbgp)
        else
            albedo_type_f % from_file = .false.
        end if

        !-- Read albedo parameters and related information
        if (check_existence(var_names, 'albedo_pars')) then
            albedo_pars_f % from_file = .true.
            call get_attribute(id_surf, char_fill, albedo_pars_f % fill, .false., 'albedo_pars')
            !--    Inquire number of albedo parameters
            call get_dimension_length(id_surf, albedo_pars_f % np, 'nalbedo_pars')
            !--    Allocate input array for albedo parameters and read them
            allocate (albedo_pars_f % pars_xy(0:albedo_pars_f % np - 1, nysg:nyng, nxlg:nxrg))
            call get_variable(id_surf, 'albedo_pars', albedo_pars_f % pars_xy, &
                              nxl, nxr, nys, nyn, 0, albedo_pars_f % np - 1, nbgp=nbgp)
        else
            albedo_pars_f % from_file = .false.
        end if

        !-- Read pavement parameters and related information
        if (check_existence(var_names, 'pavement_pars')) then
            pavement_pars_f % from_file = .true.
            call get_attribute(id_surf, char_fill, pavement_pars_f % fill, .false., 'pavement_pars')
            !--    Inquire number of pavement parameters
            call get_dimension_length(id_surf, pavement_pars_f % np, 'npavement_pars')
            !--    Allocate input array for pavement parameters and read them
            allocate (pavement_pars_f % pars_xy(0:pavement_pars_f % np - 1, nysg:nyng, nxlg:nxrg))
            call get_variable(id_surf, 'pavement_pars', pavement_pars_f % pars_xy, &
                              nxl, nxr, nys, nyn, 0, pavement_pars_f % np - 1, nbgp=nbgp)
        else
            pavement_pars_f % from_file = .false.
        end if

        !-- Read pavement subsurface parameters and related information
        if (check_existence(var_names, 'pavement_subsurface_pars')) then
            pavement_subsurface_pars_f % from_file = .true.
            call get_attribute(id_surf, char_fill, pavement_subsurface_pars_f % fill, .false., &
                               'pavement_subsurface_pars')
            !--    Inquire number of parameters
            call get_dimension_length(id_surf, pavement_subsurface_pars_f % np, &
                                      'npavement_subsurface_pars')
            !--    Inquire number of soil layers
            call get_dimension_length(id_surf, pavement_subsurface_pars_f % nz, 'zsoil')
            !--    Allocate input array for pavement parameters and read them
            allocate (pavement_subsurface_pars_f % pars_xyz(0:pavement_subsurface_pars_f % np - 1, &
                                                            0:pavement_subsurface_pars_f % nz - 1, &
                                                            nysg:nyng, nxlg:nxrg))
            call get_variable(id_surf, 'pavement_subsurface_pars', pavement_subsurface_pars_f % pars_xyz, &
                              nxl, nxr, nys, nyn, 0, pavement_subsurface_pars_f % nz - 1, &
                              0, pavement_subsurface_pars_f % np - 1, nbgp=nbgp)
        else
            pavement_subsurface_pars_f % from_file = .false.
        end if

        !-- Read vegetation parameters and related information
        if (check_existence(var_names, 'vegetation_pars')) then
            vegetation_pars_f % from_file = .true.
            call get_attribute(id_surf, char_fill, vegetation_pars_f % fill, .false., 'vegetation_pars')
            !--    Inquire number of vegetation parameters
            call get_dimension_length(id_surf, vegetation_pars_f % np, 'nvegetation_pars')
            !--    Allocate input array for surface fractions and read them
            allocate (vegetation_pars_f % pars_xy(0:vegetation_pars_f % np - 1, nysg:nyng, nxlg:nxrg))
            call get_variable(id_surf, 'vegetation_pars', vegetation_pars_f % pars_xy, &
                              nxl, nxr, nys, nyn, 0, vegetation_pars_f % np - 1, nbgp=nbgp)
        else
            vegetation_pars_f % from_file = .false.
        end if

        !-- Read root parameters/distribution and related information
        if (check_existence(var_names, 'soil_pars')) then
            soil_pars_f % from_file = .true.
            call get_attribute(id_surf, char_fill, soil_pars_f % fill, .false., 'soil_pars')
            call get_attribute(id_surf, char_lod, soil_pars_f % lod, .false., 'soil_pars')

            !--    Inquire number of soil parameters
            call get_dimension_length(id_surf, soil_pars_f % np, 'nsoil_pars')

            !--    In case of level of detail 2, also inquire number of vertical soil layers.
            if (soil_pars_f % lod == 2) then
                call get_dimension_length(id_surf, soil_pars_f % nz, 'zsoil')
            end if

            !--    Read soil parameters, depending on level of detail
            if (soil_pars_f % lod == 1) then
                allocate (soil_pars_f % pars_xy(0:soil_pars_f % np - 1, nysg:nyng, nxlg:nxrg))

                call get_variable(id_surf, 'soil_pars', soil_pars_f % pars_xy, &
                                  nxl, nxr, nys, nyn, 0, soil_pars_f % np - 1, nbgp=nbgp)

            elseif (soil_pars_f % lod == 2) then
                allocate (soil_pars_f % pars_xyz(0:soil_pars_f % np - 1, 0:soil_pars_f % nz - 1, nysg:nyng, nxlg:nxrg))
                call get_variable(id_surf, 'soil_pars', soil_pars_f % pars_xyz, &
                                  nxl, nxr, nys, nyn, 0, soil_pars_f % nz - 1, 0, soil_pars_f % np - 1, &
                                  nbgp=nbgp)
            end if
        else
            soil_pars_f % from_file = .false.
        end if

        !-- Read water parameters and related information
        if (check_existence(var_names, 'water_pars')) then
            water_pars_f % from_file = .true.
            call get_attribute(id_surf, char_fill, water_pars_f % fill, .false., 'water_pars')
            !--    Inquire number of water parameters
            call get_dimension_length(id_surf, water_pars_f % np, 'nwater_pars')
            !--    Allocate input array for water parameters and read them
            allocate (water_pars_f % pars_xy(0:water_pars_f % np - 1, nysg:nyng, nxlg:nxrg))
            call get_variable(id_surf, 'water_pars', water_pars_f % pars_xy, &
                              nxl, nxr, nys, nyn, 0, water_pars_f % np - 1, nbgp=nbgp)
        else
            water_pars_f % from_file = .false.
        end if
        !-- Read root area density - parametrized vegetation
        if (check_existence(var_names, 'root_area_dens_s')) then
            root_area_density_lsm_f % from_file = .true.
            call get_attribute(id_surf, char_fill, root_area_density_lsm_f % fill, .false., &
                               'root_area_dens_s')
            !--    Obtain number of soil layers from file and allocate variable
            call get_dimension_length(id_surf, root_area_density_lsm_f % nz, 'zsoil')
            allocate (root_area_density_lsm_f % var(0:root_area_density_lsm_f % nz - 1, nysg:nyng, nxlg:nxrg))

            !--    Read root-area density
            call get_variable(id_surf, 'root_area_dens_s', root_area_density_lsm_f % var, &
                              nxl, nxr, nys, nyn, 0, root_area_density_lsm_f % nz - 1, nbgp=nbgp)

        else
            root_area_density_lsm_f % from_file = .false.
        end if
        !-- Read street type and street crossing
        if (check_existence(var_names, 'street_type')) then
            street_type_f % from_file = .true.
            call get_attribute(id_surf, char_fill, street_type_f % fill, .false., 'street_type')
            allocate (street_type_f % var(nys:nyn, nxl:nxr))
            call get_variable(id_surf, 'street_type', street_type_f % var, nxl, nxr, nys, nyn, nbgp=0)
        else
            street_type_f % from_file = .false.
        end if

        if (check_existence(var_names, 'street_crossing')) then
            street_crossing_f % from_file = .true.
            call get_attribute(id_surf, char_fill, street_crossing_f % fill, .false., 'street_crossing')
            allocate (street_crossing_f % var(nys:nyn, nxl:nxr))
            call get_variable(id_surf, 'street_crossing', street_crossing_f % var, nxl, nxr, nys, nyn, &
                              nbgp=0)
        else
            street_crossing_f % from_file = .false.
        end if


        !-- Still missing: root_resolved and building_surface_pars.
        !-- Will be implemented as soon as they are available.
        !-- Finally, close input file
        call close_input_file(id_surf)
        ! #endif
        !-- End of CPU measurement
        ! call cpu_log(log_point_s(82), 'NetCDF input', 'stop')

        !-- Exchange ghost points for surface variables. Therefore, resize variables.
        if (albedo_type_f % from_file) then
            call exchange_horiz_2d_byte(albedo_type_f % var, nys, nyn, nxl, nxr, nbgp)
            call set_lateral_neumann_bc(albedo_type_f % var)
        end if
        if (pavement_type_f % from_file) then
            call exchange_horiz_2d_byte(pavement_type_f % var, nys, nyn, nxl, nxr, nbgp)
            call set_lateral_neumann_bc(pavement_type_f % var)
        end if
        if (soil_type_f % from_file .and. allocated(soil_type_f % var_2d)) then
            call exchange_horiz_2d_byte(soil_type_f % var_2d, nys, nyn, nxl, nxr, nbgp)
            call set_lateral_neumann_bc(soil_type_f % var_2d)
        end if
        if (vegetation_type_f % from_file) then
            call exchange_horiz_2d_byte(vegetation_type_f % var, nys, nyn, nxl, nxr, nbgp)
            call set_lateral_neumann_bc(vegetation_type_f % var)
        end if
        if (water_type_f % from_file) then
            call exchange_horiz_2d_byte(water_type_f % var, nys, nyn, nxl, nxr, nbgp)
            call set_lateral_neumann_bc(water_type_f % var)
        end if

        !-- Exchange ghost points for 3/4-D variables. For the sake of simplicity, loop further dimensions
        !-- to use 2D exchange routines. This is preferred to introducing new MPI-data types that are
        !-- required just here.
        !-- In order to avoid compiler warnings due to non-contiguous array-arguments, use temporary
        !-- 2d-arrays.
        if (soil_type_f % from_file .and. allocated(soil_type_f % var_3d)) then
            allocate (tmp_2d_byte(nysg:nyng, nxlg:nxrg))
            do k = 0, nz_soil
                tmp_2d_byte(:, :) = soil_type_f % var_3d(k, :, :)
                call exchange_horiz_2d_byte(tmp_2d_byte, nys, nyn, nxl, nxr, nbgp)
                call set_lateral_neumann_bc(tmp_2d_byte)
                soil_type_f % var_3d(k, :, :) = tmp_2d_byte(:, :)
            end do
            deallocate (tmp_2d_byte)
        end if

        allocate (tmp_2d(nysg:nyng, nxlg:nxrg))

        if (surface_fraction_f % from_file) then
            do k = 0, surface_fraction_f % nf - 1
                tmp_2d(:, :) = surface_fraction_f % frac(k, :, :)
                call exchange_horiz_2d(tmp_2d)
                call set_lateral_neumann_bc(tmp_2d)
                surface_fraction_f % frac(k, :, :) = tmp_2d(:, :)
            end do
        end if

        if (albedo_pars_f % from_file) then
            do k = 0, albedo_pars_f % np - 1
                tmp_2d(:, :) = albedo_pars_f % pars_xy(k, :, :)
                call exchange_horiz_2d(tmp_2d)
                call set_lateral_neumann_bc(tmp_2d)
                albedo_pars_f % pars_xy(k, :, :) = tmp_2d(:, :)
            end do
        end if

        if (pavement_pars_f % from_file) then
            do k = 0, pavement_pars_f % np - 1
                tmp_2d(:, :) = pavement_pars_f % pars_xy(k, :, :)
                call exchange_horiz_2d(tmp_2d)
                call set_lateral_neumann_bc(tmp_2d)
                pavement_pars_f % pars_xy(k, :, :) = tmp_2d(:, :)
            end do
        end if

        if (vegetation_pars_f % from_file) then
            do k = 0, vegetation_pars_f % np - 1
                tmp_2d(:, :) = vegetation_pars_f % pars_xy(k, :, :)
                call exchange_horiz_2d(tmp_2d)
                call set_lateral_neumann_bc(tmp_2d)
                vegetation_pars_f % pars_xy(k, :, :) = tmp_2d(:, :)
            end do
        end if

        if (water_pars_f % from_file) then
            do k = 0, water_pars_f % np - 1
                tmp_2d(:, :) = water_pars_f % pars_xy(k, :, :)
                call exchange_horiz_2d(tmp_2d)
                call set_lateral_neumann_bc(tmp_2d)
                water_pars_f % pars_xy(k, :, :) = tmp_2d(:, :)
            end do
        end if

        if (root_area_density_lsm_f % from_file) then
            do k = 0, root_area_density_lsm_f % nz - 1
                tmp_2d(:, :) = root_area_density_lsm_f % var(k, :, :)
                call exchange_horiz_2d(tmp_2d)
                call set_lateral_neumann_bc(tmp_2d)
                root_area_density_lsm_f % var(k, :, :) = tmp_2d(:, :)
            end do
        end if

        if (soil_pars_f % from_file) then
            if (soil_pars_f % lod == 1) then
                do k = 0, soil_pars_f % np - 1
                    tmp_2d(:, :) = soil_pars_f % pars_xy(k, :, :)
                    call exchange_horiz_2d(tmp_2d)
                    call set_lateral_neumann_bc(tmp_2d)
                    soil_pars_f % pars_xy(k, :, :) = tmp_2d(:, :)
                end do

            elseif (soil_pars_f % lod == 2) then
                if (myid == 0) print *, '### netcdf_data_input_surface_data for soil_pars_f%from_file #3"'
                do k2 = 0, soil_pars_f % nz - 1
                    do k = 0, soil_pars_f % np - 1
                        tmp_2d(:, :) = soil_pars_f % pars_xyz(k, k2, :, :)
                        call exchange_horiz_2d(tmp_2d)
                        call set_lateral_neumann_bc(tmp_2d)
                        soil_pars_f % pars_xyz(k, k2, :, :) = tmp_2d(:, :)
                    end do
                end do
            end if
        end if

        if (pavement_subsurface_pars_f % from_file) then
            do k2 = 0, pavement_subsurface_pars_f % nz - 1
                do k = 0, pavement_subsurface_pars_f % np - 1
                    tmp_2d(:, :) = pavement_subsurface_pars_f % pars_xyz(k, k2, :, :)
                    call exchange_horiz_2d(tmp_2d)
                    call set_lateral_neumann_bc(tmp_2d)
                    pavement_subsurface_pars_f % pars_xyz(k, k2, :, :) = tmp_2d(:, :)
                end do
            end do
        end if

        deallocate (tmp_2d)

    end subroutine netcdf_data_input_surface_data

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Reads initialization data of u, v, w, pt, q, geostrophic wind components, as well as soil
    !> moisture and soil temperature, derived from larger-scale model data.
    !--------------------------------------------------------------------------------------------------!
    subroutine netcdf_data_input_init_3d

        use arrays_3d, &
            only: pt, &
                  q, &
                  u, &
                  v, &
                  w, &
                  zu, &
                  zw

        implicit none

        character(LEN=100), dimension(:), allocatable :: var_names !< variables on dynamic input file
        character(LEN=500) :: variable_not_found = '' !< string to gather information of missing dimensions and variables

        integer(iwp) :: id_dynamic !< NetCDF id of dynamic input file
        integer(iwp) :: k !< running index along z-direction
        integer(iwp) :: n !< running index for chemistry variables
        integer(iwp) :: num_vars !< number of variables in netcdf input file

        logical :: check_passed !< flag indicating if a check passed
        logical :: dynamic_3d = .true. !< flag indicating that 3D data is read from dynamic file
        logical :: trigger_error_message !< control flag to check whether all required input data is available

        !-- Skip routine if no input file with dynamic input data is available.
        if (.not. input_pids_dynamic) return
        !-- Please note, the dynamic input file provides initial data for u and v for the prognostic grid
        !-- points in case of lateral Dirichlet conditions.
        !-- This means that the dynamic input file provides data from nxlu:nxr (for u) and from
        !-- nysv:nyn (for v) at the left and south domain boundary, respectively.
        !-- However, as work-around for the moment, PALM will run with cyclic conditions and will be
        !-- initialized with data provided by the dynamic input in case of Dirichlet boundaries.
        !-- Hence, simply set set nxlu/nysv to 1 (will be reset to its original value at the end of this
        !-- routine).
        if (bc_lr_cyc .and. nxl == 0) nxlu = 1
        if (bc_ns_cyc .and. nys == 0) nysv = 1

        !-- CPU measurement
        ! call cpu_log(log_point_s(85), 'NetCDF input init', 'start')

        ! #if defined ( __netcdf )
        !-- Open file in read-only mode
        call open_read_file(trim(input_file_dynamic)//trim(coupling_char), id_dynamic)

        !-- At first, inquire all variable names.
        call inquire_num_variables(id_dynamic, num_vars)
        !-- Allocate memory to store variable names.
        allocate (var_names(1:num_vars))
        call inquire_variable_names(id_dynamic, var_names)
        !-- Read vertical dimension of scalar und w grid.
        !-- For each dimension read, check if dimension is given in file.
        trigger_error_message = .false.
        variable_not_found = ''

        if (check_existence(var_names, 'z')) then
            call get_dimension_length(id_dynamic, init_3d % nzu, 'z')
        else
            trigger_error_message = .true.
            variable_not_found = trim('z')
        end if

        if (check_existence(var_names, 'zw')) then
            call get_dimension_length(id_dynamic, init_3d % nzw, 'zw')
        else
            trigger_error_message = .true.
            variable_not_found = trim(variable_not_found) // ', ' // 'zw'
        end if
        !-- Read also the horizontal dimensions. These are used just used for checking the compatibility
        !-- with the PALM grid before reading.
        if (check_existence(var_names, 'x')) then
            call get_dimension_length(id_dynamic, init_3d % nx, 'x')
        else
            trigger_error_message = .true.
            variable_not_found = trim(variable_not_found) // ', ' // 'x'
        end if

        if (check_existence(var_names, 'xu')) then
            call get_dimension_length(id_dynamic, init_3d % nxu, 'xu')
        else
            trigger_error_message = .true.
            variable_not_found = trim(variable_not_found) // ', ' // 'xu'
        end if

        if (check_existence(var_names, 'y')) then
            call get_dimension_length(id_dynamic, init_3d % ny, 'y')
        else
            trigger_error_message = .true.
            variable_not_found = trim(variable_not_found) // ', ' // 'y'
        end if

        if (check_existence(var_names, 'yv')) then
            call get_dimension_length(id_dynamic, init_3d % nyv, 'yv')
        else
            trigger_error_message = .true.
            variable_not_found = trim(variable_not_found) // ', ' // 'yv'
        end if

        if (trigger_error_message) then
            print *, 'initialization with data from dynamic driver - dimension(s): &"' // &
                             trim(variable_not_found) // '" not found in dynamic driver'
            ! ! call message('netcdf_data_input_mod', 'DRV0003', 1, 2, 0, 6, 0)
        end if

        !-- Check for correct horizontal and vertical dimension. Please note, checks are performed directly
        !-- here and not called from check_parameters as some varialbes are still not allocated there.
        !-- Moreover, please note, u- and v-grid has 1 grid point less than the dynamic file grid.
        if (init_3d % nx - 1 /= nx) then
            write (*, '(A,I5,A,I5,A)') &
                'number of grid points along x in dynamic input file (=', init_3d % nx - 1, &
                ') does not match the number of grid points in this run (nx=', nx, ')'
            ! ! call message('netcdf_data_input_mod', 'DRV0004', 1, 2, 0, 6, 0)
        end if
        if (init_3d % nxu - 1 /= nx - 1) then
            write (*, '(A,I5,A,I5,A)') &
                'number of grid points along x in dynamic input file (nxu=', init_3d % nxu - 1, &
                ') does not match the number of grid points in this run (nx-1=', nx - 1, ')'
            ! ! call message('netcdf_data_input_mod', 'DRV0004', 1, 2, 0, 6, 0)
        end if
        if (init_3d % ny - 1 /= ny) then
            write (*, '(A,I5,A,I5,A)') &
                'number of grid points along y in dynamic input file (=', init_3d % ny - 1, &
                ') does not match the number of grid points in this run (ny=', ny, ')'
            ! ! call message('netcdf_data_input_mod', 'DRV0004', 1, 2, 0, 6, 0)
        end if
        if (init_3d % nyv - 1 /= ny - 1) then
            write (*, '(A,I5,A,I5,A)') &
                'number of grid points along y in dynamic input file (nyv=', init_3d % nyv - 1, &
                ') does not match the number of grid points in this run (ny-1=', ny - 1, ')'
            ! ! call message('netcdf_data_input_mod', 'DRV0004', 1, 2, 0, 6, 0)
        end if

        if (init_3d % nzu /= nz) then
            write (*, '(A,I5,A,I5,A)') &
                'number of grid points along z in dynamic input file (nzu=', init_3d % nzu, &
                ') does not match the number of grid points in this run (nz=', nz, ')'
            ! ! call message('netcdf_data_input_mod', 'DRV0004', 1, 2, 0, 6, 0)
        end if
        if (init_3d % nzw /= nz - 1) then
            write (*, '(A,I5,A,I5,A)') &
                'number of grid points along z in dynamic input file (nzw=', init_3d % nzw, &
                ') does not match the number of grid points in this run (nz-1=', nz - 1, ')'
            ! ! call message('netcdf_data_input_mod', 'DRV0004', 1, 2, 0, 6, 0)
        end if
        !-- Read vertical dimensions. Later, these are required for eventual inter- and extrapolations of
        !-- the initialization data.
        if (check_existence(var_names, 'z')) then
            allocate (init_3d % zu_atmos(1:init_3d % nzu))
            call get_variable(id_dynamic, 'z', init_3d % zu_atmos)
        end if
        if (check_existence(var_names, 'zw')) then
            allocate (init_3d % zw_atmos(1:init_3d % nzw))
            call get_variable(id_dynamic, 'zw', init_3d % zw_atmos)
        end if
        !-- Check for consistency between vertical coordinates in dynamic driver and numeric grid.
        !-- Please note, depending on compiler options both may be  equal up to a certain threshold, and
        !-- differences between the numeric grid and vertical coordinate in the driver can built-up to
        !-- 10E-1-10E-0 m. For this reason, the check allows for a tolerance. In order to take into account
        !-- the general spatial resolution of the setup, the tolerance is defined as 10% of the smallest
        !-- vertical grid spacing.
        !-- Note, the check is performed separately for 'z' and 'zw'.
        if (any(abs(zu(1:nzt) - init_3d % zu_atmos(1:init_3d % nzu)) > &
                0.1_wp * minval(dz(1:number_dz)))) then
            k = 1
            do while (k <= nzt)
                if (abs(init_3d % zu_atmos(k) - zu(k)) > 0.1_wp * minval(dz(1:number_dz))) exit
                k = k + 1
            end do
            write (*, *) 'vertical coordinate z in dynamic driver does not match the ', &
                'computational grid (zu) in this run at z(', k, ')'
            ! ! call message('netcdf_data_input_mod', 'DRV0005', 1, 2, 0, 6, 0)
        end if
        if (any(abs(zw(1:nzt - 1) - init_3d % zw_atmos(1:init_3d % nzw)) > &
                0.1_wp * minval(dz(1:number_dz)))) then
            k = 1
            do while (k <= nzt)
                if (abs(init_3d % zw_atmos(k) - zw(k)) > 0.1_wp * minval(dz(1:number_dz))) exit
                k = k + 1
            end do
            write (*, *) 'vertical coordinate zw in dynamic driver does not match the ', &
                'computational grid (zw) in this run at zw(', k, ')'
            ! ! call message('netcdf_data_input_mod', 'DRV0005', 1, 2, 0, 6, 0)
        end if
        !-- Read initial geostrophic wind components at t = 0 (index 1 in file).
        if (check_existence(var_names, 'ls_forcing_ug')) then
            allocate (init_3d % ug_init(nzb:nzt + 1))
            init_3d % ug_init = 0.0_wp
            call get_variable_pr(id_dynamic, 'ls_forcing_ug', 1, init_3d % ug_init(1:nzt))
            !--    Set top-boundary condition (Neumann)
            init_3d % ug_init(nzt + 1) = init_3d % ug_init(nzt)
            init_3d % from_file_ug = .true.
        else
            init_3d % from_file_ug = .false.
        end if

        if (check_existence(var_names, 'ls_forcing_vg')) then
            allocate (init_3d % vg_init(nzb:nzt + 1))
            init_3d % vg_init = 0.0_wp
            call get_variable_pr(id_dynamic, 'ls_forcing_vg', 1, init_3d % vg_init(1:nzt))
            !--    Set top-boundary condition (Neumann)
            init_3d % vg_init(nzt + 1) = init_3d % vg_init(nzt)
            init_3d % from_file_vg = .true.
        else
            init_3d % from_file_vg = .false.
        end if
        !-- Read inital 3D data of u, v, w, pt and q. Read PE-wise yz-slices.
        !-- Please note, the u-, v- and w-component are defined on different grids with one element less in
        !-- the x-, y-, and z-direction, respectively. Hence, reading is subdivided into separate loops.
        !-- Read u-component
        if (check_existence(var_names, 'init_atmosphere_u')) then
            !--    Read attributes for the fill value and level-of-detail.
            call get_attribute(id_dynamic, char_fill, init_3d % fill_u, .false., 'init_atmosphere_u')
            call get_attribute(id_dynamic, char_lod, init_3d % lod_u, .false., 'init_atmosphere_u')
            !--    level-of-detail 1 - read initialization profile.
            if (init_3d % lod_u == 1) then
                allocate (init_3d % u_init(nzb:nzt + 1))
                init_3d % u_init = 0.0_wp
                call get_variable(id_dynamic, 'init_atmosphere_u', init_3d % u_init(nzb + 1:nzt))
                !--       Set top-boundary condition (Neumann).
                init_3d % u_init(nzt + 1) = init_3d % u_init(nzt)
            !--    level-of-detail 2 - read 3D initialization data. Before data input, u-component is
            !--    initialized with 0. This is required to have meaningful values at ghost points at the
            !--    lateral domain boundaries in case of non-cyclic boundary conditions.
            elseif (init_3d % lod_u == 2) then
                u = 0.0_wp
                call get_variable(id_dynamic, 'init_atmosphere_u', u(nzb + 1:nzt, nys:nyn, nxlu:nxr), &
                                  nxlu, nys + 1, nzb + 1, nxr - nxlu + 1, nyn - nys + 1, init_3d % nzu, dynamic_3d)
                !--       Set value at leftmost model grid point nxl = 0. This is because the dynamic file provides
                !--       data only from 1:nx-1.
                if (nxl == 0) u(nzb + 1:nzt, nys:nyn, nxl) = u(nzb + 1:nzt, nys:nyn, nxlu)
                !--       Overwrite all fill values that are found with meaningful data.
                call overwrite_fill_values(u(nzb + 1:nzt, nys:nyn, nxlu:nxr), nxlu, nxr, nys, nyn, nzb + 1, &
                                           nzt, init_3d % fill_u, 0.0_wp)
                !--       Set lateral boundary.
                if (bc_dirichlet_l .or. bc_radiation_l) u(:, :, nxlu - 1) = u(:, :, nxlu)
                if (bc_dirichlet_r .or. bc_radiation_r) u(:, :, nxr + 1) = u(:, :, nxr)
                if (bc_dirichlet_s .or. bc_radiation_s) u(:, nys - 1, :) = u(:, nys, :)
                if (bc_dirichlet_n .or. bc_radiation_n) u(:, nyn + 1, :) = u(:, nyn, :)
                !--       Set bottom and top-boundary.
                u(nzb, :, :) = u(nzb + 1, :, :)
                u(nzt + 1, :, :) = u(nzt, :, :)

            end if
            init_3d % from_file_u = .true.
        else
            print *, 'missing initial data for u-component'
            ! ! call message('netcdf_data_input_mod', 'DRV0006', 1, 2, 0, 6, 0)
        end if
        !-- Read v-component.
        if (check_existence(var_names, 'init_atmosphere_v')) then
            !--    Read attributes for the fill value and level-of-detail.
            call get_attribute(id_dynamic, char_fill, init_3d % fill_v, .false., 'init_atmosphere_v')
            call get_attribute(id_dynamic, char_lod, init_3d % lod_v, .false., 'init_atmosphere_v')
            !--    level-of-detail 1 - read initialization profile.
            if (init_3d % lod_v == 1) then
                allocate (init_3d % v_init(nzb:nzt + 1))
                init_3d % v_init = 0.0_wp
                call get_variable(id_dynamic, 'init_atmosphere_v', init_3d % v_init(nzb + 1:nzt))
                !--       Set top-boundary condition (Neumann).
                init_3d % v_init(nzt + 1) = init_3d % v_init(nzt)
            !--    level-of-detail 2 - read 3D initialization data. Before data input, v-component is
            !--    initialized with 0. This is required to have meaningful values at ghost points at the
            !--    lateral domain boundaries in case of non-cyclic boundary conditions.
            elseif (init_3d % lod_v == 2) then
                v = 0.0_wp
                call get_variable(id_dynamic, 'init_atmosphere_v', v(nzb + 1:nzt, nysv:nyn, nxl:nxr), &
                                  nxl + 1, nysv, nzb + 1, nxr - nxl + 1, nyn - nysv + 1, init_3d % nzu, dynamic_3d)
                !--       Set value at southmost model grid point nys = 0. This is because the dynamic file
                !--       provides data only from 1:ny-1.
                if (nys == 0) v(nzb + 1:nzt, nys, nxl:nxr) = v(nzb + 1:nzt, nysv, nxl:nxr)
                !--       Overwrite all fill values that are found with meaningful data.
                call overwrite_fill_values(v(nzb + 1:nzt, nysv:nyn, nxl:nxr), nxl, nxr, nysv, nyn, nzb + 1, &
                                           nzt, init_3d % fill_v, 0.0_wp)
                !--       Set lateral boundary.
                if (bc_dirichlet_l .or. bc_radiation_l) v(:, :, nxl - 1) = v(:, :, nxl)
                if (bc_dirichlet_r .or. bc_radiation_r) v(:, :, nxr + 1) = v(:, :, nxr)
                if (bc_dirichlet_s .or. bc_radiation_s) v(:, nysv - 1, :) = v(:, nysv, :)
                if (bc_dirichlet_n .or. bc_radiation_n) v(:, nyn + 1, :) = v(:, nyn, :)
                !--       Set bottom and top-boundary.
                v(nzb, :, :) = v(nzb + 1, :, :)
                v(nzt + 1, :, :) = v(nzt, :, :)
            end if
            init_3d % from_file_v = .true.
        else
            print *, 'missing initial data for v-component'
            ! ! call message('netcdf_data_input_mod', 'DRV0006', 1, 2, 0, 6, 0)
        end if
        !-- Read w-component
        if (check_existence(var_names, 'init_atmosphere_w')) then
            !--    Read attributes for the fill value and level-of-detail.
            call get_attribute(id_dynamic, char_fill, init_3d % fill_w, .false., 'init_atmosphere_w')
            call get_attribute(id_dynamic, char_lod, init_3d % lod_w, .false., 'init_atmosphere_w')
            !--    level-of-detail 1 - read initialization profile.
            if (init_3d % lod_w == 1) then
                allocate (init_3d % w_init(nzb:nzt + 1))
                init_3d % w_init = 0.0_wp
                call get_variable(id_dynamic, 'init_atmosphere_w', init_3d % w_init(nzb + 1:nzt - 1))
                !--       Set top-boundary condition (Neumann).
                init_3d % w_init(nzt:nzt + 1) = init_3d % w_init(nzt - 1)
            !--    level-of-detail 2 - read 3D initialization data. Before data input, w-component is
            !--    initialized with 0. This is required to have meaningful values at ghost points at the
            !--    lateral domain boundaries in case of non-cyclic boundary conditions.
            elseif (init_3d % lod_w == 2) then
                w = 0.0_wp
                call get_variable(id_dynamic, 'init_atmosphere_w', w(nzb + 1:nzt - 1, nys:nyn, nxl:nxr), &
                                  nxl + 1, nys + 1, nzb + 1, nxr - nxl + 1, nyn - nys + 1, init_3d % nzw, dynamic_3d)
                !--       Overwrite all fill values that are found with meaningful data.
                call overwrite_fill_values(w(nzb + 1:nzt, nys:nyn, nxl:nxr), nxl, nxr, nys, nyn, nzb + 1, &
                                           nzt, init_3d % fill_w, 0.0_wp)
                !--       Set lateral boundary.
                if (bc_dirichlet_l .or. bc_radiation_l) w(:, :, nxl - 1) = w(:, :, nxl)
                if (bc_dirichlet_r .or. bc_radiation_r) w(:, :, nxr + 1) = w(:, :, nxr)
                if (bc_dirichlet_s .or. bc_radiation_s) w(:, nys - 1, :) = w(:, nys, :)
                if (bc_dirichlet_n .or. bc_radiation_n) w(:, nyn + 1, :) = w(:, nyn, :)
                !--       Set bottom and top-boundary.
                w(nzb, :, :) = 0.0_wp
                w(nzt, :, :) = w(nzt - 1, :, :)
                w(nzt + 1, :, :) = w(nzt - 1, :, :)
            end if
            init_3d % from_file_w = .true.
        else
            print *, 'missing initial data for w-component'
            ! ! call message('netcdf_data_input_mod', 'DRV0006', 1, 2, 0, 6, 0)
        end if
        !-- Read potential temperature
        if (.not. neutral) then
            if (check_existence(var_names, 'init_atmosphere_pt')) then
                !--       Read attributes for the fill value and level-of-detail
                call get_attribute(id_dynamic, char_fill, init_3d % fill_pt, .false., 'init_atmosphere_pt')
                call get_attribute(id_dynamic, char_lod, init_3d % lod_pt, .false., 'init_atmosphere_pt')
                !--       level-of-detail 1 - read initialization profile
                if (init_3d % lod_pt == 1) then
                    allocate (init_3d % pt_init(nzb:nzt + 1))
                    call get_variable(id_dynamic, 'init_atmosphere_pt', init_3d % pt_init(nzb + 1:nzt))
                    !--          Set Neumann top and surface boundary condition for initial profil
                    init_3d % pt_init(nzb) = init_3d % pt_init(nzb + 1)
                    init_3d % pt_init(nzt + 1) = init_3d % pt_init(nzt)
                !--       level-of-detail 2 - read 3D initialization data
                elseif (init_3d % lod_pt == 2) then

                    call get_variable(id_dynamic, 'init_atmosphere_pt', pt(nzb + 1:nzt, nys:nyn, nxl:nxr), &
                                      nxl + 1, nys + 1, nzb + 1, nxr - nxl + 1, nyn - nys + 1, init_3d % nzu, dynamic_3d)
                    !--          Overwrite all fill values that are found with meaningful data.
                    call overwrite_fill_values(pt(nzb + 1:nzt, nys:nyn, nxl:nxr), nxl, nxr, nys, nyn, nzb + 1, &
                                               nzt, init_3d % fill_pt, pt_surface)
                    !--          Set lateral boundary.
                    if (bc_dirichlet_l .or. bc_radiation_l) pt(:, :, nxl - 1) = pt(:, :, nxl)
                    if (bc_dirichlet_r .or. bc_radiation_r) pt(:, :, nxr + 1) = pt(:, :, nxr)
                    if (bc_dirichlet_s .or. bc_radiation_s) pt(:, nys - 1, :) = pt(:, nys, :)
                    if (bc_dirichlet_n .or. bc_radiation_n) pt(:, nyn + 1, :) = pt(:, nyn, :)
                    !--          Set bottom and top-boundary
                    pt(nzb, :, :) = pt(nzb + 1, :, :)
                    pt(nzt + 1, :, :) = pt(nzt, :, :)
                end if
                init_3d % from_file_pt = .true.
            else
                print *, 'missing initial data for potential temperature'
                ! ! call message('netcdf_data_input_mod', 'DRV0006', 1, 2, 0, 6, 0)
            end if
        end if

        !-- Read mixing ratio
        if (humidity) then
            if (check_existence(var_names, 'init_atmosphere_qv')) then
                !--       Read attributes for the fill value and level-of-detail
                call get_attribute(id_dynamic, char_fill, init_3d % fill_q, .false., 'init_atmosphere_qv')
                call get_attribute(id_dynamic, char_lod, init_3d % lod_q, .false., 'init_atmosphere_qv')
                !--       level-of-detail 1 - read initialization profile
                if (init_3d % lod_q == 1) then
                    allocate (init_3d % q_init(nzb:nzt + 1))
                    call get_variable(id_dynamic, 'init_atmosphere_qv', init_3d % q_init(nzb + 1:nzt))
                    !--          Set bottom and top boundary condition (Neumann)
                    init_3d % q_init(nzb) = init_3d % q_init(nzb + 1)
                    init_3d % q_init(nzt + 1) = init_3d % q_init(nzt)
                !--       level-of-detail 2 - read 3D initialization data
                elseif (init_3d % lod_q == 2) then

                    call get_variable(id_dynamic, 'init_atmosphere_qv', q(nzb + 1:nzt, nys:nyn, nxl:nxr), &
                                      nxl + 1, nys + 1, nzb + 1, nxr - nxl + 1, nyn - nys + 1, init_3d % nzu, dynamic_3d)
                    !--          Overwrite all fill values that are found with meaningful data.
                    call overwrite_fill_values(q(nzb + 1:nzt, nys:nyn, nxl:nxr), nxl, nxr, nys, nyn, nzb + 1, &
                                               nzt, init_3d % fill_q, q_surface)
                    !--          Set lateral boundary.
                    if (bc_dirichlet_l .or. bc_radiation_l) q(:, :, nxl - 1) = q(:, :, nxl)
                    if (bc_dirichlet_r .or. bc_radiation_r) q(:, :, nxr + 1) = q(:, :, nxr)
                    if (bc_dirichlet_s .or. bc_radiation_s) q(:, nys - 1, :) = q(:, nys, :)
                    if (bc_dirichlet_n .or. bc_radiation_n) q(:, nyn + 1, :) = q(:, nyn, :)
                    !--          Set bottom and top-boundary
                    q(nzb, :, :) = q(nzb + 1, :, :)
                    q(nzt + 1, :, :) = q(nzt, :, :)
                end if
                init_3d % from_file_q = .true.
            else
                print *, 'missing initial data for mixing ratio'
                ! ! call message('netcdf_data_input_mod', 'DRV0006', 1, 2, 0, 6, 0)
            end if
        end if

        !-- Read chemistry variables.
        !-- Please note, for the moment, only LOD=1 is allowed
        if (air_chemistry) then
            !--    Allocate chemistry input profiles, as well as arrays for fill values and LOD's.
            allocate (init_3d % chem_init(nzb:nzt + 1, 1:ubound(init_3d % var_names_chem, 1)))
            allocate (init_3d % fill_chem(1:ubound(init_3d % var_names_chem, 1)))
            allocate (init_3d % lod_chem(1:ubound(init_3d % var_names_chem, 1)))

            do n = 1, ubound(init_3d % var_names_chem, 1)
                if (check_existence(var_names, trim(init_3d % var_names_chem(n)))) then
                    !--          Read attributes for the fill value and level-of-detail
                    call get_attribute(id_dynamic, char_fill, init_3d % fill_chem(n), .false., &
                                       trim(init_3d % var_names_chem(n)))
                    call get_attribute(id_dynamic, char_lod, init_3d % lod_chem(n), .false., &
                                       trim(init_3d % var_names_chem(n)))
                    !--          Give message that only LOD=1 is allowed.
                    if (init_3d % lod_chem(n) /= 1) then
                        print *, 'for chemistry variables only LOD=1 is allowed'
                        ! ! call message('netcdf_data_input_mod', 'DRV0007', 1, 2, 0, 6, 0)
                    end if
                    !--          level-of-detail 1 - read initialization profile
                    call get_variable(id_dynamic, trim(init_3d % var_names_chem(n)), &
                                      init_3d % chem_init(nzb + 1:nzt, n))
                    !--          Set bottom and top boundary condition (Neumann)
                    init_3d % chem_init(nzb, n) = init_3d % chem_init(nzb + 1, n)
                    init_3d % chem_init(nzt + 1, n) = init_3d % chem_init(nzt, n)
                    init_3d % from_file_chem(n) = .true.
                end if
            end do
        end if
        !-- Close input file
        call close_input_file(id_dynamic)
        ! #endif
        !-- End of CPU measurement
        ! call cpu_log(log_point_s(85), 'NetCDF input init', 'stop')
        !-- Finally, check if the input data has any fill values. Please note, checks depend on the LOD of
        !-- the input data.
        if (init_3d % from_file_u) then
            check_passed = .true.
            if (init_3d % lod_u == 1) then
                if (any(init_3d % u_init(nzb + 1:nzt + 1) == init_3d % fill_u)) check_passed = .false.
            elseif (init_3d % lod_u == 2) then
                if (any(u(nzb + 1:nzt + 1, nys:nyn, nxlu:nxr) == init_3d % fill_u)) check_passed = .false.
            end if
            if (.not. check_passed) then
                print *, 'netCDF input for init_atmosphere_u must not contain any _FillValues'
                ! ! call message('netcdf_data_input_mod', 'DRV0008', 2, 2, 0, 6, 0)
            end if
        end if

        if (init_3d % from_file_v) then
            check_passed = .true.
            if (init_3d % lod_v == 1) then
                if (any(init_3d % v_init(nzb + 1:nzt + 1) == init_3d % fill_v)) check_passed = .false.
            elseif (init_3d % lod_v == 2) then
                if (any(v(nzb + 1:nzt + 1, nysv:nyn, nxl:nxr) == init_3d % fill_v)) check_passed = .false.
            end if
            if (.not. check_passed) then
                print *, 'netCDF input for init_atmosphere_v must not contain any _FillValues'
                ! ! call message('netcdf_data_input_mod', 'DRV0008', 2, 2, 0, 6, 0)
            end if
        end if

        if (init_3d % from_file_w) then
            check_passed = .true.
            if (init_3d % lod_w == 1) then
                if (any(init_3d % w_init(nzb + 1:nzt) == init_3d % fill_w)) check_passed = .false.
            elseif (init_3d % lod_w == 2) then
                if (any(w(nzb + 1:nzt, nys:nyn, nxl:nxr) == init_3d % fill_w)) check_passed = .false.
            end if
            if (.not. check_passed) then
                print *, 'NetCDF input for init_atmosphere_w must not contain any _FillValues'
                ! ! call message('netcdf_data_input_mod', 'DRV0008', 2, 2, 0, 6, 0)
            end if
        end if

        if (init_3d % from_file_pt) then
            check_passed = .true.
            if (init_3d % lod_pt == 1) then
                if (any(init_3d % pt_init(nzb + 1:nzt + 1) == init_3d % fill_pt)) check_passed = .false.
            elseif (init_3d % lod_pt == 2) then
                if (any(pt(nzb + 1:nzt + 1, nys:nyn, nxl:nxr) == init_3d % fill_pt)) check_passed = .false.
            end if
            if (.not. check_passed) then
                print *, 'NetCDF input for init_atmosphere_pt must not contain any _FillValues'
                ! ! call message('netcdf_data_input_mod', 'DRV0008', 2, 2, 0, 6, 0)
            end if
        end if

        if (init_3d % from_file_q) then
            check_passed = .true.
            if (init_3d % lod_q == 1) then
                if (any(init_3d % q_init(nzb + 1:nzt + 1) == init_3d % fill_q)) check_passed = .false.
            elseif (init_3d % lod_q == 2) then
                if (any(q(nzb + 1:nzt + 1, nys:nyn, nxl:nxr) == init_3d % fill_q)) check_passed = .false.
            end if
            if (.not. check_passed) then
                print *, 'NetCDF input for init_atmosphere_q must not contain any _FillValues'
                ! ! call message('netcdf_data_input_mod', 'DRV0008', 2, 2, 0, 6, 0)
            end if
        end if
        !-- Set ghost boundaries.
        call exchange_horiz(u, nbgp)
        call exchange_horiz(v, nbgp)
        call exchange_horiz(w, nbgp)
        if (.not. neutral) call exchange_horiz(pt, nbgp)
        if (humidity) call exchange_horiz(q, nbgp)

        !-- Workaround for cyclic conditions. Please see above for further explanation.
        if (bc_lr_cyc .and. nxl == 0) nxlu = nxl
        if (bc_ns_cyc .and. nys == 0) nysv = nys

    end subroutine netcdf_data_input_init_3d

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Checks input file for consistency and minimum requirements.
    !--------------------------------------------------------------------------------------------------!
    subroutine netcdf_data_input_check_dynamic

        implicit none
        !-- Dynamic input file must be also present if initialization via read_from_file is prescribed.
        if (.not. input_pids_dynamic .and. index(initializing_actions, 'read_from_file') /= 0) &
            then
            print *, 'initializing_actions = "read_from_file" requires dynamic driver file "' // &
                             trim(input_file_dynamic) // trim(coupling_char) // '"'
            ! call message('netcdf_data_input_mod', 'DRV0009', 1, 2, 0, 6, 0)
        end if

    end subroutine netcdf_data_input_check_dynamic

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Checks input file for consistency and minimum requirements.
    !--------------------------------------------------------------------------------------------------!
    subroutine netcdf_data_input_check_static

        implicit none

        integer(iwp) :: i !< loop index along x-direction
        integer(iwp) :: j !< loop index along y-direction
        integer(iwp) :: n_surf !< number of different surface types at given location

        logical :: check_passed !< flag indicating if a check passed

        !-- Return if no static input file is available
        if (.not. input_pids_static) return
        !-- Check for a proper valid range of global attributes. No checks for origin_x and origin_y
        !-- so far, because these values refer to the used UTM reference zone, so that the values can
        !-- be even negative. Except for the virtual measurement module, there are also no other side
        !-- effects. The same applies for origin_z, which only sets a reference height. Since PALM uses
        !-- a relative height coordinate, there is no impact on the numerical solution (expect it is not
        !-- NaN, which has been already checked before).
        !-- Longitude and latitude need to be within a geographically reasonable range.
        if (abs(input_file_atts % origin_lon) > 180.0_wp .or. &
            abs(input_file_atts % origin_lat) > 90.0_wp) then
            write (*, *) 'static driver coordinates "origin_lon" = ', &
                input_file_atts % origin_lon, ' and/or "origin_lat" = ', &
                input_file_atts % origin_lat, &
                ' are outside the allowed range of values'
            ! call message('netcdf_data_input_init', 'DRV0001', 1, 2, 0, 6, 0)
        end if

        if (input_file_atts % rotation_angle < 0.0_wp .or. &
            input_file_atts % rotation_angle > 360.0_wp) then
            print *, 'static driver attribute "rotation_angle" ' // &
                             'is outside the allowed range of 0-360 degrees'
            ! call message('netcdf_data_input_init', 'DRV0041', 1, 2, 0, 6, 0)
        end if
        !-- Now, check for proper settings of specific input dimensions and variables.
        !-- Check for correct dimension of surface_fractions, should run from 0-2.
        if (surface_fraction_f % from_file) then
            if (surface_fraction_f % nf - 1 > 2) then
                write (*, *) 'nsurface_fraction = ', surface_fraction_f % nf, &
                    ' must not be larger than 3'
                ! call message('netcdf_data_input_mod', 'DRV0010', 1, 2, 0, 6, 0)
            end if
        end if
        !-- Skip further checks concerning buildings and natural surface properties if no urban surface and
        !-- land surface model are applied.
        if (.not. land_surface .and. .not. urban_surface) return
        !-- Check for minimum requirement of surface-classification data in case static input file is used.
        if ((.not. vegetation_type_f % from_file .or. &
             .not. pavement_type_f % from_file .or. &
             .not. water_type_f % from_file .or. &
             .not. soil_type_f % from_file) .or. &
            (urban_surface .and. .not. building_type_f % from_file)) then
            print *, 'minimum requirement for surface classification is not fulfilled'
            ! call message('netcdf_data_input_mod', 'DRV0011', 1, 2, 0, 6, 0)
        end if
        !-- Check for general availability of input variables.
        !-- If vegetation_type is 0 at any location, vegetation_pars as well as root_area_dens_s are
        !-- required.
        if (vegetation_type_f % from_file) then
            if (any(vegetation_type_f % var == 0)) then
                if (.not. vegetation_pars_f % from_file) then
                    print *, 'if vegetation_type = 0 at any location, ' // &
                                     'vegetation_pars is required'
                    ! call message('netcdf_data_input_mod', 'DRV0012', 2, 2, myid, 6, 0)
                end if
                if (.not. root_area_density_lsm_f % from_file) then
                    print *, 'if vegetation_type = 0 at any location, ' // &
                                     'root_area_dens_s is required'
                    ! call message('netcdf_data_input_mod', 'DRV0012', 2, 2, myid, 6, 0)
                end if
            end if
        end if
        !-- If soil_type is zero at any location, soil_pars is required.
        if (soil_type_f % from_file) then
            check_passed = .true.
            if (allocated(soil_type_f % var_2d)) then
                if (any(soil_type_f % var_2d == 0)) then
                    if (.not. soil_pars_f % from_file) check_passed = .false.
                end if
            else
                if (any(soil_type_f % var_3d == 0)) then
                    if (.not. soil_pars_f % from_file) check_passed = .false.
                end if
            end if
            if (.not. check_passed) then
                print *, 'if soil_type = 0 at any location, soil_pars is required'
                ! call message('netcdf_data_input_mod', 'DRV0013', 2, 2, myid, 6, 0)
            end if
        end if
        !-- Buildings require a type in case of urban-surface model.
        if (buildings_f % from_file .and. .not. building_type_f % from_file) then
            print *, 'if buildings are provided, also building_type is required'
            ! call message('netcdf_data_input_mod', 'DRV0014', 2, 2, 0, 6, 0)
        end if
        !-- Buildings require an ID.
        if (buildings_f % from_file .and. .not. building_id_f % from_file) then
            print *, 'if buildings are provided, also building_id is required'
            ! call message('netcdf_data_input_mod', 'DRV0014', 2, 2, 0, 6, 0)
        end if
        !-- If building_type is provided, also building_id is needed (due to the filtering algorithm).
        if (building_type_f % from_file .and. .not. building_id_f % from_file) then
            print *, 'if building_type is provided, also building_id is required'
            ! call message('netcdf_data_input_mod', 'DRV0016', 2, 2, 0, 6, 0)
        end if
        !-- If albedo_type is zero at any location, albedo_pars is required.
        if (albedo_type_f % from_file) then
            if (any(albedo_type_f % var == 0)) then
                if (.not. albedo_pars_f % from_file) then
                    print *, 'if albedo_type = 0 at any location, albedo_pars is required'
                    ! call message('netcdf_data_input_mod', 'DRV0017', 2, 2, myid, 6, 0)
                end if
            end if
        end if
        !-- If pavement_type is zero at any location, pavement_pars is required.
        if (pavement_type_f % from_file) then
            if (any(pavement_type_f % var == 0)) then
                if (.not. pavement_pars_f % from_file) then
                    print *, 'if pavement_type = 0 at any location, pavement_pars is required'
                    ! call message('netcdf_data_input_mod', 'DRV0018', 2, 2, myid, 6, 0)
                end if
            end if
        end if
        !-- If pavement_type is zero at any location, also pavement_subsurface_pars is required.
        if (pavement_type_f % from_file) then
            if (any(pavement_type_f % var == 0)) then
                if (.not. pavement_subsurface_pars_f % from_file) then
                    print *, 'if pavement_type = 0 at any location, ' // &
                                     'pavement_subsurface_pars is required'
                    ! call message('netcdf_data_input_mod', 'DRV0019', 2, 2, myid, 6, 0)
                end if
            end if
        end if
        !-- If water_type is zero at any location, water_pars is required.
        if (water_type_f % from_file) then
            if (any(water_type_f % var == 0)) then
                if (.not. water_pars_f % from_file) then
                    print *, 'if water_type = 0 at any location, water_pars is required'
                    ! call message('netcdf_data_input_mod', 'DRV0020', 2, 2, myid, 6, 0)
                end if
            end if
        end if
        !-- Check for local consistency of the input data.
        if (.not. cut_cell_topography) then
            do i = nxl, nxr
                do j = nys, nyn
                    !--          For each (y,x)-location at least one of the parameters vegetation_type, pavement_type,
                    !--          building_type, or water_type must be set to a non­missing value.
                    if (land_surface .and. .not. urban_surface) then
                        if (vegetation_type_f % var(j, i) == vegetation_type_f % fill .and. &
                            pavement_type_f % var(j, i) == pavement_type_f % fill .and. &
                            water_type_f % var(j, i) == water_type_f % fill) then
                            write (*, *) 'at least one of the parameters ' // &
                                'vegetation_type, pavement_type, ' // &
                                'or water_type must be set ' // &
                                'to a non-missing value at grid point (i,j) = (', i, &
                                ',', j, ')'
                            ! call message('netcdf_data_input_mod', 'DRV0021', 2, 2, myid, 6, 0)
                        end if
                    elseif (land_surface .and. urban_surface) then
                        if (vegetation_type_f % var(j, i) == vegetation_type_f % fill .and. &
                            pavement_type_f % var(j, i) == pavement_type_f % fill .and. &
                            building_type_f % var(j, i) == building_type_f % fill .and. &
                            water_type_f % var(j, i) == water_type_f % fill) then
                            write (*, *) 'at least one of the parameters ' // &
                                'vegetation_type, pavement_type, ' // &
                                'building_type, or water_type must be set ' // &
                                'to a non-missing value at grid point (i,j) = (', i, &
                                ',', j, ')'
                            ! call message('netcdf_data_input_mod', 'DRV0022', 2, 2, myid, 6, 0)
                        end if
                    end if

                    !--          Note that a soil_type is required for each location (y,x) where either vegetation_type
                    !--          or pavement_type is a non-­missing value.
                    if ((vegetation_type_f % var(j, i) /= vegetation_type_f % fill .or. &
                         pavement_type_f % var(j, i) /= pavement_type_f % fill)) then
                        check_passed = .true.
                        if (allocated(soil_type_f % var_2d)) then
                            if (soil_type_f % var_2d(j, i) == soil_type_f % fill) check_passed = .false.
                        else
                            if (any(soil_type_f % var_3d(:, j, i) == soil_type_f % fill)) check_passed = .false.
                        end if

                        if (.not. check_passed) then
                            print *, 'soil_type is required for each ' // &
                                             'location (x,y) where vegetation_type or ' // &
                                             'pavement_type is a non-missing value'
                            ! call message('netcdf_data_input_mod', 'DRV0023', 2, 2, myid, 6, 0)
                        end if
                    end if
                    !--          Check for consistency of given types. At the moment, only one of vegetation, pavement, or
                    !--          water-type can be set. This is because no tile approach is yet implemented in the
                    !--          land-surface model. Later, when this is possible, surface fraction needs to be given and
                    !--          the sum must not be larger than 1. Please note, in case more than one type is given at a
                    !--          pixel, an error message will be given.
                    n_surf = 0
                    if (vegetation_type_f % var(j, i) /= vegetation_type_f % fill) n_surf = n_surf + 1
                    if (water_type_f % var(j, i) /= water_type_f % fill) n_surf = n_surf + 1
                    if (pavement_type_f % var(j, i) /= pavement_type_f % fill) n_surf = n_surf + 1

                    if (n_surf > 1) then
                        write (*, *) 'more than one surface type (vegetation, pavement, water)', &
                            'is given at a location (i,j) = (', i, ',', j, ')'
                        ! call message('netcdf_data_input_mod', 'DRV0024', 2, 2, myid, 6, 0)

                    !                 IF ( .NOT. surface_fraction_f%from_file )  THEN
                    !                    print *, 'More than one surface type (vegetation '//                     &
                    !                                     'pavement, water) is given at a location. '//                   &
                    !                                     'Please note, this is not possible at ' //                      &
                    !                                     'the moment as no tile approach is yet ' //                     &
                    !                                     'implemented.'
                    !                    print *, 'If more than one surface type is ' //                          &
                    !                                     'given at a location, surface_fraction ' //                     &
                    !                                     'must be provided.'
                    !                    ! call message( 'netcdf_data_input_mod', 'DRVxxx', 2, 2, myid, 6, 0 )
                    !                 ELSEIF ( ANY ( surface_fraction_f%frac(:,j,i) ==  surface_fraction_f%fill ) )  THEN
                    !                    print *, 'If more than one surface type is ' //                          &
                    !                                     'given at a location, surface_fraction ' //                     &
                    !                                     'must be provided.'
                    !                    ! call message( 'netcdf_data_input_mod', 'DRVxxxx', 2, 2, myid, 6, 0 )
                    !                 ENDIF
                    end if
                    !--          Check for further mismatches. e.g. relative fractions exceed 1 or vegetation_type is set
                    !--          but surface vegetation fraction is zero, etc..
                    if (surface_fraction_f % from_file) then
                        !--             If surface fractions is given, also check that only one type is given.
                        if (sum(merge(1, 0, surface_fraction_f % frac(:, j, i) /= 0.0_wp .and. &
                                      surface_fraction_f % frac(:, j, i) /= surface_fraction_f % fill) &
                                ) > 1) &
                            then
                            write (*, *) 'surface_fraction is given for more than one type ', &
                                'at location (i,j) = (', i, ',', j, ')'
                            ! call message('netcdf_data_input_mod', 'DRV0025', 2, 2, myid, 6, 0)
                        end if
                        !--             Sum of relative fractions must be 1. Note, attributed to type conversions due to
                        !--             reading, the sum of surface fractions might be not exactly 1. Hence, the sum is check
                        !--             with a tolerance. Later, in the land-surface model, the relative fractions are
                        !--             normalized to one. Actually, surface fractions shall be _FillValue at building grid
                        !--             points, however, in order to relax this requirement and allow that surface-fraction can
                        !--             also be zero at these grid points, only perform this check at locations where some
                        !--             vegetation, pavement or water is defined.
                        if (vegetation_type_f % var(j, i) /= vegetation_type_f % fill .or. &
                            pavement_type_f % var(j, i) /= pavement_type_f % fill .or. &
                            water_type_f % var(j, i) /= water_type_f % fill) &
                            then
                            if (sum(surface_fraction_f % frac(0:2, j, i)) > 1.0_wp + 1e-8_wp .or. &
                                sum(surface_fraction_f % frac(0:2, j, i)) < 1.0_wp - 1e-8_wp) &
                                then
                                write (*, *) 'sum of all land-surface fractions does not ', &
                                    'equal 1 at location (i,j) = (', i, ',', j, ')'
                                ! call message('netcdf_data_input_mod', 'DRV0026', 2, 2, myid, 6, 0)
                            end if
                        end if
                        !--             Relative fraction for a type must not be zero at locations where this type is set.
                        if ((vegetation_type_f % var(j, i) /= vegetation_type_f % fill .and. &
                             (surface_fraction_f % frac(0, j, i) == 0.0_wp .or. &
                              surface_fraction_f % frac(0, j, i) == surface_fraction_f % fill) &
                             ) .or. &
                            (pavement_type_f % var(j, i) /= pavement_type_f % fill .and. &
                             (surface_fraction_f % frac(1, j, i) == 0.0_wp .or. &
                              surface_fraction_f % frac(1, j, i) == surface_fraction_f % fill) &
                             ) .or. &
                            (water_type_f % var(j, i) /= water_type_f % fill .and. &
                             (surface_fraction_f % frac(2, j, i) == 0.0_wp .or. &
                              surface_fraction_f % frac(2, j, i) == surface_fraction_f % fill) &
                             )) &
                            then
                            write (*, *) 'mismatch in setting of surface_fraction at ', &
                                'location (i,j) = (', i, ',', j, ')'
                            ! call message('netcdf_data_input_mod', 'DRV0027', 2, 2, myid, 6, 0)
                        end if
                        !--             Relative fraction for a type must not contain non-zero values if this type is not
                        !--             set.
                        if ((vegetation_type_f % var(j, i) == vegetation_type_f % fill .and. &
                             (surface_fraction_f % frac(0, j, i) /= 0.0_wp .and. &
                              surface_fraction_f % frac(0, j, i) /= surface_fraction_f % fill) &
                             ) .or. &
                            (pavement_type_f % var(j, i) == pavement_type_f % fill .and. &
                             (surface_fraction_f % frac(1, j, i) /= 0.0_wp .and. &
                              surface_fraction_f % frac(1, j, i) /= surface_fraction_f % fill) &
                             ) .or. &
                            (water_type_f % var(j, i) == water_type_f % fill .and. &
                             (surface_fraction_f % frac(2, j, i) /= 0.0_wp .and. &
                              surface_fraction_f % frac(2, j, i) /= surface_fraction_f % fill) &
                             )) &
                            then
                            write (*, *) 'mismatch in setting of surface_fraction at ', &
                                'location (i,j) = (', i, ',', j, ')'
                            ! call message('netcdf_data_input_mod', 'DRV0028', 2, 2, myid, 6, 0)
                        end if
                    end if
                    !--          Check vegetation_pars. If vegetation_type is 0, all parameters need to be set, otherwise,
                    !--          single parameters set by vegetation_type can be overwritten.
                    if (vegetation_type_f % from_file) then
                        if (vegetation_type_f % var(j, i) == 0) then
                            if (any(vegetation_pars_f % pars_xy(:, j, i) == vegetation_pars_f % fill)) then
                                write (*, *) 'vegetation_pars missing at location (i,j) = (', &
                                    i, ',', j, ')'
                                ! call message('netcdf_data_input_mod', 'DRV0029', 2, 2, myid, 6, 0)
                            end if
                        end if
                    end if
                    !--          Check root distribution. If vegetation_type is 0, all levels must be set.
                    if (vegetation_type_f % from_file) then
                        if (vegetation_type_f % var(j, i) == 0) then
                            if (any(root_area_density_lsm_f % var(:, j, i) == root_area_density_lsm_f % fill)) &
                                then
                                write (*, *) 'root_area_dens_s missing at location (i,j) = (', &
                                    i, ',', j, ')'
                                ! call message('netcdf_data_input_mod', 'DRV0030', 2, 2, myid, 6, 0)
                            end if
                        end if
                    end if
                    !--          Check soil parameters. If soil_type is 0, all parameters must be set.
                    if (soil_type_f % from_file) then
                        check_passed = .true.
                        if (allocated(soil_type_f % var_2d)) then
                            if (soil_type_f % var_2d(j, i) == 0) then
                                if (any(soil_pars_f % pars_xy(:, j, i) == soil_pars_f % fill)) &
                                    check_passed = .false.
                            end if
                        else
                            if (any(soil_type_f % var_3d(:, j, i) == 0)) then
                                if (any(soil_pars_f % pars_xy(:, j, i) == soil_pars_f % fill)) &
                                    check_passed = .false.
                            end if
                        end if
                        if (.not. check_passed) then
                            write (*, *) 'soil_pars missing at location (i,j) = (', &
                                i, ',', j, ')'
                            ! call message('netcdf_data_input_mod', 'DRV0031', 2, 2, myid, 6, 0)
                        end if
                    end if
                    !--          Check if building_type is set at each building and vice versa.
                    !--          Please note, buildings are already processed and filtered.
                    !--          For this reason, consistency checks are based on topo_flags rather than
                    !--          buildings_f (buildings are represented by bit 6 in topo_flags).
                    if (building_type_f % from_file .and. buildings_f % from_file) then
                        if (any(btest(topo_flags(:, j, i), 6)) .and. &
                            building_type_f % var(j, i) == building_type_f % fill .or. &
                            .not. any(btest(topo_flags(:, j, i), 6)) .and. &
                            building_type_f % var(j, i) /= building_type_f % fill) &
                            then
                            write (*, *) 'building_type missing at location (i,j) = (', &
                                i, ',', j, ')'
                            ! call message('netcdf_data_input_mod', 'DRV0033', 2, 2, myid, 6, 0)
                        end if
                    end if
                    !--          Check if at each location where a building is present also an ID is set and vice versa.
                    if (buildings_f % from_file) then
                        if (any(btest(topo_flags(:, j, i), 6)) .and. &
                            building_id_f % var(j, i) == building_id_f % fill .or. &
                            .not. any(btest(topo_flags(:, j, i), 6)) .and. &
                            building_id_f % var(j, i) /= building_id_f % fill) &
                            then
                            write (*, *) 'building_id missing at location (i,j) = (', &
                                i, ',', j, ')'
                            ! call message('netcdf_data_input_mod', 'DRV0034', 2, 2, myid, 6, 0)
                        end if
                    end if
                    !--          Check if building ID is set where a bulding is defined.
                    if (buildings_f % from_file) then
                        if (any(btest(topo_flags(:, j, i), 6)) .and. &
                            building_id_f % var(j, i) == building_id_f % fill) then
                            write (*, *) 'building_id missing at location (i,j) = (', &
                                i, ',', j, ')'
                            ! call message('netcdf_data_input_mod', 'DRV0034', 2, 2, myid, 6, 0)
                        end if
                    end if
                    !--          Check albedo parameters. If albedo_type is 0, all parameters must be set.
                    if (albedo_type_f % from_file) then
                        if (albedo_type_f % var(j, i) == 0) then
                            if (any(albedo_pars_f % pars_xy(:, j, i) == albedo_pars_f % fill)) then
                                write (*, *) 'albedo_pars missing at location (i,j) = (', &
                                    i, ',', j, ')'
                                ! call message('netcdf_data_input_mod', 'DRV0035', 2, 2, myid, 6, 0)
                            end if
                        end if
                    end if

                    !--          Check pavement parameters. If pavement_type is 0, all parameters of pavement_pars must
                    !--          be set at this location.
                    if (pavement_type_f % from_file) then
                        if (pavement_type_f % var(j, i) == 0) then
                            if (any(pavement_pars_f % pars_xy(:, j, i) == pavement_pars_f % fill)) then
                                write (*, *) 'pavement_pars missing at location (i,j) = (', &
                                    i, ',', j, ')'
                                ! call message('netcdf_data_input_mod', 'DRV0036', 2, 2, myid, 6, 0)
                            end if
                        end if
                    end if
                    !--          Check pavement-subsurface parameters. If pavement_type is 0, all parameters of
                    !--          pavement_subsurface_pars must be set at this location.
                    if (pavement_type_f % from_file) then
                        if (pavement_type_f % var(j, i) == 0) then
                            if (any(pavement_subsurface_pars_f % pars_xyz(:, :, j, i) == &
                                    pavement_subsurface_pars_f % fill)) &
                                then
                                write (*, *) 'pavement_subsurface_pars missing at location ', &
                                    '(i,j) = (', i, ',', j, ')'
                                ! call message('netcdf_data_input_mod', 'DRV0037', 2, 2, myid, 6, 0)
                            end if
                        end if
                    end if

                    !--          Check water parameters. If water_type is 0, all parameters must be set  at this
                    !--          location.
                    if (water_type_f % from_file) then
                        if (water_type_f % var(j, i) == 0) then
                            if (any(water_pars_f % pars_xy(:, j, i) == water_pars_f % fill)) then
                                write (*, *) 'water_pars missing at location (i,j) = (', &
                                    i, ',', j, ')'
                                ! call message('netcdf_data_input_mod', 'DRV0038', 2, 2, myid, 6, 0)
                            end if
                        end if
                    end if

                end do
            end do
        end if

    end subroutine netcdf_data_input_check_static

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Checks if a given variables is on file
    !--------------------------------------------------------------------------------------------------!
    function check_existence(vars_in_file, var_name)

        implicit none

        character(LEN=*) :: var_name !< variable to be checked
        character(LEN=*), dimension(:) :: vars_in_file !< list of variables in file

        integer(iwp) :: i !< loop variable

        logical :: check_existence !< flag indicating whether a variable exist or not - actual return value

        i = 1
        check_existence = .false.
        do while (i <= size(vars_in_file))
            check_existence = trim(vars_in_file(i)) == trim(var_name) .or. check_existence
            i = i + 1
        end do

        return

    end function check_existence

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Closes an existing netCDF file.
    !--------------------------------------------------------------------------------------------------!
    subroutine close_input_file(id)

        implicit none

        integer(iwp), intent(INOUT) :: id !< file id

        ! #if defined( __netcdf )
        nc_stat = NF90_CLOSE(id)
        call handle_error('close', 540)

    ! #endif
    end subroutine close_input_file

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Opens an existing netCDF file for reading only and returns its id.
    !--------------------------------------------------------------------------------------------------!
    subroutine open_read_file(filename, id)

        implicit none

        character(LEN=*), intent(IN) :: filename !< filename
        integer(iwp), intent(INOUT) :: id !< file id
        logical :: file_exists !< true if file exists

        ! #if defined( __netcdf )
        !-- Check if requested file exists
        inquire (FILE=filename, EXIST=file_exists)

        if (.not. file_exists) then
            write (*, *) 'required input file "' // filename // '" does not exist'
            stop
            ! ! call message('open_read_file', 'DRV0039', 2, 2, 0, 6, 1)
        end if

        ! #if defined( __netcdf4_parallel )
        !         !-- If __netcdf4_parallel is defined, parrallel NetCDF will be used.
        !         nc_stat = NF90_OPEN(filename, ior(NF90_NOWRITE, NF90_MPIIO), id, COMM=comm2d, &
        !                             INFO=MPI_INFO_NULL)
        !         !-- In case the previous open call fails, check for possible Netcdf 3 file, and open it. However,
        !         !-- this case, disable parallel access.
        !         if (nc_stat /= NF90_NOERR) then
        !             nc_stat = NF90_OPEN(filename, NF90_NOWRITE, id)
        !             collective_read = .false.
        !         else
        ! ! #if defined( __nec )
        !             ! collective_read = .false. ! collective read causes hang situations on NEC Aurora
        ! ! #else
        !             collective_read = .true.
        ! ! #endif
        !         end if
        ! #else
        !-- All MPI processes open the file and read it (but not in parallel).
        nc_stat = NF90_OPEN(filename, NF90_NOWRITE, id)

        ! #endif
        call handle_error('open_read_file', 539)

    ! #endif
    end subroutine open_read_file

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Reads global or variable-related attributes of type INTEGER (32-bit)
    !--------------------------------------------------------------------------------------------------!
    subroutine get_attribute_int32(id, attribute_name, value, global, variable_name, ignore_error)

        implicit none

        character(LEN=*) :: attribute_name !< attribute name
        character(LEN=*), optional :: variable_name !< variable name

        integer(iwp), intent(IN) :: id !< file id
        integer(iwp) :: id_var !< variable id
        integer(iwp), intent(INOUT) :: value !< read value

        logical :: check_error !< flag indicating if handle_error shall be checked
        logical, intent(IN) :: global !< flag indicating global attribute
        logical, intent(IN), optional :: ignore_error !< flag indicating if errors should be checked or not

        ! #if defined( __netcdf )
        if (present(ignore_error)) then
            check_error = .not. ignore_error
        else
            check_error = .true.
        end if
        !-- Read global attribute
        if (global) then
            nc_stat = NF90_GET_ATT(id, NF90_GLOBAL, trim(attribute_name), value)
            if (check_error) call handle_error('get_attribute_int32 global', 522, attribute_name)
        !-- Read attributes referring to a single variable. Therefore, first inquire variable id
        else
            nc_stat = NF90_INQ_VARID(id, trim(variable_name), id_var)
            if (check_error) call handle_error('get_attribute_int32', 522, attribute_name)
            nc_stat = NF90_GET_ATT(id, id_var, trim(attribute_name), value)
            if (check_error) call handle_error('get_attribute_int32', 522, attribute_name)
        end if

    ! #endif
    end subroutine get_attribute_int32

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Reads global or variable-related attributes of type INTEGER (8-bit)
    !--------------------------------------------------------------------------------------------------!
    subroutine get_attribute_int8(id, attribute_name, value, global, variable_name, ignore_error)

        implicit none

        character(LEN=*) :: attribute_name !< attribute name
        character(LEN=*), optional :: variable_name !< variable name

        integer(iwp), intent(IN) :: id !< file id
        integer(iwp) :: id_var !< variable id
        integer(ibp), intent(INOUT) :: value !< read value

        logical :: check_error !< flag indicating if handle_error shall be checked
        logical, intent(IN) :: global !< flag indicating global attribute
        logical, intent(IN), optional :: ignore_error !< flag indicating if errors should be checked or not

        ! #if defined( __netcdf )
        if (present(ignore_error)) then
            check_error = .not. ignore_error
        else
            check_error = .true.
        end if
        !-- Read global attribute
        if (global) then
            nc_stat = NF90_GET_ATT(id, NF90_GLOBAL, trim(attribute_name), value)
            if (check_error) call handle_error('get_attribute_int8 global', 523, attribute_name)
        !-- Read attributes referring to a single variable. Therefore, first inquire variable id
        else
            nc_stat = NF90_INQ_VARID(id, trim(variable_name), id_var)
            if (check_error) call handle_error('get_attribute_int8', 523, attribute_name)
            nc_stat = NF90_GET_ATT(id, id_var, trim(attribute_name), value)
            if (check_error) call handle_error('get_attribute_int8', 523, attribute_name)
        end if

    ! #endif
    end subroutine get_attribute_int8

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Reads global or variable-related attributes of type REAL
    !--------------------------------------------------------------------------------------------------!
    subroutine get_attribute_real(id, attribute_name, value, global, variable_name, ignore_error)

        implicit none

        character(LEN=*) :: attribute_name !< attribute name
        character(LEN=*), optional :: variable_name !< variable name

        integer(iwp), intent(IN) :: id !< file id
        integer(iwp) :: id_var !< variable id

        logical :: check_error !< flag indicating if handle_error shall be checked
        logical, intent(IN) :: global !< flag indicating global attribute
        logical, intent(IN), optional :: ignore_error !< flag indicating if errors should be checked or not

        real(wp), intent(INOUT) :: value !< read value

        ! #if defined( __netcdf )
        if (present(ignore_error)) then
            check_error = .not. ignore_error
        else
            check_error = .true.
        end if
        !-- Read global attribute
        if (global) then
            nc_stat = NF90_GET_ATT(id, NF90_GLOBAL, trim(attribute_name), value)
            if (check_error) call handle_error('get_attribute_real global', 524, attribute_name)
        !-- Read attributes referring to a single variable. Therefore, first inquire variable id
        else
            nc_stat = NF90_INQ_VARID(id, trim(variable_name), id_var)
            if (check_error) call handle_error('get_attribute_real', 524, attribute_name)
            nc_stat = NF90_GET_ATT(id, id_var, trim(attribute_name), value)
            if (check_error) call handle_error('get_attribute_real', 524, attribute_name)
        end if

    ! #endif
    end subroutine get_attribute_real

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Reads global or variable-related attributes of type CHARACTER
    !> Remark: reading attributes of type NF_STRING return an error code 56 -
    !> Attempt to convert between text & numbers.
    !--------------------------------------------------------------------------------------------------!
    subroutine get_attribute_string(id, attribute_name, value, global, variable_name, ignore_error)

        implicit none

        character(LEN=*) :: attribute_name !< attribute name
        character(LEN=*), intent(INOUT) :: value !< read value
        character(LEN=*), optional :: variable_name !< variable name

        integer(iwp), intent(IN) :: id !< file id
        integer(iwp) :: id_var !< variable id

        logical :: check_error !< flag indicating if handle_error shall be checked
        logical, intent(IN) :: global !< flag indicating global attribute
        logical, intent(IN), optional :: ignore_error !< flag indicating if errors should be checked or not

        ! #if defined( __netcdf )
        if (present(ignore_error)) then
            check_error = .not. ignore_error
        else
            check_error = .true.
        end if
        !-- Read global attribute
        if (global) then
            nc_stat = NF90_GET_ATT(id, NF90_GLOBAL, trim(attribute_name), value)
            if (check_error) call handle_error('get_attribute_string global', 525, attribute_name)
        !-- Read attributes referring to a single variable. Therefore, first inquire variable id
        else
            nc_stat = NF90_INQ_VARID(id, trim(variable_name), id_var)
            if (check_error) call handle_error('get_attribute_string', 525, attribute_name)

            nc_stat = NF90_GET_ATT(id, id_var, trim(attribute_name), value)
            if (check_error) call handle_error('get_attribute_string', 525, attribute_name)

        end if

    ! #endif
    end subroutine get_attribute_string

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Get dimension array for a given dimension
    !> Used.
    !--------------------------------------------------------------------------------------------------!
    subroutine get_dimension_length(id, dim_len, variable_name)

        implicit none

        character(LEN=100) :: dum !< dummy variable to receive return character
        character(LEN=*) :: variable_name !< dimension name

        integer(iwp) :: dim_len !< dimension size
        integer(iwp), intent(IN) :: id !< file id
        integer(iwp) :: id_dim !< dimension id

        ! #if defined( __netcdf )
        !-- First, inquire dimension ID
        nc_stat = NF90_INQ_DIMID(id, trim(variable_name), id_dim)
        call handle_error('get_dimension_length', 526, variable_name)
        !-- Inquire dimension length
        nc_stat = NF90_INQUIRE_DIMENSION(id, id_dim, dum, LEN=dim_len)
        call handle_error('get_dimension_length', 526, variable_name)


    ! #endif
    end subroutine get_dimension_length

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Routine for reading-in a character string from the chem emissions netcdf input file.
    !--------------------------------------------------------------------------------------------------!
    subroutine get_variable_string(id, variable_name, var_string, names_number)

        implicit none

        character(LEN=*) :: variable_name !< variable name
        character(LEN=25), allocatable, dimension(:), intent(INOUT) :: var_string
        character(LEN=1), allocatable, dimension(:, :) :: tmp_var_string !< variable to be read

        integer(iwp) :: i, j !< index to go through the length of the dimensions
        integer(iwp), intent(IN) :: id !< file id
        integer(iwp) :: id_var !< variable id
        integer(iwp), intent(IN) :: names_number !< number of names
        integer(iwp) :: max_string_length = 25 !< this is both the maximum length of a name, but also the number of the
        !< components of the first dimensions (rows)

        ! #if defined( __netcdf )
        allocate (tmp_var_string(max_string_length, names_number))

        allocate (var_string(names_number))

        !-- Inquire variable id
        nc_stat = NF90_INQ_VARID(id, trim(variable_name), id_var)

        !-- Get variable
        !-- Start cycle over the emission species
        do i = 1, names_number
            !--    Read the first letter of each component
            nc_stat = NF90_GET_VAR(id, id_var, var_string(i), start=(/1, i/), count=(/1, 1/))
            call handle_error('get_variable_string', 701)
            !--    Start cycle over charachters
            do j = 1, max_string_length
                !--       Read the rest of the components of the name
                nc_stat = NF90_GET_VAR(id, id_var, tmp_var_string(j, i), start=(/j, i/), &
                                       count=(/1, 1/))
                call handle_error('get_variable_string', 702)

                if (iachar(tmp_var_string(j, i)) == 0) then
                    tmp_var_string(j, i) = ''
                end if

                if (j > 1) then
                    !--          Concatenate first letter of the name and the others
                    var_string(i) = trim(var_string(i)) // trim(tmp_var_string(j, i))
                end if
            end do
        end do


    ! #endif
    end subroutine get_variable_string

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Generalized routine for reading strings from a netCDF variable to replace existing
    !> get_variable_string ( )
    !>
    !> Improvements:
    !>   - Expanded string size limit from 25 to 512
    !>   - No longer removes spaces between text magically (this seems to have been aimed at a very
    !>     specific application, but I don't know what)
    !>   - Simplified implementation
    !>
    !> Caveats:
    !>   - Somehow I could not get the subroutine to work with str_array(:,:) so I reverted to a
    !>     hard-coded str_array(:,512), hopefully large enough for most general applications. This also
    !>     means the character variable used for str_array must be of size (:,512)
    !--------------------------------------------------------------------------------------------------!
    subroutine get_variable_string_generic(id, var_name, str_array, num_str, str_len)

        implicit none

        character(LEN=*), intent(IN) :: var_name !< netCDF variable name
        character(LEN=512), allocatable, intent(INOUT) :: str_array(:) !< output string array

        integer(iwp) :: buf_len !< string buffer size
        integer(iwp) :: id_var !< netCDF variable ID
        integer(iwp), intent(IN) :: id !< netCDF file ID
        integer(iwp) :: k !< generic counter
        integer(iwp), intent(IN) :: num_str !< number of string elements in array
        integer(iwp), intent(IN) :: str_len !< size of each string element


        ! #if defined( __netcdf )
        !-- Set buffer length to up to hard-coded string size
        buf_len = min(abs(str_len), 512)

        !-- Allocate necessary memories for string variables
        allocate (str_array(num_str))
        !-- Get variable id
        nc_stat = NF90_INQ_VARID(id, trim(var_name), id_var)
        !-- Extract string variables
        do k = 1, num_str
            str_array(k) = ''
            nc_stat = NF90_GET_VAR(id, id_var, str_array(k), start=(/1, k/), &
                                   count=(/buf_len, 1/))
            call handle_error('get_variable_string_generic', 702)
        end do


    ! #endif
    end subroutine get_variable_string_generic

    !--------------------------------------------------------------------------------------------------!
    !    Description:
    !    ------------
    !>   Reads a character variable in a 1D array
    !--------------------------------------------------------------------------------------------------!
    subroutine get_variable_1d_char(id, variable_name, var)

        implicit none

        character(LEN=*) :: variable_name !< variable name

        character(LEN=*), dimension(:), intent(INOUT) :: var !< variable to be read

        integer(iwp) :: i !< running index over variable dimension
        integer(iwp), intent(IN) :: id !< file id
        integer(iwp) :: id_var !< dimension id

        integer(iwp), dimension(2) :: dimid !< dimension IDs
        integer(iwp), dimension(2) :: dimsize !< dimension size


        ! #if defined( __netcdf )
        !-- First, inquire variable ID
        nc_stat = NF90_INQ_VARID(id, trim(variable_name), id_var)
        call handle_error('get_variable_1d_int', 527, variable_name)
        !-- Inquire dimension IDs
        nc_stat = NF90_INQUIRE_VARIABLE(id, id_var, dimids=dimid(1:2))
        call handle_error('get_variable_1d_char', 527, variable_name)
        !-- Read dimesnion length
        nc_stat = NF90_INQUIRE_DIMENSION(id, dimid(1), LEN=dimsize(1))
        nc_stat = NF90_INQUIRE_DIMENSION(id, dimid(2), LEN=dimsize(2))

        !-- Read character array. Note, each element is read individually, in order to better separate
        !-- single strings.
        do i = 1, dimsize(2)
            nc_stat = NF90_GET_VAR(id, id_var, var(i), &
                                   start=(/1, i/), &
                                   count=(/dimsize(1), 1/))
            call handle_error('get_variable_1d_char', 527, variable_name)
        end do


    ! #endif
    end subroutine get_variable_1d_char

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Reads a 1D integer variable from file.
    !--------------------------------------------------------------------------------------------------!
    subroutine get_variable_1d_int(id, variable_name, var)

        implicit none

        character(LEN=*) :: variable_name !< variable name

        integer(iwp), intent(IN) :: id !< file id
        integer(iwp) :: id_var !< dimension id

        integer(iwp), dimension(:), intent(INOUT) :: var !< variable to be read

        ! #if defined( __netcdf )
        !-- First, inquire variable ID
        nc_stat = NF90_INQ_VARID(id, trim(variable_name), id_var)
        call handle_error('get_variable_1d_int', 527, variable_name)
        !-- Read variable
        nc_stat = NF90_GET_VAR(id, id_var, var)
        call handle_error('get_variable_1d_int', 527, variable_name)


    ! #endif
    end subroutine get_variable_1d_int

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Reads a 1D float variable from file.
    !--------------------------------------------------------------------------------------------------!
    subroutine get_variable_1d_real(id, variable_name, var, is, count_elements)

        implicit none

        character(LEN=*) :: variable_name !< variable name

        integer(iwp), intent(IN), optional :: count_elements !< number of elements to be read
        integer(iwp), intent(IN) :: id !< file id
        integer(iwp) :: id_var !< dimension id
        integer(iwp), intent(IN), optional :: is !< start index

        real(wp), dimension(:), intent(INOUT) :: var !< variable to be read
        ! #if defined( __netcdf )
        !-- First, inquire variable ID
        nc_stat = NF90_INQ_VARID(id, trim(variable_name), id_var)
        call handle_error('get_variable_1d_real', 528, variable_name)
        !-- Read variable
        if (present(is)) then
            nc_stat = NF90_GET_VAR(id, id_var, var, start=(/is/), count=(/count_elements/))
            call handle_error('get_variable_1d_real', 528, variable_name)
        else
            nc_stat = NF90_GET_VAR(id, id_var, var)
            call handle_error('get_variable_1d_real', 528, variable_name)
        end if


    ! #endif
    end subroutine get_variable_1d_real

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Reads a time-dependent 1D float variable from file.
    !--------------------------------------------------------------------------------------------------!
    subroutine get_variable_pr(id, variable_name, t, var)

        implicit none

        character(LEN=*) :: variable_name !< variable name

        integer(iwp), intent(IN) :: id !< file id
        integer(iwp), dimension(1:2) :: id_dim !< dimension ids
        integer(iwp) :: id_var !< dimension id
        integer(iwp) :: n_file !< number of data-points in file along z dimension
        integer(iwp), intent(IN) :: t !< timestep number

        real(wp), dimension(:), intent(INOUT) :: var !< variable to be read

        ! #if defined( __netcdf )
        !-- First, inquire variable ID
        nc_stat = NF90_INQ_VARID(id, trim(variable_name), id_var)
        !-- Inquire dimension size of vertical dimension
        nc_stat = NF90_INQUIRE_VARIABLE(id, id_var, DIMIDS=id_dim)
        nc_stat = NF90_INQUIRE_DIMENSION(id, id_dim(1), LEN=n_file)
        !-- Read variable.
        nc_stat = NF90_GET_VAR(id, id_var, var, &
                               start=(/1, t/), &
                               count=(/n_file, 1/))
        call handle_error('get_variable_pr', 529, variable_name)

    ! #endif
    end subroutine get_variable_pr

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Reads a per-surface pars variable from file. Because all surfaces are stored as flat 1-D array,
    !> each PE has to scan the data and find the surface indices belonging to its subdomain. During this
    !> scan, it also builds a necessary (j,i) index.
    !--------------------------------------------------------------------------------------------------!
    subroutine get_variable_surf(id, variable_name, surf)

        use grid_variables, &
            only: ddx, &
                  ddy

        use basic_constants_and_equations_mod, &
            only: pi

        implicit none

        integer(iwp), parameter :: nsurf_pars_read = 2**15 !< read buffer size (value > 10^15 makes problem with
        !< ifort)

        character(LEN=*) :: variable_name !< variable name

        integer(iwp) :: i !< grid index in x-direction
        integer(iwp) :: j !< grid index in y-direction
        integer(iwp), intent(IN) :: id !< file id
        integer(iwp) :: id_azimuth !< azimuth variable id
        integer(iwp) :: id_var !< variable id
        integer(iwp) :: id_xs !< xs variable id
        integer(iwp) :: id_ys !< ys variable id
        integer(iwp) :: id_zs !< zs variable id
        integer(iwp) :: id_zenith !< zeith variable id
        integer(iwp) :: is !< local surface index
        integer(iwp) :: is0, isc !< read surface start and count
        integer(iwp) :: isurf !< netcdf surface index
        integer(iwp) :: nsurf !< total number of surfaces in file

        integer(iwp), dimension(6) :: coords !< integer coordinates of surface location
        integer(iwp), dimension(2) :: id_dim !< dimension ids

        integer(iwp), dimension(:, :), allocatable :: nsurf_ji !< numbers of surfaces by coords

        real(wp) :: oro_max_l !< maximum terrain height under building

        real(wp), dimension(:), allocatable :: azimuth !< read buffer for azimuth(s)
        real(wp), dimension(:), allocatable :: zenith !< read buffer for zenith(s)
        real(wp), dimension(:), allocatable :: xs !< surface coordinate array of x-dimension
        real(wp), dimension(:), allocatable :: ys !< surface coordinate array of y-dimension
        real(wp), dimension(:), allocatable :: zs !< surface coordinate array of z-dimension

        real(wp), dimension(:, :), allocatable :: pars_read !< read buffer for the building parameters

        type(pars_surf) :: surf !< parameters variable to be loaded

        ! #if defined( __netcdf )
        !-- First, inquire variable ID's
        nc_stat = NF90_INQ_VARID(id, trim(variable_name), id_var)
        nc_stat = NF90_INQ_VARID(id, 'zs', id_zs)
        nc_stat = NF90_INQ_VARID(id, 'ys', id_ys)
        nc_stat = NF90_INQ_VARID(id, 'xs', id_xs)
        nc_stat = NF90_INQ_VARID(id, 'zenith', id_zenith)
        nc_stat = NF90_INQ_VARID(id, 'azimuth', id_azimuth)
        !-- Inquire dimension sizes for the number of surfaces and parameters given
        nc_stat = NF90_INQUIRE_VARIABLE(id, id_var, DIMIDS=id_dim)
        nc_stat = NF90_INQUIRE_DIMENSION(id, id_dim(1), LEN=nsurf)
        nc_stat = NF90_INQUIRE_DIMENSION(id, id_dim(2), LEN=surf % np)

        allocate (pars_read(nsurf_pars_read, surf % np), &
                  zs(nsurf_pars_read), ys(nsurf_pars_read), &
                  xs(nsurf_pars_read), zenith(nsurf_pars_read), &
                  azimuth(nsurf_pars_read), &
                  nsurf_ji(nys:nyn, nxl:nxr))

        nsurf_ji(:, :) = 0
        !-- Scan surface coordinates, count locally
        is0 = 1
        do
            isc = min(nsurf_pars_read, nsurf - is0 + 1)
            if (isc <= 0) exit

            nc_stat = NF90_GET_VAR(id, id_ys, ys, &
                                   start=(/is0/), &
                                   count=(/isc/))
            nc_stat = NF90_GET_VAR(id, id_xs, xs, &
                                   start=(/is0/), &
                                   count=(/isc/))
            nc_stat = NF90_GET_VAR(id, id_zenith, zenith, &
                                   start=(/is0/), &
                                   count=(/isc/))
            nc_stat = NF90_GET_VAR(id, id_azimuth, azimuth, &
                                   start=(/is0/), &
                                   count=(/isc/))
            call handle_error('get_variable_surf', 682, 'azimuth')

            do isurf = 1, isc
                !--       Parse coordinates, detect if belongs to subdomain
                coords = transform_coords(xs(isurf), ys(isurf), zenith(isurf), azimuth(isurf))
                if (coords(2) < nys .or. coords(2) > nyn .or. &
                    coords(3) < nxl .or. coords(3) > nxr) cycle

                nsurf_ji(coords(2), coords(3)) = nsurf_ji(coords(2), coords(3)) + 1
            end do
            is0 = is0 + isc
        end do
        !-- Populate reverse index from surface counts
        allocate (surf % index_ji(2, nys:nyn, nxl:nxr))
        isurf = 1
        do j = nys, nyn
            do i = nxl, nxr
                surf % index_ji(:, j, i) = (/isurf, isurf + nsurf_ji(j, i) - 1/)
                isurf = isurf + nsurf_ji(j, i)
            end do
        end do

        surf % nsurf = isurf - 1
        allocate (surf % pars(0:surf % np - 1, surf % nsurf), &
                  surf % coords(6, surf % nsurf))
        !-- Scan surfaces again, saving pars into allocated structures
        nsurf_ji(:, :) = 0
        is0 = 1
        do
            isc = min(nsurf_pars_read, nsurf - is0 + 1)
            if (isc <= 0) exit

            nc_stat = NF90_GET_VAR(id, id_var, pars_read(1:isc, 1:surf % np), &
                                   start=(/is0, 1/), &
                                   count=(/isc, surf % np/))
            call handle_error('get_variable_surf', 683, variable_name)

            nc_stat = NF90_GET_VAR(id, id_zs, zs, &
                                   start=(/is0/), &
                                   count=(/isc/))
            nc_stat = NF90_GET_VAR(id, id_ys, ys, &
                                   start=(/is0/), &
                                   count=(/isc/))
            nc_stat = NF90_GET_VAR(id, id_xs, xs, &
                                   start=(/is0/), &
                                   count=(/isc/))
            nc_stat = NF90_GET_VAR(id, id_zenith, zenith, &
                                   start=(/is0/), &
                                   count=(/isc/))
            nc_stat = NF90_GET_VAR(id, id_azimuth, azimuth, &
                                   start=(/is0/), &
                                   count=(/isc/))

            do isurf = 1, isc
                !--       Parse coordinates, detect if belongs to subdomain
                coords = transform_coords(xs(isurf), ys(isurf), zenith(isurf), azimuth(isurf))
                if (coords(2) < nys .or. coords(2) > nyn .or. &
                    coords(3) < nxl .or. coords(3) > nxr) cycle
                !--       Determine maximum terrain under building (base z-coordinate). Using normal vector to
                !--       locate building inner coordinates.
                oro_max_l = buildings_f % oro_max(coords(2) - coords(5), coords(3) - coords(6))
                if (oro_max_l == buildings_f % fill1) then
                    write (*, *) 'found building surface on ' // &
                        'non-building coordinates (xs, ys, zenith, azimuth): ', &
                        xs(isurf), ys(isurf), zenith(isurf), azimuth(isurf)
                    ! call message('get_variable_surf', 'DRV0040', 2, 2, myid, 6, 0)
                end if
                !--       Urban layer has no stretching, therefore using dz(1) instead of linear searching through
                !--       zu/zw.
                coords(1) = nint((zs(isurf) + oro_max_l) / dz(1) + 0.5_wp + 0.5_wp * coords(4), &
                                 KIND=iwp)
                !--       Save surface entry
                is = surf % index_ji(1, coords(2), coords(3)) + nsurf_ji(coords(2), coords(3))
                surf % pars(:, is) = pars_read(isurf, :)
                surf % coords(:, is) = coords(:)

                nsurf_ji(coords(2), coords(3)) = nsurf_ji(coords(2), coords(3)) + 1
            end do

            is0 = is0 + isc
        end do

        deallocate (pars_read, zs, ys, xs, zenith, azimuth, nsurf_ji)

    contains

        pure function transform_coords(x, y, zenith, azimuth)

            integer(iwp), dimension(6) :: transform_coords !< (k,j,i,norm_z,norm_y,norm_x)

            real(wp), intent(in) :: azimuth !< surface normal azimuth angle in degrees
            real(wp), intent(in) :: x !< surface centre coordinate in x in metres from origin
            real(wp), intent(in) :: y !< surface centre coordinate in y in metres from origin
            real(wp), intent(in) :: zenith !< surface normal zenith angle in degrees

            transform_coords(4) = nint(cos(zenith * pi / 180.0_wp), KIND=iwp)
            if (transform_coords(4) == 0) then
                transform_coords(5) = nint(cos(azimuth * pi / 180.0_wp), KIND=iwp)
                transform_coords(6) = nint(sin(azimuth * pi / 180.0_wp), KIND=iwp)
            else
                transform_coords(5) = 0.0_wp
                transform_coords(6) = 0.0_wp
            end if

            transform_coords(1) = -999.0_wp ! not calculated here
            transform_coords(2) = nint(y * ddy - 0.5_wp + 0.5_wp * transform_coords(5), KIND=iwp)
            transform_coords(3) = nint(x * ddx - 0.5_wp + 0.5_wp * transform_coords(6), KIND=iwp)

        end function transform_coords

    ! #endif
    end subroutine get_variable_surf

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Reads a 2D REAL variable from a file. Reading is done processor-wise, i.e. each core reads its
    !> own domain in slices along x.
    !--------------------------------------------------------------------------------------------------!
    subroutine get_variable_2d_real(id, variable_name, var, is, ie, js, je, nbgp)

        implicit none

        character(LEN=*) :: variable_name !< variable name

        integer(iwp), intent(IN) :: id !< file id
        integer(iwp), intent(IN) :: ie !< start index for subdomain input along x direction
        integer(iwp), intent(IN) :: is !< end index for subdomain input along x direction
        integer(iwp), intent(IN) :: je !< start index for subdomain input along y direction
        integer(iwp), intent(IN) :: js !< end index for subdomain input along y direction
        integer(iwp), intent(IN) :: nbgp !< number of ghost layers of variable to be read

        integer(iwp) :: i !< running index along x direction
        integer(iwp) :: id_var !< variable id
        integer(iwp) :: j !< running index along y direction

        real(wp), dimension(:, :), allocatable :: tmp !< temporary variable to read data from file according
        !< to its reverse memory access
        real(wp), dimension(js - nbgp:je + nbgp, is - nbgp:ie + nbgp), intent(INOUT) :: var !< variable to be read

        ! #if defined( __netcdf )
        !-- Inquire variable id
        nc_stat = NF90_INQ_VARID(id, trim(variable_name), id_var)
        !-- Check for collective read-operation and set respective NetCDF flags if required.
        if (collective_read) then
        ! #if defined( __netcdf4_parallel )
        !             nc_stat = NF90_VAR_PAR_ACCESS(id, id_var, NF90_COLLECTIVE)
        ! #endif
        end if

        !-- Allocate temporary variable according to memory access on file.
        allocate (tmp(is:ie, js:je))
        !-- Get variable
        nc_stat = NF90_GET_VAR(id, id_var, tmp, start=(/is + 1, js + 1/), &
                               count=(/ie - is + 1, je - js + 1/))
        call handle_error('get_variable_2d_real', 530, variable_name)
        !-- Resort data. Please note, dimension subscripts of var all start at 1.
        !    IF ( myid == 0 )  PRINT*, '### get_variable_2d_real: "', TRIM( variable_name )
        do i = is, ie
            do j = js, je
                !          var(j-js+1,i-is+1) = tmp(i,j)
                var(j, i) = tmp(i, j)
            end do
        end do

        deallocate (tmp)

    ! #endif
    end subroutine get_variable_2d_real

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Reads a 2D 32-bit INTEGER variable from file. Reading is done processor-wise, i.e. each core
    !> reads its own domain in slices along x.
    !--------------------------------------------------------------------------------------------------!
    subroutine get_variable_2d_int32(id, variable_name, var, is, ie, js, je, nbgp)

        implicit none

        character(LEN=*) :: variable_name !< variable name

        integer(iwp), intent(IN) :: id !< file id
        integer(iwp), intent(IN) :: ie !< start index for subdomain input along x direction
        integer(iwp), intent(IN) :: is !< end index for subdomain input along x direction
        integer(iwp), intent(IN) :: je !< start index for subdomain input along y direction
        integer(iwp), intent(IN) :: js !< end index for subdomain input along y direction
        integer(iwp), intent(IN) :: nbgp !< number of ghost layers of variable to be read

        integer(iwp) :: i !< running index along x direction
        integer(iwp) :: id_var !< variable id
        integer(iwp) :: j !< running index along y direction

        integer(iwp), dimension(:, :), allocatable :: tmp !< temporary variable to read data from file according
        !< to its reverse memory access
        !    INTEGER(iwp), DIMENSION(:,:), INTENT(INOUT) ::  var  !< variable to be read
        integer(iwp), dimension(js - nbgp:je + nbgp, is - nbgp:ie + nbgp), intent(INOUT) :: var !< variable to be read

        ! #if defined( __netcdf )
        !-- Inquire variable id
        nc_stat = NF90_INQ_VARID(id, trim(variable_name), id_var)
        !-- Check for collective read-operation and set respective NetCDF flags if required.
        if (collective_read) then
        ! #if defined( __netcdf4_parallel )
        !             nc_stat = NF90_VAR_PAR_ACCESS(id, id_var, NF90_COLLECTIVE)
        ! #endif
        end if
        !-- Allocate temporary variable according to memory access on file.
        allocate (tmp(is:ie, js:je))
        !-- Get variable
        nc_stat = NF90_GET_VAR(id, id_var, tmp, start=(/is + 1, js + 1/), &
                               count=(/ie - is + 1, je - js + 1/))

        call handle_error('get_variable_2d_int32', 531, variable_name)
        !-- Resort data. Please note, dimension subscripts of var all start at 1.
        !    IF ( myid == 0 )  PRINT*, '### get_variable_2d_int32: "', TRIM( variable_name )
        do i = is, ie
            do j = js, je
                !          var(j-js+1,i-is+1) = tmp(i,j)
                var(j, i) = tmp(i, j)
            end do
        end do

        deallocate (tmp)

    ! #endif
    end subroutine get_variable_2d_int32

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Reads a 2D 8-bit INTEGER variable from file. Reading is done processor-wise, i.e. each core reads
    !> its own domain in slices along x.
    !--------------------------------------------------------------------------------------------------!
    subroutine get_variable_2d_int8(id, variable_name, var, is, ie, js, je, nbgp)

        implicit none

        character(LEN=*) :: variable_name !< variable name

        integer(iwp), intent(IN) :: id !< file id
        integer(iwp), intent(IN) :: ie !< start index for subdomain input along x direction
        integer(iwp), intent(IN) :: je !< start index for subdomain input along y direction
        integer(iwp), intent(IN) :: js !< end index for subdomain input along y direction
        integer(iwp), intent(IN) :: is !< end index for subdomain input along x direction
        integer(iwp), intent(IN) :: nbgp !< number of ghost layers of variable to be read

        integer(iwp) :: i !< running index along x direction
        integer(iwp) :: id_var !< variable id
        integer(iwp) :: j !< running index along y direction

        integer(ibp), dimension(:, :), allocatable :: tmp !< temporary variable to read data from file according
        !< to its reverse memory access
        !    INTEGER(ibp), DIMENSION(:,:), INTENT(INOUT) ::  var  !< variable to be read
        integer(ibp), dimension(js - nbgp:je + nbgp, is - nbgp:ie + nbgp), intent(INOUT) :: var !< variable to be read

        ! #if defined( __netcdf )
        !-- Inquire variable id
        nc_stat = NF90_INQ_VARID(id, trim(variable_name), id_var)
        !-- Check for collective read-operation and set respective NetCDF flags if required.
        if (collective_read) then
        ! #if defined( __netcdf4_parallel )
        !             nc_stat = NF90_VAR_PAR_ACCESS(id, id_var, NF90_COLLECTIVE)
        ! #endif
        end if
        !-- Allocate temporary variable according to memory access on file.
        allocate (tmp(is:ie, js:je))
        !-- Get variable
        nc_stat = NF90_GET_VAR(id, id_var, tmp, &
                               start=(/is + 1, js + 1/), &
                               count=(/ie - is + 1, je - js + 1/))

        call handle_error('get_variable_2d_int8', 532, variable_name)
        !-- Resort data. Please note, dimension subscripts of var all start at 1.
        !    IF ( myid == 0 )  PRINT*, '### get_variable_2d_int8: "', TRIM( variable_name )
        do i = is, ie
            do j = js, je
                !          var(j-js+1,i-is+1) = tmp(i,j)
                var(j, i) = tmp(i, j)
            end do
        end do

        deallocate (tmp)

    ! #endif
    end subroutine get_variable_2d_int8

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Reads a 2D REAL variable from a file. This case, the variable is no spatial data but an
    !> arbitrary 2D list. Only dimension size is passed.
    !--------------------------------------------------------------------------------------------------!
    subroutine get_variable_2d_list_int(id, variable_name, var, n1, n2)

        use indices
        use pegrid

        implicit none

        character(LEN=*) :: variable_name !< variable name

        integer(iwp) :: i !< running index along x direction
        integer(iwp), intent(IN) :: id !< file id
        integer(iwp) :: id_var !< variable id
        integer(iwp) :: j !< running index along y direction
        integer(iwp) :: n1 !< length of first dimension
        integer(iwp) :: n2 !< length of second dimension

        integer(iwp), dimension(:, :), allocatable :: tmp !< temporary variable to read data from file according
        !< to its reverse memory access
        integer(iwp), dimension(:, :), intent(INOUT) :: var !< variable to be read

        ! #if defined( __netcdf )
        !-- Inquire variable id
        nc_stat = NF90_INQ_VARID(id, trim(variable_name), id_var)

        !-- Allocate temporary variable according to memory access on file.
        allocate (tmp(1:n2, 1:n1))
        !-- Get variable
        nc_stat = NF90_GET_VAR(id, id_var, tmp, start=(/1, 1/), count=(/n2, n1/))
        call handle_error('get_variable_2d_list_int', 530, variable_name)
        !-- Resort data. Please note, dimension subscripts of var all start at 1.
        do i = 1, n2
            do j = 1, n1
                var(j, i) = tmp(i, j)
            end do
        end do

        deallocate (tmp)

    ! #endif
    end subroutine get_variable_2d_list_int

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Reads a 2D REAL variable from a file. This case, the variable is no spatial data but an
    !> arbitrary 2D list. Only dimension size is passed.
    !--------------------------------------------------------------------------------------------------!
    subroutine get_variable_2d_list_real(id, variable_name, var, n1, n2)

        use indices
        use pegrid

        implicit none

        character(LEN=*) :: variable_name !< variable name

        integer(iwp) :: i !< running index along x direction
        integer(iwp), intent(IN) :: id !< file id
        integer(iwp) :: id_var !< variable id
        integer(iwp) :: j !< running index along y direction
        integer(iwp) :: n1 !< length of first dimension
        integer(iwp) :: n2 !< length of second dimension

        real(wp), dimension(:, :), allocatable :: tmp !< temporary variable to read data from file according
        !< to its reverse memory access
        real(wp), dimension(:, :), intent(INOUT) :: var !< variable to be read

        ! #if defined( __netcdf )
        !-- Inquire variable id
        nc_stat = NF90_INQ_VARID(id, trim(variable_name), id_var)

        !-- Allocate temporary variable according to memory access on file.
        allocate (tmp(1:n2, 1:n1))
        !-- Get variable
        nc_stat = NF90_GET_VAR(id, id_var, tmp, start=(/1, 1/), count=(/n2, n1/))
        call handle_error('get_variable_2d_list_real', 530, variable_name)
        !-- Resort data. Please note, dimension subscripts of var all start at 1.
        do i = 1, n2
            do j = 1, n1
                var(j, i) = tmp(i, j)
            end do
        end do

        deallocate (tmp)

    ! #endif
    end subroutine get_variable_2d_list_real

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Reads a 3D 8-bit INTEGER variable from file.
    !--------------------------------------------------------------------------------------------------!
    subroutine get_variable_3d_int8(id, variable_name, var, is, ie, js, je, ks, ke, nbgp)

        implicit none

        character(LEN=*) :: variable_name !< variable name

        integer(iwp), intent(IN) :: id !< file id
        integer(iwp), intent(IN) :: ie !< start index for subdomain input along x direction
        integer(iwp), intent(IN) :: is !< end index for subdomain input along x direction
        integer(iwp), intent(IN) :: je !< start index for subdomain input along y direction
        integer(iwp), intent(IN) :: js !< end index for subdomain input along y direction
        integer(iwp), intent(IN) :: ke !< start index of 3rd dimension
        integer(iwp), intent(IN) :: ks !< end index of 3rd dimension
        integer(iwp), intent(IN) :: nbgp !< number of ghost layers of variable to be read

        integer(iwp) :: i !< index along x direction
        integer(iwp) :: id_var !< variable id
        integer(iwp) :: j !< index along y direction
        integer(iwp) :: k !< index along any 3rd dimension

        integer(ibp), dimension(:, :, :), allocatable :: tmp !< temporary variable to read data from file according
        !< to its reverse memory access

        integer(ibp), dimension(ks:ke, js - nbgp:je + nbgp, is - nbgp:ie + nbgp), intent(INOUT) :: var !< variable to be read

        ! #if defined( __netcdf )
        !-- Inquire variable id
        nc_stat = NF90_INQ_VARID(id, trim(variable_name), id_var)
        !-- Check for collective read-operation and set respective NetCDF flags if required.
        if (collective_read) then
        ! #if defined( __netcdf4_parallel )
        !             nc_stat = NF90_VAR_PAR_ACCESS(id, id_var, NF90_COLLECTIVE)
        ! #endif
        end if
        !-- Allocate temporary variable according to memory access on file.
        allocate (tmp(is:ie, js:je, ks:ke))
        !-- Get variable
        nc_stat = NF90_GET_VAR(id, id_var, tmp, start=(/is + 1, js + 1, ks + 1/), &
                               count=(/ie - is + 1, je - js + 1, ke - ks + 1/))

        call handle_error('get_variable_3d_int8', 533, variable_name)
        !-- Resort data. Please note, dimension subscripts of var all start at 1.
        do i = is, ie
            do j = js, je
                do k = ks, ke
                    var(k, j, i) = tmp(i, j, k)
                end do
            end do
        end do

        deallocate (tmp)

    ! #endif
    end subroutine get_variable_3d_int8

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Reads a 3D float variable from file.
    !--------------------------------------------------------------------------------------------------!
    subroutine get_variable_3d_real(id, variable_name, var, is, ie, js, je, ks, ke, nbgp)

        implicit none

        character(LEN=*) :: variable_name !< variable name

        integer(iwp), intent(IN) :: id !< file id
        integer(iwp), intent(IN) :: ie !< start index for subdomain input along x direction
        integer(iwp), intent(IN) :: is !< end index for subdomain input along x direction
        integer(iwp), intent(IN) :: je !< start index for subdomain input along y direction
        integer(iwp), intent(IN) :: js !< end index for subdomain input along y direction
        integer(iwp), intent(IN) :: ke !< start index of 3rd dimension
        integer(iwp), intent(IN) :: ks !< end index of 3rd dimension
        integer(iwp), intent(IN) :: nbgp !< number of ghost layers of variable to be read

        integer(iwp) :: i !< index along x direction
        integer(iwp) :: id_var !< variable id
        integer(iwp) :: j !< index along y direction
        integer(iwp) :: k !< index along any 3rd dimension

        real(wp), dimension(:, :, :), allocatable :: tmp !< temporary variable to read data from file according
        !< to its reverse memory access

        real(wp), dimension(ks:ke, js - nbgp:je + nbgp, is - nbgp:ie + nbgp), intent(INOUT) :: var !< variable to be read

        ! #if defined( __netcdf )
        !-- Inquire variable id
        nc_stat = NF90_INQ_VARID(id, trim(variable_name), id_var)
        !-- Check for collective read-operation and set respective NetCDF flags if required.
        if (collective_read) then
        ! #if defined( __netcdf4_parallel )
        !             nc_stat = NF90_VAR_PAR_ACCESS(id, id_var, NF90_COLLECTIVE)
        ! #endif
        end if
        !-- Allocate temporary variable according to memory access on file.
        allocate (tmp(is:ie, js:je, ks:ke))
        !-- Get variable
        nc_stat = NF90_GET_VAR(id, id_var, tmp, start=(/is + 1, js + 1, ks + 1/), &
                               count=(/ie - is + 1, je - js + 1, ke - ks + 1/))

        call handle_error('get_variable_3d_real', 534, variable_name)
        !-- Resort data.
        do i = is, ie
            do j = js, je
                do k = ks, ke
                    var(k, j, i) = tmp(i, j, k)
                end do
            end do
        end do

        deallocate (tmp)

    ! #endif
    end subroutine get_variable_3d_real

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Reads a 4D float variable from file.
    !--------------------------------------------------------------------------------------------------!
    subroutine get_variable_4d_real(id, variable_name, var, is, ie, js, je, k1s, k1e, k2s, k2e, nbgp)

        implicit none

        character(LEN=*) :: variable_name !< variable name

        integer(iwp), intent(IN) :: id !< file id
        integer(iwp), intent(IN) :: ie !< start index for subdomain input along x direction
        integer(iwp), intent(IN) :: is !< end index for subdomain input along x direction
        integer(iwp), intent(IN) :: je !< start index for subdomain input along y direction
        integer(iwp), intent(IN) :: js !< end index for subdomain input along y direction
        integer(iwp), intent(IN) :: k1e !< start index for 3rd dimension
        integer(iwp), intent(IN) :: k1s !< end index for 3rd dimension
        integer(iwp), intent(IN) :: k2e !< start index for 4th dimension
        integer(iwp), intent(IN) :: k2s !< end index for 4th dimension
        integer(iwp), intent(IN) :: nbgp !< number of ghost layers of variable to be read

        integer(iwp) :: i !< index along x direction
        integer(iwp) :: id_var !< variable id
        integer(iwp) :: j !< index along y direction
        integer(iwp) :: k1 !< index along 3rd direction
        integer(iwp) :: k2 !< index along 4th direction

        real(wp), dimension(:, :, :, :), allocatable :: tmp !< temporary variable to read data from file according
        !< to its reverse memory access
        real(wp), dimension(k2s:k2e, k1s:k1e, js - nbgp:je + nbgp, is - nbgp:ie + nbgp), intent(INOUT) :: var !< variable to be read

        ! #if defined( __netcdf )
        !-- Inquire variable id
        nc_stat = NF90_INQ_VARID(id, trim(variable_name), id_var)
        !-- Check for collective read-operation and set respective NetCDF flags if required.
        if (collective_read) then
        ! #if defined( __netcdf4_parallel )
        !             nc_stat = NF90_VAR_PAR_ACCESS(id, id_var, NF90_COLLECTIVE)
        ! #endif
        end if

        !-- Allocate temporary variable according to memory access on file.
        allocate (tmp(is:ie, js:je, k1s:k1e, k2s:k2e))
        !-- Get variable
        nc_stat = NF90_GET_VAR(id, id_var, tmp, start=(/is + 1, js + 1, k1s + 1, k2s + 1/), &
                               count=(/ie - is + 1, je - js + 1, k1e - k1s + 1, k2e - k2s + 1/))

        call handle_error('get_variable_4d_real', 535, variable_name)
        !-- Resort data. Please note, dimension subscripts of var all start at 1.
        do i = is, ie
            do j = js, je
                do k1 = k1s, k1e
                    do k2 = k2s, k2e
                        !                var(k2-k2s+1,k1-k1s+1,j-js+1,i-is+1) = tmp(i,j,k1,k2)
                        var(k2, k1, j, i) = tmp(i, j, k1, k2)
                    end do
                end do
            end do
        end do

        deallocate (tmp)

    ! #endif
    end subroutine get_variable_4d_real

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Reads a 4D float variable from file and store it to a 3-d variable.
    !--------------------------------------------------------------------------------------------------!
    subroutine get_variable_4d_to_3d_real(id, variable_name, var, ns, is, ie, js, je, ks, ke, nbgp)

        implicit none

        character(LEN=*) :: variable_name !< variable name

        integer(iwp) :: i !< index along x direction
        integer(iwp), intent(IN) :: id !< file id
        integer(iwp) :: id_var !< variable id
        integer(iwp) :: ie !< end index for subdomain input along x direction
        integer(iwp) :: is !< start index for subdomain input along x direction
        integer(iwp) :: j !< index along y direction
        integer(iwp) :: je !< end index for subdomain input along y direction
        integer(iwp) :: js !< start index for subdomain input along y direction
        integer(iwp) :: k !< index along any 4th dimension
        integer(iwp) :: ke !< end index of 4th dimension
        integer(iwp) :: ks !< start index of 4th dimension
        integer(iwp) :: ns !< start index for subdomain input along n dimension
        integer(iwp), intent(IN) :: nbgp !< number of ghost layers of variable to be read

        real(wp), dimension(:, :, :), allocatable :: tmp !< temporary variable to read data from file according
        !< to its reverse memory access

        real(wp), dimension(:, :, :), intent(INOUT) :: var !< variable where the read data have to be stored:
        !< one dimension is reduced in the process

        ! #if defined( __netcdf )
        !-- To avoid compiler warning about unused variable.
        if (nbgp /= 0) continue
        !-- Inquire variable id
        nc_stat = NF90_INQ_VARID(id, trim(variable_name), id_var)
        !-- Check for collective read-operation and set respective NetCDF flags if required.
        if (collective_read) then
            nc_stat = NF90_VAR_PAR_ACCESS(id, id_var, NF90_COLLECTIVE)
        end if
        !-- Allocate temporary variable according to memory access on file.
        allocate (tmp(is:ie, js:je, ks:ke))
        !-- Get variable
        nc_stat = NF90_GET_VAR(id, id_var, tmp, start=(/is + 1, js + 1, ks + 1, ns + 1/), &
                               count=(/ie - is + 1, je - js + 1, ke - ks + 1, 1/))

        call handle_error('get_variable_4d_to_3d_real', 536, variable_name)
        !-- Resort data. Please note, dimension subscripts of var all start at 1.
        do i = is, ie
            do j = js, je
                do k = ks, ke
                    var(k - ks + 1, j - js + 1, i - is + 1) = tmp(i, j, k)
                end do
            end do
        end do

        deallocate (tmp)

    ! #endif
    end subroutine get_variable_4d_to_3d_real

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Reads a 3D float variables from dynamic driver with the last dimension only having 1 entry
    !> (time,z). Please note, the passed arguments are start indices and number of elements in each
    !> dimension, which is in contrast to the other 3d versions where start- and end indices are passed.
    !> The different handling compared to get_variable_2d_real is due to its different start-index
    !> treatment.
    !--------------------------------------------------------------------------------------------------!
    subroutine get_variable_2d_real_dynamic(id, variable_name, var, i1s, i2s, count_1, count_2)

        implicit none

        character(LEN=*) :: variable_name !< variable name

        integer(iwp) :: count_1 !< number of elements to be read along 1st dimension (with respect to file)
        integer(iwp) :: count_2 !< number of elements to be read along 2nd dimension (with respect to file)
        integer(iwp) :: i1 !< running index along 1st dimension on file
        integer(iwp) :: i1s !< start index for subdomain input along 1st dimension (with respect to file)
        integer(iwp) :: i2 !< running index along 2nd dimension on file
        integer(iwp) :: i2s !< start index for subdomain input along 2nd dimension (with respect to file)
        integer(iwp), intent(IN) :: id !< file id
        integer(iwp) :: id_var !< variable id
        integer(iwp) :: lb1 !< lower bound of 1st dimension (with respect to file)
        integer(iwp) :: lb2 !< lower bound of 2nd dimension (with respect to file)
        integer(iwp) :: ub1 !< upper bound of 1st dimension (with respect to file)
        integer(iwp) :: ub2 !< upper bound of 2nd dimension (with respect to file)

        real(wp), dimension(:, :), allocatable :: tmp !< temporary variable to read data from file according to its reverse memory
        !< access

        real(wp), dimension(:, :, :), intent(INOUT) :: var !< input variable

        ! #if defined( __netcdf )
        !-- Inquire variable id.
        nc_stat = NF90_INQ_VARID(id, trim(variable_name), id_var)
        !-- Allocate temporary variable according to memory access on file.
        !-- Therefore, determine dimension bounds of input array.
        lb1 = lbound(var, 2)
        ub1 = ubound(var, 2)
        lb2 = lbound(var, 1)
        ub2 = ubound(var, 1)

        allocate (tmp(lb1:ub1, lb2:ub2))
        !-- Get variable
        nc_stat = NF90_GET_VAR(id, id_var, tmp, start=(/i1s, i2s/), &
                               count=(/count_1, count_2/))

        call handle_error('get_variable_2d_real_dynamic', 537, variable_name)
        !-- Resort data. Please note, dimension subscripts of var all start at 1.
        do i2 = lb2, ub2
            do i1 = lb1, ub1
                var(i2, i1, 1) = tmp(i1, i2)
            end do
        end do

        deallocate (tmp)

    ! #endif
    end subroutine get_variable_2d_real_dynamic

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Reads time slice from a 2D float variable (t,k) in file
    !--------------------------------------------------------------------------------------------------!
    subroutine get_variable_2d_real_time_slice(ncid, var_name, var_data, t0, ndata, init_data)

        implicit none

        character(LEN=*), intent(IN) :: var_name !< variable name

        integer, intent(IN) :: ncid !< netCDF file handle
        integer, intent(IN) :: ndata !< slice size
        integer, intent(IN) :: t0 !< index for time slice
        integer :: varid !< variable ID

        logical, intent(IN) :: init_data !< if data is to be initiated

        real(KIND=wp), dimension(:), intent(INOUT) :: var_data !< variable data (number of positions of a given time slice)

        if (init_data) var_data = 0.0

        ! #if defined( __netcdf )
        nc_stat = NF90_INQ_VARID(ncid, trim(var_name), varid)
        nc_stat = NF90_GET_VAR(ncid, varid, var_data, start=(/1, t0/), count=(/ndata, 1/))
        call handle_error('get_variable_2d_time_slice', 555, var_name)

    ! #endif
    end subroutine get_Variable_2d_real_time_slice

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Reads a 3D float variables from dynamic driver, such as time-dependent xy-, xz- or yz-boundary
    !> data as well as 3D initialization data. Please note, the passed arguments are start indices and
    !> number of elements in each dimension, which is in contrast to the other 3d versions where start-
    !> and end indices are passed. The different handling of 3D dynamic variables is due to its
    !> asymmetry for the u- and v component.
    !--------------------------------------------------------------------------------------------------!
    subroutine get_variable_3d_real_dynamic(id, variable_name, var, i1s, i2s, i3s, &
                                            count_1, count_2, count_3, par_access)

        implicit none

        character(LEN=*) :: variable_name !< variable name

        integer(iwp) :: count_1 !< number of elements to be read along 1st dimension (with respect to file)
        integer(iwp) :: count_2 !< number of elements to be read along 2nd dimension (with respect to file)
        integer(iwp) :: count_3 !< number of elements to be read along 3rd dimension (with respect to file)
        integer(iwp) :: i1 !< running index along 1st dimension on file
        integer(iwp) :: i1s !< start index for subdomain input along 1st dimension (with respect to file)
        integer(iwp) :: i2 !< running index along 2nd dimension on file
        integer(iwp) :: i2s !< start index for subdomain input along 2nd dimension (with respect to file)
        integer(iwp) :: i3 !< running index along 3rd dimension on file
        integer(iwp) :: i3s !< start index of 3rd dimension, in dynamic file this is either time
        !<(2D boundary) or z (3D)
        integer(iwp), intent(IN) :: id !< file id
        integer(iwp) :: id_var !< variable id
        integer(iwp) :: lb1 !< lower bound of 1st dimension (with respect to file)
        integer(iwp) :: lb2 !< lower bound of 2nd dimension (with respect to file)
        integer(iwp) :: lb3 !< lower bound of 3rd dimension (with respect to file)
        integer(iwp) :: ub1 !< upper bound of 1st dimension (with respect to file)
        integer(iwp) :: ub2 !< upper bound of 2nd dimension (with respect to file)
        integer(iwp) :: ub3 !< upper bound of 3rd dimension (with respect to file)

        logical :: par_access !< additional flag indicating whether parallel read operations should be
        !< performed or not

        real(wp), dimension(:, :, :), allocatable :: tmp !< temporary variable to read data from file according
        !< to its reverse memory access

        real(wp), dimension(:, :, :), intent(INOUT) :: var !< input variable

        ! #if defined( __netcdf )
        !-- Inquire variable id.
        nc_stat = NF90_INQ_VARID(id, trim(variable_name), id_var)
        !-- Check for collective read-operation and set respective NetCDF flags if required.
        !-- Please note, in contrast to the other input routines where each PEs reads its subdomain data,
        !-- dynamic input data not by all PEs, only by those which encompass lateral model boundaries.
        !-- Hence, collective read operations are only enabled for top-boundary data.
        if (collective_read .and. par_access) then
        ! #if defined( __netcdf4_parallel )
        !             nc_stat = NF90_VAR_PAR_ACCESS(id, id_var, NF90_COLLECTIVE)
        ! #endif
        end if
        !-- Allocate temporary variable according to memory access on file.
        !-- Therefore, determine dimension bounds of input array.
        lb1 = lbound(var, 3)
        ub1 = ubound(var, 3)
        lb2 = lbound(var, 2)
        ub2 = ubound(var, 2)
        lb3 = lbound(var, 1)
        ub3 = ubound(var, 1)
        allocate (tmp(lb1:ub1, lb2:ub2, lb3:ub3))
        !-- Get variable
        nc_stat = NF90_GET_VAR(id, id_var, tmp, start=(/i1s, i2s, i3s/), &
                               count=(/count_1, count_2, count_3/))

        call handle_error('get_variable_3d_real_dynamic', 537, variable_name)
        !-- Resort data. Please note, dimension subscripts of var all start at 1.
        do i3 = lb3, ub3
            do i2 = lb2, ub2
                do i1 = lb1, ub1
                    var(i3, i2, i1) = tmp(i1, i2, i3)
                end do
            end do
        end do

        deallocate (tmp)

    ! #endif
    end subroutine get_variable_3d_real_dynamic

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Reads a 4D float variable from dynamic driver, such as time-dependent 3d-data. Please note,
    !> the passed arguments are start indices and number of elements in each dimension, which is in
    !> contrast to some of the 3d versions.
    !--------------------------------------------------------------------------------------------------!
    subroutine get_variable_4d_real_dynamic(id, variable_name, var, i1s, i2s, i3s, i4s, &
                                            count_1, count_2, count_3, count_4, par_access)

        implicit none

        character(LEN=*) :: variable_name !< variable name

        integer(iwp) :: count_1 !< number of elements to be read along 1st dimension (with respect to file)
        integer(iwp) :: count_2 !< number of elements to be read along 2nd dimension (with respect to file)
        integer(iwp) :: count_3 !< number of elements to be read along 3rd dimension (with respect to file)
        integer(iwp) :: count_4 !< number of elements to be read along 4th dimension (with respect to file)
        integer(iwp) :: i1 !< running index along 1st dimension on file
        integer(iwp) :: i1s !< start index for subdomain input along 1st dimension (with respect to file)
        integer(iwp) :: i2 !< running index along 2nd dimension on file
        integer(iwp) :: i2s !< start index for subdomain input along 2nd dimension (with respect to file)
        integer(iwp) :: i3 !< running index along 3rd dimension on file
        integer(iwp) :: i3s !< start index of 3rd dimension
        integer(iwp) :: i4 !< running index along 3rd dimension on file
        integer(iwp) :: i4s !< start index of 3rd dimension
        integer(iwp), intent(IN) :: id !< file id
        integer(iwp) :: id_var !< variable id
        integer(iwp) :: lb1 !< lower bound of 1st dimension (with respect to file)
        integer(iwp) :: lb2 !< lower bound of 2nd dimension (with respect to file)
        integer(iwp) :: lb3 !< lower bound of 3rd dimension (with respect to file)
        integer(iwp) :: lb4 !< lower bound of 4th dimension (with respect to file)
        integer(iwp) :: ub1 !< upper bound of 1st dimension (with respect to file)
        integer(iwp) :: ub2 !< upper bound of 2nd dimension (with respect to file)
        integer(iwp) :: ub3 !< upper bound of 3rd dimension (with respect to file)
        integer(iwp) :: ub4 !< lower bound of 4th dimension (with respect to file)

        logical :: par_access !< additional flag indicating whether parallel read operations should be
        !< performed or not

        real(wp), dimension(:, :, :, :), allocatable :: tmp !< temporary variable to read data from file according
        !< to its reverse memory access

        real(wp), dimension(:, :, :, :), intent(INOUT) :: var !< input variable

        ! #if defined( __netcdf )
        !-- Inquire variable id.
        nc_stat = NF90_INQ_VARID(id, trim(variable_name), id_var)
        !-- Check for collective read-operation and set respective NetCDF flags if required.
        !-- Please note, in contrast to the other input routines where each PEs reads its subdomain data,
        !-- dynamic input data not by all PEs, only by those which encompass lateral model boundaries.
        !-- Hence, collective read operations are only enabled for top-boundary data.
        if (collective_read .and. par_access) then
        ! #if defined( __netcdf4_parallel )
        !             nc_stat = NF90_VAR_PAR_ACCESS(id, id_var, NF90_COLLECTIVE)
        ! #endif
        end if
        !-- Allocate temporary variable according to memory access on file.
        !-- Therefore, determine dimension bounds of input array.
        lb1 = lbound(var, 4)
        ub1 = ubound(var, 4)
        lb2 = lbound(var, 3)
        ub2 = ubound(var, 3)
        lb3 = lbound(var, 2)
        ub3 = ubound(var, 2)
        lb4 = lbound(var, 1)
        ub4 = ubound(var, 1)
        allocate (tmp(lb1:ub1, lb2:ub2, lb3:ub3, lb4:ub4))
        !-- Get variable
        nc_stat = NF90_GET_VAR(id, id_var, tmp, start=(/i1s, i2s, i3s, i4s/), &
                               count=(/count_1, count_2, count_3, count_4/))

        call handle_error('get_variable_4d_real_dynamic', 537, variable_name)
        !-- Resort data. Please note, dimension subscripts of var all start at 1.
        do i4 = lb4, ub4
            do i3 = lb3, ub3
                do i2 = lb2, ub2
                    do i1 = lb1, ub1
                        var(i4, i3, i2, i1) = tmp(i1, i2, i3, i4)
                    end do
                end do
            end do
        end do

        deallocate (tmp)

    ! #endif
    end subroutine get_variable_4d_real_dynamic

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Reads a 5D float variable from file and store it to a 4-d variable.
    !--------------------------------------------------------------------------------------------------!
    subroutine get_variable_5d_to_4d_real(id, variable_name, var, ns, ts, te, is, ie, js, je, ks, &
                                          ke, nbgp)

        implicit none

        character(LEN=*) :: variable_name !< variable name

        integer(iwp) :: i !< index along x direction
        integer(iwp), intent(IN) :: id !< file id
        integer(iwp) :: id_var !< variable id
        integer(iwp) :: ie !< end index for subdomain input along x direction
        integer(iwp) :: is !< start index for subdomain input along x direction
        integer(iwp) :: j !< index along y direction
        integer(iwp) :: je !< end index for subdomain input along y direction
        integer(iwp) :: js !< start index for subdomain input along y direction
        integer(iwp) :: k !< index along any 5th dimension
        integer(iwp) :: ke !< end index of 5th dimension
        integer(iwp) :: ks !< start index of 5th dimension
        integer(iwp) :: ns !< start index for subdomain input along n dimension: ns coincides here with
        !< ne, since, we select only one value along the 1st dimension n
        integer(iwp) :: t !< index along t direction
        integer(iwp) :: te !< end index for subdomain input along t direction
        integer(iwp) :: ts !< start index for subdomain input along t direction
        integer(iwp), intent(IN) :: nbgp !< number of ghost layers of variable to be read

        real(wp), dimension(:, :, :, :), allocatable :: tmp !< temporary variable to read data from file according to its reverse
        !< memory access
        real(wp), dimension(:, :, :, :), intent(INOUT) :: var !< variable to be read

        ! #if defined( __netcdf )
        !-- To avoid compiler warning about unused variable.
        if (nbgp /= 0) continue
        !-- Inquire variable id
        nc_stat = NF90_INQ_VARID(id, trim(variable_name), id_var)
        !-- Check for collective read-operation and set respective NetCDF flags if required.
        if (collective_read) then
            nc_stat = NF90_VAR_PAR_ACCESS(id, id_var, NF90_COLLECTIVE)
        end if
        !-- Allocate temporary variable according to memory access on file.
        allocate (tmp(ks:ke, js:je, is:is, ts:te))
        !-- Get variable
        nc_stat = NF90_GET_VAR(id, id_var, tmp, start=(/ks + 1, js + 1, is + 1, ts + 1, ns/), &
                               count=(/ke - ks + 1, je - js + 1, ie - is + 1, te - ts + 1, 1/))

        call handle_error('get_variable_5d_to_4d_real', 538, variable_name)
        !-- Resort data. Please note, dimension subscripts of var all start at 1.
        do t = ts, te
            do i = is, ie
                do j = js, je
                    do k = ks, ke
                        var(t - ts + 1, i - is + 1, j - js + 1, k - ks + 1) = tmp(k, j, i, t)
                    end do
                end do
            end do
        end do

        deallocate (tmp)

    ! #endif
    end subroutine get_variable_5d_to_4d_real

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Reads a 5D float variable from file.
    !> Note - This subroutine is used specific for reading NC variable emission_values having a "z"
    !>        dimension. Mentioned dimension is to be removed in the future and this subroutine shall
    !>        be depreciated accordingly.
    !--------------------------------------------------------------------------------------------------!
    subroutine get_variable_5d_real(id, variable_name, var, is, ie, js, je, k1s, k1e, k2s, k2e, k3s, &
                                    k3e)

        implicit none

        character(LEN=*) :: variable_name !< variable name

        integer(iwp) :: i !< i index
        integer(iwp), intent(IN) :: id !< netCDF file ID (ncid)
        integer(iwp) :: id_var !< netCDF variable ID (varid)
        integer(iwp) :: ie !< i index start
        integer(iwp) :: is !< i index end
        integer(iwp) :: j !< j index
        integer(iwp) :: je !< j index start
        integer(iwp) :: js !< j index end
        integer(iwp) :: k1 !< k1 index
        integer(iwp) :: k1e !< k1 index start
        integer(iwp) :: k1s !< k1 index end
        integer(iwp) :: k2 !< k2 index
        integer(iwp) :: k2e !< k2 index start
        integer(iwp) :: k2s !< k2 index end
        integer(iwp) :: k3 !< k3 index
        integer(iwp) :: k3e !< k3 index start
        integer(iwp) :: k3s !< k3 index end

        real(wp), dimension(:, :, :, :, :), allocatable :: tmp !< temp array to read data from file
        real(wp), dimension(:, :, :, :, :), intent(INOUT) :: var !< variable to be read


        ! #if defined( __netcdf )
        !-- Inquire variable id
        nc_stat = NF90_INQ_VARID(id, trim(variable_name), id_var)

        !-- Check for collective read-operation and set respective NetCDF flags if required.
        if (collective_read) then


        ! #if defined( __netcdf4_parallel )
        !             nc_stat = NF90_VAR_PAR_ACCESS(id, id_var, NF90_COLLECTIVE)
        ! #endif
        end if

        !-- Allocate temporary variable according to memory access on file.
        allocate (tmp(is:ie, js:je, k1s:k1e, k2s:k2e, k3s:k3e))

        !-- Get variable from file
        nc_stat = NF90_GET_VAR(id, id_var, tmp, start=(/is + 1, js + 1, k1s + 1, k2s + 1, k3s + 1/), &
                               count=(/ie - is + 1, je - js + 1, k1e - k1s + 1, k2e - k2s + 1, k3e - k3s + 1/))

        call handle_error('get_variable_5d_real', 535, variable_name)

        !-- Resort (reverse index order) and standardize (from 1 to N) output array
        do i = is, ie
            do j = js, je
                do k1 = k1s, k1e
                    do k2 = k2s, k2e
                        do k3 = k3s, k3e
                            var(k3 - k3s + 1, k2 - k2s + 1, k1 - k1s + 1, j - js + 1, i - is + 1) = tmp(i, j, k1, k2, k3)
                        end do
                    end do
                end do
            end do
        end do

        deallocate (tmp)

    ! #endif
    end subroutine get_variable_5d_real

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Reads a 5D float variables for chemistry initializat. Please note, the passed arguments are
    !> start indices and number of elements in each dimension, which is in contrast to the other 3d
    !> versions where start- and end indices are passed.
    !> Note(1) - This subroutine is more flexible than get_variable_xd_real as it provides much better
    !>           control over starting and count indices.
    !> Note(2) - This subroutine is used specific for reading NC variable emission_values having a "z"
    !>           dimension. Mentioned dimension is to be removed in the future and this subroutine shall
    !>           be depreciated accordingly
    !--------------------------------------------------------------------------------------------------!
    subroutine get_variable_5d_real_dynamic(id, variable_name, var, i1s, i2s, i3s, i4s, i5s, &
                                            count_1, count_2, count_3, count_4, count_5, par_access)

        implicit none

        character(LEN=*) :: variable_name !< variable name

        integer(iwp) :: count_1 !< # elements read in dimension 1 wrt file
        integer(iwp) :: count_2 !< # elements read in dimension 2 wrt file
        integer(iwp) :: count_3 !< # elements read in dimension 3 wrt file
        integer(iwp) :: count_4 !< # elements read in dimension 4 wrt file
        integer(iwp) :: count_5 !< # elements read in dimension 5 wrt file
        integer(iwp) :: i1 !< index for dimension 1 on file
        integer(iwp) :: i1s !< starting index for dimension 1 hyperslab
        integer(iwp) :: i2 !< index for dimension 2 on file
        integer(iwp) :: i2s !< starting index for dimension 2 hyperslab
        integer(iwp) :: i3 !< index for dimension 3 on file
        integer(iwp) :: i3s !< starting index for dimension 3 hyperslab
        integer(iwp) :: i4 !< index for dimension 4 on file
        integer(iwp) :: i4s !< starting index for dimension 4 hyperslab
        integer(iwp) :: i5 !< index for dimension 5 on file
        integer(iwp) :: i5s !< starting index for dimension 5 hyperslab
        integer(iwp), intent(IN) :: id !< netCDF file id (ncid)
        integer(iwp) :: id_var !< netCDF variable id (varid)
        integer(iwp) :: lb1 !< lower bound of dimension 1 wrt file
        integer(iwp) :: lb2 !< lower bound of dimension 2 wrt file
        integer(iwp) :: lb3 !< lower bound of dimension 3 wrt file
        integer(iwp) :: lb4 !< lower bound of dimension 4 wrt file
        integer(iwp) :: lb5 !< lower bound of dimension 5 wrt file
        integer(iwp) :: ub1 !< upper bound of dimension 1 wrt file
        integer(iwp) :: ub2 !< upper bound of dimension 2 wrt file
        integer(iwp) :: ub3 !< upper bound of dimension 3 wrt file
        integer(iwp) :: ub4 !< upper bound of dimension 4 wrt file
        integer(iwp) :: ub5 !< upper bound of dimension 5 wrt file

        logical :: par_access !< additional flag indicating parallel read

        real(wp), dimension(:, :, :, :, :), allocatable :: tmp !< temporary variable to read data from file according to its reverse
        !< array index order
        real(wp), dimension(:, :, :, :, :), intent(INOUT) :: var !< input variable

        ! #if defined( __netcdf )
        !-- Inquire variable id.
        nc_stat = NF90_INQ_VARID(id, trim(variable_name), id_var)

        !-- Check for collective read-operation and set respective NetCDF flags if required.
        !-- Please note, in contrast to the other input routines where each PEs reads its subdomain data,
        !-- dynamic input data not by all PEs, only by those which encompass lateral model boundaries.
        !-- Hence, collective read operations are only enabled for top-boundary data.
        if (collective_read .and. par_access) then


        ! #if defined( __netcdf4_parallel )
        !             nc_stat = NF90_VAR_PAR_ACCESS(id, id_var, NF90_COLLECTIVE)
        ! #endif
        end if

        !-- Allocate temporary variable according to memory access on file.
        !-- Therefore, determine dimension bounds of input array.
        lb1 = lbound(var, 5)
        ub1 = ubound(var, 5)
        lb2 = lbound(var, 4)
        ub2 = ubound(var, 4)
        lb3 = lbound(var, 3)
        ub3 = ubound(var, 3)
        lb4 = lbound(var, 2)
        ub4 = ubound(var, 2)
        lb5 = lbound(var, 1)
        ub5 = ubound(var, 1)
        allocate (tmp(lb1:ub1, lb2:ub2, lb3:ub3, lb4:ub4, lb5:ub5))

        !-- Get variable
        nc_stat = NF90_GET_VAR(id, id_var, tmp, start=(/i1s, i2s, i3s, i4s, i5s/), &
                               count=(/count_1, count_2, count_3, count_4, count_5/))

        call handle_error('get_variable_3d_real_dynamic', 537, variable_name)

        !-- Assign temp array to output.  Note reverse index order
        do i5 = lb5, ub5
            do i4 = lb4, ub4
                do i3 = lb3, ub3
                    do i2 = lb2, ub2
                        do i1 = lb1, ub1
                            var(i5, i4, i3, i2, i1) = tmp(i1, i2, i3, i4, i5)
                        end do
                    end do
                end do
            end do
        end do

        deallocate (tmp)

    ! #endif
    end subroutine get_variable_5d_real_dynamic

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Inquires the number of variables in a file
    !--------------------------------------------------------------------------------------------------!
    subroutine inquire_num_variables(id, num_vars)

        implicit none

        integer(iwp), intent(IN) :: id !< file id
        integer(iwp), intent(INOUT) :: num_vars !< number of variables in a file

        ! #if defined( __netcdf )
        nc_stat = NF90_INQUIRE(id, NVARIABLES=num_vars)
        call handle_error('inquire_num_variables', 539)

    ! #endif
    end subroutine inquire_num_variables

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Inquires the variable names belonging to a file.
    !--------------------------------------------------------------------------------------------------!
    subroutine inquire_variable_names(id, var_names)

        implicit none

        character(LEN=*), dimension(:), intent(INOUT) :: var_names !< return variable - variable names

        integer(iwp) :: i !< loop variable
        integer(iwp), intent(IN) :: id !< file id
        integer(iwp) :: num_vars !< number of variables (unused return parameter)

        integer(iwp), dimension(:), allocatable :: varids !< dummy array to strore variable ids temporarily

        ! #if defined( __netcdf )
        allocate (varids(1:size(var_names)))
        nc_stat = NF90_INQ_VARIDS(id, NVARS=num_vars, VARIDS=varids)
        call handle_error('inquire_variable_names', 540)

        do i = 1, size(var_names)
            nc_stat = NF90_INQUIRE_VARIABLE(id, varids(i), NAME=var_names(i))
            call handle_error('inquire_variable_names', 540)
        end do

        deallocate (varids)

    ! #endif
    end subroutine inquire_variable_names

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Inquires the _FillValue settings of an integer variable.
    !--------------------------------------------------------------------------------------------------!
    subroutine inquire_fill_value_int(id, var_name, no_fill, fill_value)

        character(LEN=*), intent(IN) :: var_name !< variable name

        integer(iwp), intent(IN) :: id !< file id
        integer(iwp) :: id_var !< netCDF variable id (varid)
        integer(iwp) :: fill_value !< fill value
        integer(iwp) :: no_fill !< flag indicating whether fill values are set or not

        ! #if defined( __netcdf )
        nc_stat = NF90_INQ_VARID(id, trim(var_name), id_var)
        nc_stat = NF90_INQ_VAR_FILL(id, id_var, no_fill, fill_value)
        ! #endif
        !-- Further line is just to avoid compiler warnings. no_fill might be used in future.
        if (no_fill == 0 .or. no_fill /= 0) continue

    end subroutine inquire_fill_value_int

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Inquires the _FillValue settings of a real variable.
    !--------------------------------------------------------------------------------------------------!
    subroutine inquire_fill_value_real(id, var_name, no_fill, fill_value)

        character(LEN=*), intent(IN) :: var_name !< variable name

        integer(iwp), intent(IN) :: id !< file id
        integer(iwp) :: id_var !< netCDF variable id (varid)
        integer(iwp) :: no_fill !< flag indicating whether fill values are set or not

        ! #if defined( __imuk_old )
        !         integer(iwp) :: fill_value_int !< fill value workaround
        ! #endif
        real(wp), intent(OUT) :: fill_value !< fill value

        ! #if defined( __netcdf )
        nc_stat = NF90_INQ_VARID(id, trim(var_name), id_var)
        ! #if defined( __imuk_old )
        !         nc_stat = NF90_INQ_VAR_FILL(id, id_var, no_fill, fill_value_int)
        !         fill_value = fill_value_int
        ! #else
        nc_stat = NF90_INQ_VAR_FILL(id, id_var, no_fill, fill_value)
        ! #endif
        ! #endif
        !-- Further line is just to avoid compiler warnings. no_fill might be used in future.
        if (no_fill == 0 .or. no_fill /= 0) continue

    end subroutine inquire_fill_value_real

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Modify a 3D array to that all fill values that are found are overwritten with meaningful data.
    !> This ignores the presence of any topography (topo_flags) and fills the 3D array completely.
    !--------------------------------------------------------------------------------------------------!
    subroutine overwrite_fill_values(var, is, ie, js, je, ks, ke, fill_value, default_valid_value)

        implicit none

        integer(iwp) :: i !< index along x direction
        integer(iwp), intent(IN) :: ie !< start index for subdomain input along x direction
        integer(iwp), intent(IN) :: is !< end index for subdomain input along x direction
        integer(iwp) :: j !< index along y direction
        integer(iwp), intent(IN) :: je !< start index for subdomain input along y direction
        integer(iwp), intent(IN) :: js !< end index for subdomain input along y direction
        integer(iwp) :: k !< index along any 3rd dimension
        integer(iwp), intent(IN) :: ke !< start index of 3rd dimension
        integer(iwp), intent(IN) :: ks !< end index of 3rd dimension

        real(wp), intent(IN) :: fill_value !< fill value to look for and eliminate
        real(wp), intent(IN) :: default_valid_value !< default valid value as given by the subroutine call
        real(wp) :: nearest_valid_value !< stored valid value (taken from gridpoints above)

        real(wp), dimension(ks:ke, js:je, is:ie), intent(INOUT) :: var !< variable to be read

        nearest_valid_value = default_valid_value
        !-- Fill in data from top to bottom using the closest value above that is not a fill value.
        do i = is, ie
            do j = js, je
                do k = ke, ks, -1
                    if (var(k, j, i) /= fill_value) then
                        nearest_valid_value = var(k, j, i)
                    else
                        var(k, j, i) = nearest_valid_value
                    end if
                end do
            end do
        end do

    end subroutine overwrite_fill_values

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Prints out a text message corresponding to the current status.
    !--------------------------------------------------------------------------------------------------!
    subroutine handle_error(routine_name, errno, name)

        implicit none

        character(LEN=7) :: message_identifier !< string for the error number
        character(LEN=*), optional :: name !< name of variable where reading failed
        character(LEN=*) :: routine_name !< routine name where the error happened

        integer(iwp) :: errno

        ! #if defined( __netcdf )
        if (nc_stat /= NF90_NOERR) then

            write (*, '(''NCF'',I4.4)') errno
            print *, "Problem reading attribute/variable - " // &
                     trim(name) // ": " // trim(NF90_STRERROR(nc_stat))

            if (present(name)) then
                print *, "Problem reading attribute/variable - " // &
                         trim(name) // ": " // trim(NF90_STRERROR(nc_stat))
            else
                print *, trim(NF90_STRERROR(nc_stat))
            end if
        end if

    ! #endif
    end subroutine handle_error

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! -------------------------------------------------------------------------------------------------!
    !> Determine list of all available building IDs within the domain.
    !--------------------------------------------------------------------------------------------------!
    subroutine list_building_ids(build_ids_unique)

        integer(iwp) :: i !< running index along x-direction
        integer(iwp) :: j !< running index along y-direction
        integer(iwp) :: nr !< index variable indication maximum terrain height for respective building ID
        integer(iwp) :: num_build !< counter for number of buildings

        ! #if defined( __parallel )
        !         integer(iwp), dimension(:), allocatable :: displace_dum !< displacements of start addresses, used for MPI_ALLGATHERV
        ! #endif
        integer(iwp), dimension(:), allocatable :: build_ids !< building IDs on entire model domain
        integer(iwp), dimension(:), allocatable :: build_ids_unique !< temporary array used for resizing
        integer(iwp), dimension(:), allocatable :: build_ids_unique_tmp !< temporary array used for resizing
        integer(iwp), dimension(:), allocatable :: build_ids_l !< building IDs on local subdomain
        integer(iwp), dimension(:), allocatable :: build_ids_l_tmp !< temporary array used to resize array of building IDs

        integer(iwp), dimension(0:numprocs - 1) :: num_buildings !< number of buildings with different ID on entire model domain

        if (.not. building_id_f % from_file) return

        num_buildings = 0
        !-- Allocate at least one element for building ids and give it an inital negative value that
        !-- will be overwritten later. This, however, is necessary in case there all IDs in the model
        !-- domain are fill values.
        allocate (build_ids_l(1))
        build_ids_l = -1
        do i = nxl, nxr
            do j = nys, nyn
                if (building_id_f % var(j, i) /= building_id_f % fill) then
                    if (num_buildings(myid) > 0) then
                        if (any(building_id_f % var(j, i) == build_ids_l)) then
                            cycle
                        else
                            num_buildings(myid) = num_buildings(myid) + 1
                            !--                Resize array with different local building ids
                            allocate (build_ids_l_tmp(1:size(build_ids_l)))
                            build_ids_l_tmp = build_ids_l
                            deallocate (build_ids_l)
                            allocate (build_ids_l(1:num_buildings(myid)))
                            build_ids_l(1:num_buildings(myid) - 1) = build_ids_l_tmp(1:num_buildings(myid) - 1)
                            build_ids_l(num_buildings(myid)) = building_id_f % var(j, i)
                            deallocate (build_ids_l_tmp)
                        end if
                    !--          First occuring building id on PE.
                    else
                        num_buildings(myid) = num_buildings(myid) + 1
                        build_ids_l(1) = building_id_f % var(j, i)
                    end if
                end if
            end do
        end do
        !-- Determine number of different building ids for the entire domain.
        ! #if defined( __parallel )
        !         call MPI_ALLREDUCE(MPI_IN_PLACE, num_buildings, numprocs, MPI_INTEGER, MPI_SUM, comm2d, ierr)
        ! #endif
        !-- Gather all buildings ids on each PEs.
        !-- First, allocate array encompassing all building ids in model domain.
        allocate (build_ids(1:sum(num_buildings)))


        ! #if defined( __parallel )
        !         !-- Allocate array for displacements.
        !         !-- As each PE may has a different number of buildings, so that the block sizes send by each
        !         !-- PE may not be equal. Hence,  information about the respective displacement is required,
        !         !-- indicating the respective adress where each MPI-task writes into the receive buffer array.
        !         allocate (displace_dum(0:numprocs - 1))
        !         displace_dum(0) = 0
        !         do i = 1, numprocs - 1
        !             displace_dum(i) = displace_dum(i - 1) + num_buildings(i - 1)
        !         end do
        !         call MPI_ALLGATHERV(build_ids_l(1:num_buildings(myid)), num_buildings(myid), &
        !                             MPI_INTEGER, build_ids, num_buildings, displace_dum, MPI_INTEGER, &
        !                             comm2d, ierr)
        !         deallocate (displace_dum)
        ! #else
        build_ids = build_ids_l

        ! #endif
        deallocate (build_ids_l)
        !-- Note, in parallel mode building ids can occure mutliple times, as each PE has send its own
        !-- ids. Therefore, sort out building ids which appear more than one time.
        num_build = 0
        do nr = 1, size(build_ids)

            if (allocated(build_ids_unique)) then
                if (any(build_ids(nr) == build_ids_unique)) then
                    cycle
                else
                    num_build = num_build + 1
                    !--          Resize.
                    allocate (build_ids_unique_tmp(1:num_build))
                    build_ids_unique_tmp(1:num_build - 1) = build_ids_unique(1:num_build - 1)
                    deallocate (build_ids_unique)
                    allocate (build_ids_unique(1:num_build))
                    build_ids_unique(1:num_build - 1) = build_ids_unique_tmp(1:num_build - 1)
                    build_ids_unique(num_build) = build_ids(nr)
                    deallocate (build_ids_unique_tmp)
                end if
            else
                num_build = num_build + 1
                allocate (build_ids_unique(1:num_build))
                build_ids_unique(num_build) = build_ids(nr)
            end if
        end do

    end subroutine list_building_ids

end module netcdf_data_input_mod
