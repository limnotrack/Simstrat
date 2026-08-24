! ---------------------------------------------------------------------------------
!     Simstrat a physical 1D model for lakes and reservoirs
!
!     Developed by:  Group of Applied System Analysis
!                    Dept. of Surface Waters - Research and Management
!                    Eawag - Swiss Federal institute of Aquatic Science and Technology
!
!     Copyright (C) 2026, Eawag
!
!
!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program.  If not, see <http://www.gnu.org/licenses/>.
! ---------------------------------------------------------------------------------
!<    +---------------------------------------------------------------+
!     |  Simstrat - AED interface (libaed-api generation)
!     |
!     |  Successor of the simstrat_aed2 module: instead of driving the
!     |  old libaed2 routines directly, this module couples Simstrat to
!     |  the generic host interface provided by libaed-api. The API owns
!     |  light attenuation, mobility (via a host callback), equilibration,
!     |  state repair and the benthic/surface/pelagic kinetics including
!     |  flux integration; Simstrat keeps vertical diffusion, advection
!     |  (in strat_lateral/strat_advection) and output.
!     |
!     |  The API stores state as cc(n_vars, n_layers) with benthic state
!     |  in a separate cc_hz vector, while the rest of Simstrat expects
!     |  the historical AED2_state(n_layers, n_vars+n_vars_ben) layout.
!     |  To keep strat_outputfile/strat_lateral/strat_advection untouched
!     |  this module maintains mirror arrays in the old layout and
!     |  transposes at the interface each timestep (single column, cheap).
!<    +---------------------------------------------------------------+

module simstrat_aed
   use strat_simdata
   use strat_grid
   use strat_solver
   use utilities
   use aed_api
   use aed_core, only: aed_variable_t, aed_get_var, V_STATE, V_DIAGNOSTIC, V_EXTERNAL
   use aed_zones, only: aed_init_zones, aedZones, aed_n_zones, api_set_zone_funcs, &
                         api_zone_t, calc_zone_areas_t, copy_to_zone_t, copy_from_zone_t
   use, intrinsic :: iso_c_binding, only: C_INT32_T, C_DOUBLE
   use, intrinsic :: ieee_arithmetic

   implicit none
   private

   type, public :: SimstratAED
      class(AEDConfig), pointer :: aed_cfg
      class(StaggeredGrid), pointer :: grid

      ! Variable counts (n_aed_vars is the total before the API registers
      ! the environment globals; all state/diag variables live in 1..n_aed_vars)
      integer :: n_aed_vars, n_vars, n_vars_ben, n_vars_diag, n_vars_diag_sheet

      ! State/diagnostic arrays in API layout (var, layer)
      real(RK), pointer, dimension(:,:) :: cc            ! (n_vars, nz_grid)
      real(RK), pointer, dimension(:)   :: cc_hz         ! (n_vars_ben)
      real(RK), pointer, dimension(:,:) :: cc_diag       ! (n_vars_diag, nz_grid)
      real(RK), pointer, dimension(:)   :: cc_diag_hz    ! (n_vars_diag_sheet)

      ! Mirrors in historical Simstrat layout (layer, var), exposed via ModelState
      real(RK), pointer, dimension(:,:) :: mirror_state  ! (nz_grid, n_vars+n_vars_ben)
      real(RK), pointer, dimension(:,:) :: mirror_diag   ! (nz_grid, n_vars_diag)
      real(RK), pointer, dimension(:)   :: mirror_diag_sheet ! (n_vars_diag_sheet)

      ! Environment targets handed to the API as pointers. Scalar values that
      ! are not pointer components of ModelState are copied here every step.
      real(RK), pointer :: dt_c, yearday_c, air_temp_c, air_pres_c, longwave_c
      real(RK), pointer :: blueice_c, whiteice_c, colarea_c, sedzone_c, kw_c, lat_c, lon_c
      integer, pointer :: bot_i
      logical, pointer :: active_c

      real(RK), pointer, dimension(:) :: par_c, rad_c, pres_c, tss_c
      real(RK), pointer, dimension(:) :: nir_c, uva_c, uvb_c
      real(RK), pointer, dimension(:) :: heights_c, bioext_c, biodrag_c
      real(RK), pointer, dimension(:) :: sedzones_c, zeros_c

      ! Variable names (pointed to from ModelState for output). names holds
      ! pelagic AND benthic names concatenated (size n_vars+n_vars_ben),
      ! matching AED2_state's columns 1:1; bennames aliases its tail.
      character(len=48), pointer :: names(:)
      character(len=48), pointer :: bennames(:)
      character(len=48), pointer :: diagnames(:)
      character(len=48), pointer :: diagnames_sheet(:)

      ! Sediment zones (GLM-style depth bands, only used if benthic_mode > 1).
      ! zone_cc/zone_cc_hz/zone_diag/zone_diag_hz back aed_init_zones; their
      ! per-zone slices are aliased inside libaed-api's aedZones(:) array.
      integer :: n_zones = 0
      real(RK), pointer, dimension(:) :: zone_heights_m ! cumulative height above bottom per zone
      real(RK), pointer, dimension(:,:,:) :: zone_cc, zone_diag
      real(RK), pointer, dimension(:,:) :: zone_cc_hz, zone_diag_hz
      ! Per-zone benthic state/diagnostics for output, (n_zones, n_vars_ben) /
      ! (n_zones, n_vars_diag_sheet) -- transposed re-layout of zone_cc_hz/
      ! zone_diag_hz (which are variable-major), refreshed each step and
      ! exposed via ModelState for dedicated <name>_zone_out.dat files.
      real(RK), pointer, dimension(:,:) :: zone_state_out, zone_diag_sheet_out

   contains
      procedure, pass(self), public :: init
      procedure, pass(self), public :: update
   end type SimstratAED

contains

   ! Called once during Simstrat initialization: configures the AED models from
   ! the namelist file, hands environment and data pointers to libaed-api and
   ! sets the initial conditions.

   subroutine init(self, state, grid, model_cfg, aed_cfg, model_param)
      implicit none

      ! Arguments
      class(SimstratAED) :: self
      class(ModelState) :: state
      class(StaggeredGrid), target :: grid
      class(ModelConfig), target :: model_cfg
      class(AEDConfig), target :: aed_cfg
      class(ModelParam), target :: model_param

      ! Local variables
      type(aed_coupling_t) :: cpl
      type(aed_env_t) :: env(1)
      type(aed_data_t) :: dat(1)
      procedure(aed_mobility_fn_t), pointer :: pmob
      type(aed_variable_t), pointer :: tvar
      integer :: av, v, sv, benthic_mode

      self%grid => grid
      self%aed_cfg => aed_cfg

      associate (n_vars => self%n_vars, &
                 n_vars_ben => self%n_vars_ben, &
                 n_vars_diag => self%n_vars_diag, &
                 n_vars_diag_sheet => self%n_vars_diag_sheet)

         ! ------------------------------------------------------------------
         ! (1) Parse the AED models configuration (&aed_models namelist)
         ! ------------------------------------------------------------------
         self%n_aed_vars = aed_configure_models(aed_cfg%aed_config_file, &
                                                n_vars, n_vars_ben, n_vars_diag, n_vars_diag_sheet)

         print "(/,5X,'AED : n_aed_vars  = ',I4,' ; MaxLayers         = ',I4)", self%n_aed_vars, grid%nz_grid
         print "(  5X,'AED : n_vars      = ',I4,' ; n_vars_ben        = ',I4)", n_vars, n_vars_ben
         print "(  5X,'AED : n_vars_diag = ',I4,' ; n_vars_diag_sheet = ',I4,/)", n_vars_diag, n_vars_diag_sheet

         ! ------------------------------------------------------------------
         ! (2) Allocate coupler memory (API arrays, mirrors, env targets)
         ! ------------------------------------------------------------------
         call allocate_memory(self)

         self%air_pres_c = model_param%p_air*100._RK   ! [mbar] -> [Pa]
         self%lat_c = model_param%Lat*pi/180._RK       ! [deg]  -> [rad]
         self%lon_c = model_param%Lon*pi/180._RK       ! [deg]  -> [rad]
         self%kw_c = aed_cfg%background_extinction
         self%dt_c = state%dt
         self%yearday_c = day_of_year(state%current_year, state%current_month, state%current_day)

         ! ------------------------------------------------------------------
         ! (3) Coupling options
         ! ------------------------------------------------------------------
         benthic_mode = aed_cfg%benthic_mode
         if (benthic_mode > 1) then
            if (aed_cfg%n_zones < 1 .or. .not. allocated(aed_cfg%zone_heights)) then
               call warn('Sediment zones (benthic_mode > 1) require AEDConfig.NZones and AEDConfig.ZoneHeights. Using benthic_mode = 1.')
               benthic_mode = 1
            end if
         end if

         cpl%mobility_off = .not. aed_cfg%particle_mobility
         cpl%bioshade_feedback = aed_cfg%bioshade_feedback
         cpl%repair_state = .true.
         cpl%split_factor = 1
         cpl%benthic_mode = benthic_mode
         ! Simstrat provides the in-water light field itself (from rad_vol), as
         ! the old coupler did; the API then skips its internal Light routine.
         cpl%link_ext_par = .true.
         cpl%link_rain_loss = .false.
         cpl%link_solar_shade = .false.
         cpl%link_bottom_drag = .false.
         cpl%link_water_clarity = .false.
         cpl%do_particle_bgc = .false.
         ! Simstrat's column is GLM-like (bottom = index 1, cumulative layer
         ! areas): the GLM-style benthic path reproduces the flank distribution
         ! of benthic fluxes (benthic_mode = 1) of the old Simstrat-AED2 coupler.
         cpl%glm_style_zones = .true.
         cpl%Kw => self%kw_c

         call aed_set_coupling(cpl)

         ! ------------------------------------------------------------------
         ! (4) Environment: hand pointers into Simstrat state/grid to the API
         ! ------------------------------------------------------------------
         env(1)%n_layers = grid%nz_grid
         env(1)%top_idx => grid%nz_occupied
         env(1)%bot_idx => self%bot_i
         env(1)%active => self%active_c

         env(1)%timestep => self%dt_c
         env(1)%yearday => self%yearday_c
         env(1)%latitude => self%lat_c
         env(1)%longitude => self%lon_c

         ! Above water (meteorological) conditions
         env(1)%longwave => self%longwave_c
         env(1)%air_temp => self%air_temp_c
         env(1)%air_pres => self%air_pres_c
         env(1)%wind => state%uv10
         env(1)%rain => state%rain
         env(1)%I_0 => state%rad0

         ! Column/grid information
         env(1)%col_depth => grid%lake_level
         env(1)%col_area => self%colarea_c
         env(1)%height => self%heights_c
         env(1)%depth => grid%layer_depth(1:grid%nz_grid)
         env(1)%area => grid%Az_vol
         env(1)%dz => grid%h(1:grid%nz_grid)

         ! Water column physics
         env(1)%temp => state%T
         env(1)%salt => state%S
         env(1)%rho => state%rho
         env(1)%cvel => self%zeros_c
         env(1)%pres => self%pres_c
         env(1)%rad => self%rad_c

         ! Light
         env(1)%extc => state%absorb_vol
         env(1)%par => self%par_c
         env(1)%nir => self%nir_c
         env(1)%uva => self%uva_c
         env(1)%uvb => self%uvb_c

         ! Turbidity (no external solids model in Simstrat)
         env(1)%tss => self%tss_c
         env(1)%ss1 => self%zeros_c
         env(1)%ss2 => self%zeros_c
         env(1)%ss3 => self%zeros_c
         env(1)%ss4 => self%zeros_c

         ! Benthic environment
         env(1)%ustar_bed => self%zeros_c
         env(1)%wv_uorb => self%zeros_c
         env(1)%wv_t => self%zeros_c
         env(1)%layer_stress => state%u_taub
         env(1)%delzBlueIce => self%blueice_c
         env(1)%delzWhiteIce => self%whiteice_c
         env(1)%sed_zones => self%sedzones_c
         env(1)%sed_zone => self%sedzone_c

         ! Feedback arrays read back by Simstrat after each AED step
         env(1)%biodrag => self%biodrag_c
         env(1)%bioextc => self%bioext_c

         call aed_set_model_env(env, 1, grid%nz_grid)

         ! ------------------------------------------------------------------
         ! (4b) Sediment zones (must be registered before aed_set_model_data,
         !      which sizes its zone-dependent arrays off aed_n_zones)
         ! ------------------------------------------------------------------
         if (benthic_mode > 1) call init_zones(self, aed_cfg)

         ! ------------------------------------------------------------------
         ! (5) Data arrays (the API initialises them to module default values)
         ! ------------------------------------------------------------------
         dat(1)%cc => self%cc
         dat(1)%cc_hz => self%cc_hz
         dat(1)%cc_diag => self%cc_diag
         dat(1)%cc_diag_hz => self%cc_diag_hz

         call aed_set_model_data(dat, 1, grid%nz_grid)

         ! Verify all variable dependencies requested by the AED modules are met
         call aed_check_model_setup

         ! Register Simstrat's (GLM-derived) settling/rising scheme with the API
         pmob => mobility
         call aed_set_mobility_fn(pmob)

         ! ------------------------------------------------------------------
         ! (6) Variable names and ModelState links (output, inflow, advection)
         ! ------------------------------------------------------------------
         call assign_var_names(self)

         state%AED2_state => self%mirror_state
         state%AED2_diagnostic => self%mirror_diag
         state%AED2_diagnostic_sheet => self%mirror_diag_sheet
         state%n_AED2_state = n_vars + n_vars_ben
         state%n_AED2_diagnostic = n_vars_diag
         state%n_AED2_diagnostic_sheet = n_vars_diag_sheet
         state%AED2_state_names => self%names
         state%AED2_diagnostic_names => self%diagnames
         state%AED2_diagnostic_names_sheet => self%diagnames_sheet

         if (self%n_zones > 0) then
            state%n_AED_zones = self%n_zones
            state%AED_zone_heights => self%zone_heights_m
            state%AED_zone_state => self%zone_state_out
            state%AED_zone_diagnostic_sheet => self%zone_diag_sheet_out
            state%AED_zone_state_names => self%bennames
         end if

         ! ------------------------------------------------------------------
         ! (7) Initial conditions: default values from the namelist were set by
         !     aed_set_model_data; override them from <var>_ini.dat files where
         !     present (same convention as the old coupler).
         ! ------------------------------------------------------------------
         ! Start from the API defaults
         do v = 1, n_vars
            self%mirror_state(:, v) = self%cc(v, :)
         end do
         do sv = 1, n_vars_ben
            self%mirror_state(:, n_vars + sv) = self%cc_hz(sv)
         end do

         v = 0; sv = 0
         do av = 1, self%n_aed_vars
            if (.not. aed_get_var(av, tvar)) stop "Error getting variable info"
            if (tvar%var_type == V_STATE) then
               if (tvar%sheet) then
                  sv = sv + 1
                  call AED_InitCondition(self, self%mirror_state(:, n_vars + sv), tvar%name, tvar%initial)
               else
                  v = v + 1
                  call AED_InitCondition(self, self%mirror_state(:, v), tvar%name, tvar%initial)
               end if
            end if
         end do

         ! Push the (possibly overridden) initial state into the API arrays
         call mirror_to_api(self)

         write(*,"(/,5X,'----------  AED config : end  ----------',/)")
      end associate
   end subroutine init


   ! Called in the main loop of simstrat (in simstrat.f90) at each timestep.
   ! Updates the environment copies, runs one AED biogeochemistry step through
   ! libaed-api (mobility, light, kinetics, flux integration) and then applies
   ! vertical diffusion to the AED state variables.

   subroutine update(self, state)
      implicit none

      ! Arguments
      class(SimstratAED) :: self
      class(ModelState) :: state

      ! Local variables
      integer :: i, v
      logical :: do_surface

      associate (grid => self%grid, nz => self%grid%nz_grid)

         ! (1) Refresh per-step environment copies
         self%dt_c = state%dt
         self%yearday_c = day_of_year(state%current_year, state%current_month, state%current_day)
         self%air_temp_c = state%T_atm
         self%longwave_c = state%ha
         self%blueice_c = state%black_ice_h
         self%whiteice_c = state%white_ice_h
         self%colarea_c = grid%Az(grid%ubnd_fce)

         self%pres_c(1:nz) = -grid%z_volume(1:nz)
         self%heights_c(1) = grid%h(1)
         do i = 2, nz
            self%heights_c(i) = self%heights_c(i-1) + grid%h(i)
         end do

         ! In-water light field from Simstrat's radiation profile [W/m2]
         self%par_c(:) = state%rad_vol(:)*rho_0*cp
         self%rad_c(:) = self%par_c(:)

         ! (2) Push the Simstrat-side state (possibly modified by advection,
         !     inflow and diffusion since the last call) into the API arrays
         call mirror_to_api(self)

         ! (2b) Refresh zone environment (temperature, salinity, light, ...) at
         !      each zone's representative layer. Zone geometry and pelagic
         !      state exchange are handled by the registered zone callbacks;
         !      these physical fields aren't carried by that callback
         !      interface, so they're set directly here instead.
         if (self%n_zones > 0) call update_zone_environment(self, state)

         ! (3) Run one AED timestep (no surface exchange under ice cover)
         do_surface = .not. (state%total_ice_h > 0)
         call aed_run_model(1, grid%nz_occupied, do_surface)

         ! (4) Pull results back into the Simstrat-layout mirrors
         call api_to_mirror(self)
         if (self%n_zones > 0) call update_zone_output(self)

         ! (5) Light absorption feedback to the Simstrat temperature model
         if (self%aed_cfg%bioshade_feedback) then
            state%absorb_vol(1:grid%nz_occupied) = self%aed_cfg%background_extinction &
                                                   + self%bioext_c(1:grid%nz_occupied)
            call grid%interpolate_to_face(grid%z_volume, state%absorb_vol, grid%nz_occupied, state%absorb)
         end if

         ! (6) Diffusive transport of AED variables (advective transport is done
         !     in the usual Simstrat routines lateral/lateral_rho/advection)
         do v = 1, self%n_vars
            call diffusion_AED_state(self, state, v)
         end do
      end associate
   end subroutine update


   ! ---------------------------------------------------------------------------
   ! Layout conversion between the Simstrat mirrors and the API arrays
   ! ---------------------------------------------------------------------------

   subroutine mirror_to_api(self)
      class(SimstratAED) :: self
      integer :: v, sv

      do v = 1, self%n_vars
         self%cc(v, :) = self%mirror_state(:, v)
      end do
      do sv = 1, self%n_vars_ben
         ! Benthic (sheet) state lives in the bottom cell of the mirror column
         self%cc_hz(sv) = self%mirror_state(1, self%n_vars + sv)
      end do
   end subroutine mirror_to_api


   subroutine api_to_mirror(self)
      class(SimstratAED) :: self
      integer :: v, sv, d

      do v = 1, self%n_vars
         self%mirror_state(:, v) = self%cc(v, :)
      end do
      do sv = 1, self%n_vars_ben
         self%mirror_state(1, self%n_vars + sv) = self%cc_hz(sv)
      end do
      do d = 1, self%n_vars_diag
         self%mirror_diag(:, d) = self%cc_diag(d, :)
      end do
      self%mirror_diag_sheet(:) = self%cc_diag_hz(:)
   end subroutine api_to_mirror


   ! ---------------------------------------------------------------------------
   ! Memory allocation
   ! ---------------------------------------------------------------------------

   subroutine allocate_memory(self)
      implicit none

      ! Arguments
      class(SimstratAED) :: self

      ! Local variables
      integer :: status, nz

      nz = self%grid%nz_grid

      associate (n_vars => self%n_vars, n_vars_ben => self%n_vars_ben, &
                 n_vars_diag => self%n_vars_diag, n_vars_diag_sheet => self%n_vars_diag_sheet)

         ! API-layout arrays
         allocate (self%cc(n_vars, nz), stat=status)
         if (status /= 0) stop 'allocate_memory(): Error allocating (cc)'
         allocate (self%cc_hz(n_vars_ben), stat=status)
         if (status /= 0) stop 'allocate_memory(): Error allocating (cc_hz)'
         allocate (self%cc_diag(n_vars_diag, nz), stat=status)
         if (status /= 0) stop 'allocate_memory(): Error allocating (cc_diag)'
         allocate (self%cc_diag_hz(n_vars_diag_sheet), stat=status)
         if (status /= 0) stop 'allocate_memory(): Error allocating (cc_diag_hz)'
         self%cc = 0.0_RK
         self%cc_hz = 0.0_RK
         self%cc_diag = 0.0_RK
         self%cc_diag_hz = 0.0_RK

         ! Simstrat-layout mirrors
         allocate (self%mirror_state(nz, n_vars + n_vars_ben), stat=status)
         if (status /= 0) stop 'allocate_memory(): Error allocating (mirror_state)'
         allocate (self%mirror_diag(nz, n_vars_diag), stat=status)
         if (status /= 0) stop 'allocate_memory(): Error allocating (mirror_diag)'
         allocate (self%mirror_diag_sheet(n_vars_diag_sheet), stat=status)
         if (status /= 0) stop 'allocate_memory(): Error allocating (mirror_diag_sheet)'
         self%mirror_state = 0.0_RK
         self%mirror_diag = 0.0_RK
         self%mirror_diag_sheet = 0.0_RK

         ! Scalar environment targets
         allocate (self%dt_c, self%yearday_c, self%air_temp_c, self%air_pres_c, &
                   self%longwave_c, self%blueice_c, self%whiteice_c, self%colarea_c, &
                   self%sedzone_c, self%kw_c, self%lat_c, self%lon_c, self%bot_i, self%active_c)
         self%dt_c = 0.0_RK; self%yearday_c = 0.0_RK; self%air_temp_c = 0.0_RK
         self%air_pres_c = 0.0_RK; self%longwave_c = 0.0_RK
         self%blueice_c = 0.0_RK; self%whiteice_c = 0.0_RK; self%colarea_c = 0.0_RK
         self%sedzone_c = 0.0_RK; self%kw_c = 0.0_RK; self%lat_c = 0.0_RK; self%lon_c = 0.0_RK
         self%bot_i = 1
         self%active_c = .true.

         ! Array environment targets
         allocate (self%par_c(nz), self%rad_c(nz), self%pres_c(nz), self%tss_c(nz), &
                   self%nir_c(nz), self%uva_c(nz), self%uvb_c(nz), self%heights_c(nz), &
                   self%bioext_c(nz), self%biodrag_c(nz), self%sedzones_c(nz), self%zeros_c(nz))
         self%par_c = 0.0_RK; self%rad_c = 0.0_RK; self%pres_c = 0.0_RK; self%tss_c = 0.0_RK
         self%nir_c = 0.0_RK; self%uva_c = 0.0_RK; self%uvb_c = 0.0_RK; self%heights_c = 0.0_RK
         self%bioext_c = 0.0_RK; self%biodrag_c = 0.0_RK; self%sedzones_c = 0.0_RK
         self%zeros_c = 0.0_RK

         ! Names: names holds pelagic+benthic combined (matches AED2_state's
         ! n_vars+n_vars_ben columns); bennames aliases its benthic tail.
         allocate (self%names(n_vars + n_vars_ben), stat=status)
         if (status /= 0) stop 'allocate_memory(): Error allocating (names)'
         self%bennames => self%names(n_vars + 1:n_vars + n_vars_ben)
         allocate (self%diagnames(n_vars_diag), stat=status)
         if (status /= 0) stop 'allocate_memory(): Error allocating (diagnames)'
         allocate (self%diagnames_sheet(n_vars_diag_sheet), stat=status)
         if (status /= 0) stop 'allocate_memory(): Error allocating (diagnames_sheet)'
      end associate
   end subroutine allocate_memory


   ! ---------------------------------------------------------------------------
   ! Sediment zones (GLM-style depth bands, benthic_mode > 1). Registers the
   ! zone count/heights and geometry/state-exchange callbacks with
   ! libaed-api; zone-level benthic state (aedZones(:)%z_cc_hz) then evolves
   ! correctly inside aed_run_model. Zone environment (temperature, salinity,
   ! light, ...) is refreshed directly each step in `update`, since the
   ! API's zone callback interfaces don't carry those fields.
   ! ---------------------------------------------------------------------------

   subroutine init_zones(self, aed_cfg)
      implicit none

      class(SimstratAED) :: self
      class(AEDConfig) :: aed_cfg

      procedure(calc_zone_areas_t), pointer :: p_areas
      procedure(copy_to_zone_t), pointer :: p_to
      procedure(copy_from_zone_t), pointer :: p_from
      integer :: status, zon

      self%n_zones = aed_cfg%n_zones
      allocate (self%zone_heights_m(self%n_zones))
      self%zone_heights_m = aed_cfg%zone_heights

      allocate (self%zone_cc(self%n_vars, 1, self%n_zones), stat=status)
      if (status /= 0) stop 'init_zones(): Error allocating (zone_cc)'
      allocate (self%zone_cc_hz(self%n_vars_ben, self%n_zones + 1), stat=status)
      if (status /= 0) stop 'init_zones(): Error allocating (zone_cc_hz)'
      allocate (self%zone_diag(self%n_vars_diag, 1, self%n_zones), stat=status)
      if (status /= 0) stop 'init_zones(): Error allocating (zone_diag)'
      allocate (self%zone_diag_hz(self%n_vars_diag_sheet, self%n_zones + 1), stat=status)
      if (status /= 0) stop 'init_zones(): Error allocating (zone_diag_hz)'
      self%zone_cc = 0.0_RK
      self%zone_cc_hz = 0.0_RK
      self%zone_diag = 0.0_RK
      self%zone_diag_hz = 0.0_RK

      allocate (self%zone_state_out(self%n_zones, self%n_vars_ben), stat=status)
      if (status /= 0) stop 'init_zones(): Error allocating (zone_state_out)'
      allocate (self%zone_diag_sheet_out(self%n_zones, self%n_vars_diag_sheet), stat=status)
      if (status /= 0) stop 'init_zones(): Error allocating (zone_diag_sheet_out)'
      self%zone_state_out = 0.0_RK
      self%zone_diag_sheet_out = 0.0_RK

      call aed_init_zones(self%n_zones, 1, self%zone_cc, self%zone_cc_hz, self%zone_diag, self%zone_diag_hz)

      do zon = 1, self%n_zones
         aedZones(zon)%z_env%z_height = self%zone_heights_m(zon)
         aedZones(zon)%longitude => self%lon_c
         aedZones(zon)%latitude => self%lat_c
         ! zones are always fully wet here (no riparian/dry handling; that is benthic_mode 3)
         aedZones(zon)%z_env%z_pc_wet = 1.0_RK
      end do

      p_areas => simstrat_calc_zone_areas
      p_to => simstrat_copy_to_zone
      p_from => simstrat_copy_from_zone
      call api_set_zone_funcs(p_to, p_from, p_areas)
   end subroutine init_zones


   ! Copies each zone's own (kinetics-evolved) benthic state/diagnostics
   ! into the (n_zones, n_vars) output-facing arrays exposed via ModelState,
   ! for the dedicated <name>_zone_out.dat files. This is a transpose of
   ! zone_cc_hz/zone_diag_hz (variable-major) into zone-major layout, and is
   ! purely for output -- it does not feed back into the simulation.
   subroutine update_zone_output(self)
      implicit none
      class(SimstratAED) :: self
      integer :: zon

      do zon = 1, self%n_zones
         if (self%n_vars_ben > 0) self%zone_state_out(zon, :) = aedZones(zon)%z_cc_hz(:)
         if (self%n_vars_diag_sheet > 0) self%zone_diag_sheet_out(zon, :) = aedZones(zon)%z_cc_diag_hz(:)
      end do
   end subroutine update_zone_output


   ! Refreshes each zone's physical environment (temperature, salinity,
   ! light, meteorology, ...) from Simstrat's own state at that zone's
   ! representative layer, each timestep, ahead of the AED run.
   subroutine update_zone_environment(self, state)
      implicit none
      class(SimstratAED) :: self
      class(ModelState) :: state

      integer :: zon, lev

      do zon = 1, self%n_zones
         lev = zone_layer(self%zone_heights_m(zon), self%heights_c, self%grid%nz_occupied)
         aedZones(zon)%z_env%z_temp = state%T(lev)
         aedZones(zon)%z_env%z_salt = state%S(lev)
         aedZones(zon)%z_env%z_rho = state%rho(lev)
         aedZones(zon)%z_env%z_pres = self%pres_c(lev)
         aedZones(zon)%z_env%z_extc = state%absorb_vol(lev)
         aedZones(zon)%z_env%z_par = self%par_c(lev)
         aedZones(zon)%z_env%z_nir = self%nir_c(lev)
         aedZones(zon)%z_env%z_uva = self%uva_c(lev)
         aedZones(zon)%z_env%z_uvb = self%uvb_c(lev)
         aedZones(zon)%z_env%z_tss = self%tss_c(lev)
         aedZones(zon)%z_env%z_wind = state%uv10
         aedZones(zon)%z_env%z_air_temp = self%air_temp_c
         aedZones(zon)%z_env%z_air_pres = self%air_pres_c
         aedZones(zon)%z_env%z_rain = state%rain
         aedZones(zon)%z_env%z_I_0 = state%rad0
         aedZones(zon)%z_env%z_longwave = self%longwave_c
         aedZones(zon)%z_env%z_layer_stress = state%u_taub
         aedZones(zon)%z_env%z_sed_zone = real(zon, RK)
         aedZones(zon)%z_env%z_sed_zones = real(zon, RK)
      end do
   end subroutine update_zone_environment


   ! Representative layer for a zone: the shallowest layer whose cumulative
   ! height (from the bottom) first exceeds the zone's own upper boundary
   ! height, matching libaed-api's own internal zlev search in aed_api.F90.
   pure function zone_layer(zone_height, heights, wlev) result(lev)
      real(RK), intent(in) :: zone_height
      real(RK), dimension(:), intent(in) :: heights
      integer, intent(in) :: wlev
      integer :: lev, i

      lev = 0
      do i = 1, wlev
         if (zone_height < heights(i)) then
            lev = i
            exit
         end if
      end do
      if (lev == 0) lev = wlev
   end function zone_layer


   ! Zone footprint area/thickness from the host's per-layer geometry, at
   ! each zone's representative layer.
   subroutine simstrat_calc_zone_areas(theZones, n_zones, areas, heights, wlev)
      implicit none
      type(api_zone_t), dimension(:), intent(inout) :: theZones
      integer, intent(in) :: n_zones
      real(RK), dimension(:), pointer, intent(in) :: areas
      real(RK), dimension(:), pointer, intent(in) :: heights
      integer, intent(in) :: wlev

      integer :: zon, lev

      do zon = 1, n_zones
         lev = zone_layer(theZones(zon)%z_env%z_height, heights, wlev)
         theZones(zon)%z_env%z_area = areas(lev)
         if (lev > 1) then
            theZones(zon)%z_env%z_dz = heights(lev) - heights(lev - 1)
         else
            theZones(zon)%z_env%z_dz = heights(lev)
         end if
         theZones(zon)%z_env%z_col_area = areas(wlev)
         theZones(zon)%z_env%z_col_depth = heights(wlev)
         theZones(zon)%z_env%z_depth = heights(wlev) - heights(lev)
      end do
   end subroutine simstrat_calc_zone_areas


   ! Overlying water column pelagic state/diagnostics at each zone's
   ! representative layer, so zone benthic kinetics can read e.g. the
   ! dissolved oxygen concentration in the water immediately above it.
   subroutine simstrat_copy_to_zone(theZones, n_zones, heights, x_cc, x_cc_hz, x_diag, x_diag_hz, wlev)
      implicit none
      type(api_zone_t), dimension(:), intent(inout) :: theZones
      integer, intent(in) :: n_zones
      real(RK), dimension(:), pointer, intent(in) :: heights
      real(RK), dimension(:,:), pointer, intent(in) :: x_cc
      real(RK), dimension(:), pointer, intent(in) :: x_cc_hz
      real(RK), dimension(:,:), pointer, intent(in) :: x_diag
      real(RK), dimension(:), pointer, intent(in) :: x_diag_hz
      integer, intent(in) :: wlev

      integer :: zon, lev

      do zon = 1, n_zones
         lev = zone_layer(theZones(zon)%z_env%z_height, heights, wlev)
         if (size(x_cc, 1) > 0) theZones(zon)%z_cc(:, zon) = x_cc(:, lev)
         if (size(x_diag, 1) > 0) theZones(zon)%z_cc_diag(:, zon) = x_diag(:, lev)
      end do
   end subroutine simstrat_copy_to_zone


   ! Reports each zone's (now kinetics-evolved) benthic state back onto
   ! Simstrat's single-scalar-per-variable cc_hz/diag_hz slots, as an
   ! area-weighted average across zones. Simstrat's own state arrays are not
   ! zone-resolved, so this is a summary for output/mass-balance purposes;
   ! each zone's own mass (theZones(:)%z_cc_hz) is unaffected by this and
   ! keeps evolving correctly regardless. Pelagic arrays (x_cc/x_diag) are
   ! deliberately left untouched here: the water-column feedback from
   ! benthic fluxes is already applied via the flux-based disaggregation
   ! built into libaed-api's GLM-style zone handling -- touching x_cc here
   ! too would double-count it.
   subroutine simstrat_copy_from_zone(theZones, n_zones, heights, x_cc, x_cc_hz, x_diag, x_diag_hz, wlev)
      implicit none
      type(api_zone_t), dimension(:), intent(in) :: theZones
      integer, intent(in) :: n_zones
      real(RK), dimension(:), pointer, intent(in) :: heights
      real(RK), dimension(:,:), pointer, intent(inout) :: x_cc
      real(RK), dimension(:), pointer, intent(inout) :: x_cc_hz
      real(RK), dimension(:,:), pointer, intent(inout) :: x_diag
      real(RK), dimension(:), pointer, intent(inout) :: x_diag_hz
      integer, intent(in) :: wlev

      integer :: zon
      real(RK) :: total_area, w

      total_area = 0.0_RK
      do zon = 1, n_zones
         total_area = total_area + theZones(zon)%z_env%z_area
      end do
      if (total_area <= 0.0_RK) return

      if (size(x_cc_hz) > 0) x_cc_hz = 0.0_RK
      if (size(x_diag_hz) > 0) x_diag_hz = 0.0_RK
      do zon = 1, n_zones
         w = theZones(zon)%z_env%z_area/total_area
         if (size(x_cc_hz) > 0) x_cc_hz = x_cc_hz + w*theZones(zon)%z_cc_hz
         if (size(x_diag_hz) > 0) x_diag_hz = x_diag_hz + w*theZones(zon)%z_cc_diag_hz
      end do
   end subroutine simstrat_copy_from_zone


   ! ---------------------------------------------------------------------------
   ! Variable names, printed to screen and used by the output module
   ! ---------------------------------------------------------------------------

   subroutine assign_var_names(self)
      implicit none

      ! Arguments
      class(SimstratAED) :: self

      ! Local variables
      type(aed_variable_t), pointer :: tvar
      integer :: i, j

      print "(5X,'Configured variables to simulate:')"

      j = 0
      do i = 1, self%n_aed_vars
         if (aed_get_var(i, tvar)) then
            if (.not. tvar%sheet .and. tvar%var_type == V_STATE) then
               j = j + 1
               self%names(j) = trim(tvar%name)
               print *, "     S(", j, ") AED pelagic(3D) variable: ", trim(self%names(j))
            end if
         end if
      end do

      j = 0
      do i = 1, self%n_aed_vars
         if (aed_get_var(i, tvar)) then
            if (tvar%sheet .and. tvar%var_type == V_STATE) then
               j = j + 1
               self%bennames(j) = trim(tvar%name)
               print *, "     B(", j, ") AED benthic(2D) variable: ", trim(self%bennames(j))
            end if
         end if
      end do

      j = 0
      do i = 1, self%n_aed_vars
         if (aed_get_var(i, tvar)) then
            if (.not. tvar%sheet .and. tvar%var_type == V_DIAGNOSTIC) then
               j = j + 1
               self%diagnames(j) = trim(tvar%name)
               print *, "     D(", j, ") AED diagnostic 3D variable: ", trim(tvar%name)
            end if
         end if
      end do

      j = 0
      do i = 1, self%n_aed_vars
         if (aed_get_var(i, tvar)) then
            if (tvar%sheet .and. tvar%var_type == V_DIAGNOSTIC) then
               j = j + 1
               self%diagnames_sheet(j) = trim(tvar%name)
               print *, "     D(", j, ") AED diagnostic 2D variable: ", trim(tvar%name)
            end if
         end if
      end do
   end subroutine assign_var_names


   ! ---------------------------------------------------------------------------
   ! Initial conditions from <var>_ini.dat files (same convention as before)
   ! ---------------------------------------------------------------------------

   subroutine AED_InitCondition(self, var, varname, default_val)
      implicit none

      class(SimstratAED) :: self
      real(RK), intent(inout) :: var(1:self%grid%nz_grid) ! Vector of initial conditions
      real(RK), intent(in) :: default_val ! Depth-independent value (default from the AED namelist)
      character(len=*), intent(in) :: varname ! Identifying the variable

      real(RK) :: z_read(self%grid%max_length_input_data), var_read(self%grid%max_length_input_data)
      character(len=100) :: fname
      integer :: i, nval

      fname = trim(self%aed_cfg%path_aed_initial)//trim(varname)//'_ini.dat'
      open (14, action='read', status='unknown', err=1, file=fname) ! Opens initial conditions file
      write (6, *) 'reading initial conditions of ', trim(varname)
      read (14, *) ! Skip header
      do i = 1, self%grid%max_length_input_data ! Read initial values
         read (14, *, end=9) z_read(i), var_read(i)
      end do
9     nval = i - 1 ! Number of values
      if (nval < 0) then
         write (6, *) 'Error reading ', trim(varname), ' initial conditions file (no data found).'
         stop
      end if
      close (14)
      do i = 1, nval
         z_read(i) = abs(z_read(i)) ! Make depths positive
      end do

      call reverse_in_place(z_read(1:nval))
      z_read(1:nval) = self%grid%z_zero - z_read(1:nval)
      call reverse_in_place(var_read(1:nval))

      if (nval == 1) then
         write (6, *) '      Only one row! Water column will be initially homogeneous.'
         var(1:self%grid%nz_grid) = var_read(1)
      else
         call Interp(z_read(1:nval), var_read(1:nval), nval, self%grid%z_volume, var, self%grid%nz_grid)
      end if
      return

1     write (6, *) '   File ''', trim(fname), ''' not found. Initial conditions set to default value from the AED namelist.'
      var(1:self%grid%nz_grid) = default_val
      return
   end subroutine AED_InitCondition


   ! ---------------------------------------------------------------------------
   ! Vertical diffusion of AED state variables. Copy of the diffusion algorithm
   ! used for the Simstrat state variables (see strat_discretization).
   ! ---------------------------------------------------------------------------

   subroutine diffusion_AED_state(self, state, var_index)
      implicit none

      ! Arguments
      class(SimstratAED) :: self
      class(ModelState) :: state
      integer :: var_index

      ! Local variables
      real(RK), dimension(self%grid%ubnd_vol) :: boundaries, sources, lower_diag, main_diag, upper_diag, rhs

      boundaries = 0.
      sources = 0.

      if (var_index == state%n_pH) state%AED2_state(:, state%n_pH) = 10.**(-state%AED2_state(:, state%n_pH))
      call euleri_create_LES_MFQ_AED(self, state%AED2_state(:, var_index), state%nuh, sources, boundaries, lower_diag, main_diag, upper_diag, rhs, state%dt)
      call solve_tridiag_thomas(lower_diag, main_diag, upper_diag, rhs, state%AED2_state(:, var_index), self%grid%ubnd_vol)
      if (var_index == state%n_pH) state%AED2_state(:, state%n_pH) = -log10(state%AED2_state(:, state%n_pH))
   end subroutine diffusion_AED_state


   subroutine euleri_create_LES_MFQ_AED(self, var, nu, sources, boundaries, lower_diag, main_diag, upper_diag, rhs, dt)
      class(SimstratAED), intent(inout) :: self
      real(RK), dimension(:), intent(inout) :: var, sources, boundaries, lower_diag, upper_diag, main_diag, rhs, nu
      real(RK), intent(inout) :: dt
      integer :: n

      n = self%grid%ubnd_vol

      ! Build diagonals
      upper_diag(1) = 0.0_RK
      upper_diag(2:n) = dt*nu(2:n)*self%grid%AreaFactor_1(2:n)
      lower_diag(1:n - 1) = dt*nu(2:n)*self%grid%AreaFactor_2(1:n - 1)
      lower_diag(n) = 0.0_RK
      main_diag(1:n) = 1.0_RK - upper_diag(1:n) - lower_diag(1:n) + boundaries(1:n)*dt

      ! Calculate RHS
      ! A*phi^{n+1} = phi^{n}+dt*S^{n}
      rhs(1:n) = var(1:n) + dt*sources(1:n)
   end subroutine euleri_create_LES_MFQ_AED


   ! ---------------------------------------------------------------------------
   ! Mobility (settling/rising) callback registered with libaed-api. The API
   ! calls this once per mobile state variable and timestep. The algorithm is
   ! taken from the GLM code (http://aed.see.uwa.edu.au/research/models/GLM/),
   ! as in the old coupler, but works purely on the arguments supplied by the
   ! API (layer thickness h, layer areas A) so it matches aed_mobility_fn_t.
   !
   ! Assumptions:
   ! 1) movement direction has at most one change down the layers
   ! 2) sides of the lake slope inward (i.e. bottom is narrower than top)
   ! ---------------------------------------------------------------------------

   subroutine mobility(N, dt, h, A, ww, min_C, cc)
      implicit none

      ! Arguments (must match aed_mobility_fn_t exactly)
      integer(kind=C_INT32_T), intent(in) :: N     ! number of vertical layers
      real(kind=C_DOUBLE), intent(in) :: dt        ! time step (s)
      real(kind=C_DOUBLE), intent(in) :: h(:)      ! layer thickness (m)
      real(kind=C_DOUBLE), intent(in) :: A(:)      ! layer areas (m2)
      real(kind=C_DOUBLE), intent(in) :: ww(:)     ! vertical speed (m/s)
      real(kind=C_DOUBLE), intent(in) :: min_C     ! minimum allowed cell concentration
      real(kind=C_DOUBLE), intent(inout) :: cc(:)  ! cell concentration

      ! Local variables
      real(RK) :: dtMax, tdt, tmp
      real(RK), dimension(N) :: mins, vols, Y
      integer :: dirChng, signum, i, count

      if (N < 1 .or. dt <= 0.0_RK) return

      ! Determine mobility timestep, i.e. the maximum time step for which
      ! particles will not pass through more than one layer
      dtMax = dt
      dirChng = 0 ! layer at which direction switches from sinking to rising or vice versa
      signum = int(sign(1.0_RK, ww(1))) ! positive for rising, negative for sinking

      do i = 1, N
         ! for convenience
         vols(i) = A(i)*h(i)
         mins(i) = min_C*vols(i)
         Y(i) = cc(i)*vols(i)

         ! look for the change of direction
         if (signum /= int(sign(1.0_RK, ww(i)))) then
            signum = -signum
            dirChng = i - 1
         end if

         ! check if all movement can be from within this cell
         if (abs(ww(i)*dt) > h(i)) then
            tdt = h(i)/abs(ww(i))
            if (tdt < dtMax) dtMax = tdt
         end if

         ! check if movement can all be into the next cell.
         ! if movement is settling, next is below, otherwise next is above.
         if (ww(i) > 0.) then
            if ((i < N) .and. (abs(ww(i))*dt) > h(i + 1)) then
               tdt = h(i + 1)/abs(ww(i))
               if (tdt < dtMax) dtMax = tdt
            end if
         else if ((i > 1) .and. (abs(ww(i))*dt) > h(i - 1)) then
            tdt = h(i - 1)/abs(ww(i))
            if (tdt < dtMax) dtMax = tdt
         end if
      end do
      if (dirChng == 0 .and. ww(1) > 0.) dirChng = N ! all rising
      if (dirChng == 0 .and. ww(N) < 0.) dirChng = N ! all sinking

      ! Repeat in steps of dtMax until the full timestep dt is consumed
      tdt = dtMax
      count = 0
      do
         count = count + 1
         if (count*dtMax > dt) tdt = dt - (count - 1)*dtMax ! last (partial) substep

         ! 2 possibilities:
         ! 1) lower levels rising, upper levels sinking
         ! 2) lower levels sinking, upper levels rising
         if (ww(1) > 0.) then ! lower levels rising
            if (ww(N) < 0.) then ! top levels are sinking
               call Sinking(Y, cc, ww, vols, mins, A, tdt, N, dirChng, tmp)
               Y(dirChng) = Y(dirChng) + tmp
               cc(dirChng) = Y(dirChng)/vols(dirChng)
            end if
            call Rising(Y, cc, ww, vols, mins, A, tdt, 1, dirChng)
         else ! lower levels sinking
            call Sinking(Y, cc, ww, vols, mins, A, tdt, dirChng, 1, tmp)
            if (ww(N) > 0.) then ! top levels are rising
               call Rising(Y, cc, ww, vols, mins, A, tdt, dirChng, N)
            end if
         end if

         if (count*dtMax >= dt) exit
      end do
   end subroutine mobility


   ! Rising is the easier of the two since the slope means we don't need to look
   ! at relative areas (the cell above will always be >= the current cell);
   ! all moving matter is moved to the next cell up.

   subroutine Rising(Y, conc, settling_v, vols, mins, A, dt, start_i, end_i)
      ! Arguments
      real(RK), dimension(:), intent(inout) :: conc, Y
      real(RK), dimension(:), intent(in) :: settling_v, vols, mins, A
      real(RK) :: dt
      integer :: start_i, end_i

      ! Local variables
      real(RK) :: mov, moved
      integer :: i

      mov = 0.
      moved = 0.

      do i = start_i, end_i
         ! speed times time (=h) times area times concentration = mass to move
         mov = (settling_v(i)*dt)*A(i)*conc(i)
         ! if removing that much would bring it below min conc
         if ((Y(i) + moved - mov) < mins(i)) mov = Y(i) + moved - mins(i)

         Y(i) = Y(i) + moved - mov
         conc(i) = Y(i)/vols(i) ! return it to a concentration
         moved = mov ! for the next step
      end do
      ! nothing rises out of the top cell, but we still add that which came from below
      Y(end_i) = Y(end_i) + moved
      conc(end_i) = Y(end_i)/vols(end_i)
   end subroutine Rising


   ! For each cell: calculate how much is going to move, subtract it, add what
   ! moved from the previous (upper) cell and fix the concentration; part of the
   ! moving mass settles on the flanks (area ratio of the cell below).

   subroutine Sinking(Y, conc, settling_v, vols, mins, A, dt, start_i, end_i, moved)
      ! Arguments
      real(RK), dimension(:), intent(inout) :: conc, Y
      real(RK), dimension(:), intent(in) :: settling_v, vols, mins, A
      real(RK) :: dt
      real(RK), intent(out) :: moved
      integer :: start_i, end_i

      ! Local variables
      real(RK) :: mov
      integer :: i

      mov = 0.
      moved = 0.

      do i = start_i, end_i, -1
         ! speed times time (=h) times area times concentration = mass to move
         mov = (abs(settling_v(i))*dt)*A(i)*conc(i)
         ! if removing that much would bring it below min conc
         if ((Y(i) + moved - mov) < mins(i)) mov = Y(i) + moved - mins(i)
         Y(i) = Y(i) + moved - mov
         conc(i) = Y(i)/vols(i) ! return it to a concentration

         ! now mov holds how much has moved out of the cell, but not all
         ! of that will go into the next cell (part of it settles on the flanks)
         if (i > 1) then
            moved = mov*(A(i - 1)/A(i)) ! for the next step
         else
            moved = mov ! we are about to exit anyway
         end if
      end do
   end subroutine Sinking


   ! ---------------------------------------------------------------------------
   ! Day of year (including day fraction), leap-year aware
   ! ---------------------------------------------------------------------------

   pure function day_of_year(year, month, day) result(yd)
      integer, intent(in) :: year, month
      real(RK), intent(in) :: day
      real(RK) :: yd

      integer, parameter :: days_before(12) = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]

      yd = days_before(month) + day
      if (month > 2) then
         if (mod(year, 4) == 0 .and. (mod(year, 100) /= 0 .or. mod(year, 400) == 0)) yd = yd + 1.0_RK
      end if
   end function day_of_year

end module simstrat_aed
