module tol_reflectance
    use iso_c_binding
    use start_mod ! contains PI
    implicit none
#include <flotsam.inc>
    ! edge pressure for single layer case
    real(FLOTSAM_REAL), dimension(2) :: tol_edge_pressure = (/0, 101300/)
    ! IMPORTANT: this index vector is intended for C++, must be 0-based!
    ! in theory depends on number of atm. layers, which is 1 in case of TOL
    integer(c_int), dimension(1) :: loc = (/0/)

    real(FLOTSAM_REAL), parameter :: pi_flotsam=4.D0*ATAN(1.D0)

contains

subroutine init_flotsam_channel()
    !Real, Intent(In) :: sa, tau_ap
    !integer, intent(out) :: iband
    integer :: ichan
    !Real, Dimension(1:N_AOD_max), Intent(In) :: AOD_max_steps   ! Aerosol parameters
    !real, intent(in) :: tau_ap ! tau_ap is tau_force in start_exe
    !real :: tau_apriori

    !interface
        ! flotsam reflectance interface
   !end interface

    print*, 'inside LM_flotsam (LM)'

    ! Create channel (wavelength-dependent information)
    ! -------------------------------------------------------------------
    ! use ..._vacuum if rayleigh scattering is OFF
    ichan = flotsam_new_channel_vacuum()
    !print*, 'inside LM_flotsam (LM)'

    ! What is AOD_max_steps exactly?
    ! -> seems to be just the range of AOD values.. see aod_parameters.f90
    !print*, AOD_MAX_STEPS

    !tau_apriori = tau_ap
!   ! retrieving position from AOD steps of dynamic properties
    !kk = MINLOC(ABS(AOD_max_steps-tau_apriori), DIM=1)
    !print*, kk

end subroutine init_flotsam_channel

subroutine prepare_ref_tol_flotsam(&
    ichan, & !phase_function, &
    mu_sun, mu_view, saa, vaa, &
    sca_aerus, &
    !pf_interp, &
    iband, iprof, &
    n_sca, &
    n_pfc, n_layers, n_aero_layers)

   ! NOTE regarding allocatable vs automatic arrays
   !  automatic arrays have shorter syntax and are placed on the stack
   !    on some systems stack size can be quite limited (few kB)
   !    rule of thumb: use for small arrays
   !  allocatable arrays are placed on the heap -> use for larger arrays
   !  the placement can be adjusted by compile options.
   !  placement may be compiler dependent
   
   integer, intent(in) :: ichan, n_sca, n_pfc, n_layers, n_aero_layers
   !real, intent(in) :: phase_function(n_sca,1)
   real, intent(in) ::  saa, vaa,  mu_sun, mu_view
   real, intent(out) :: sca_aerus
   integer(c_int), intent(out) :: iband, iprof
   !real, intent(inout) ::  pf_interp
   !real(FLOTSAM_REAL), allocatable :: phase_function(:,:)
   real(FLOTSAM_REAL) :: azim, sca, mu_sun_flotsam, mu_view_flotsam
   !real(FLOTSAM_REAL) :: pfc_array(n_aero_layers, n_pfc), phase_function_flotsam(n_sca,1)
   integer(c_int) :: rstatus
   integer :: n_pf=1 ! check: what is this exactly?
   
   ! i guess this changes if no lambertian?
   ! --> indeed, might have to take ocean pixels into account here?
   !   in that case should use different refl. model, see ECRAD
   ! need to ask Xavier if we are doing that too, or just use lambertian here
   integer(c_int) :: n_albedo_components = 1
   
   !allocate(phase_function(n_sca, n_pf))
   mu_sun_flotsam = mu_sun
   mu_view_flotsam = mu_view
   !phase_function_flotsam = phase_function
   
   ! NOTE ichan coming from above is 1-based, but should be 0-based?
   !  -> pass I-1 in start_exe...?
   !print*, 'calc ref tol for ', ichan
   ! TODO: steps towards _init_band_profile
   !  is aod_pos in aerus-geo the same as aod_index in my codes?
   iband = flotsam_new_band_profile()
   !print*, iband
   
   iprof = flotsam_new_background_profile()
   
   ! if several aerosol layers, dim-2 size might be >1...?
   !  (should be n_pf instead of 1 in that case...?)
   ! phase_function = reshape(pf_smooth, (/n_sca,1/))
   ! can these 2 lines not be done in 1 step...?
   !pf_components = reshape(pf_components, (/n_pfc,1/))
   !pfc_array(1, :) = pf_components(:,1)
   
   rstatus = flotsam_set_edge_pressure(iprof, n_layers, tol_edge_pressure)
   rstatus = flotsam_init_band_profile(iband, ichan, iprof)
   
   azim = vaa - saa
   sca = flotsam_scattering_angle(&
       real(mu_sun, FLOTSAM_REAL), real(mu_view, FLOTSAM_REAL), azim)
   print*, sca
   sca_aerus = sca
   
   rstatus = flotsam_set_geometry(iband, mu_sun_flotsam, mu_view_flotsam, azim)
   
   !pf_interp = flotsam_interp_phase_func(&
   !    n_sca, real(phase_function, FLOTSAM_REAL), sca)
   !print*, pf_interp

end subroutine prepare_ref_tol_flotsam

subroutine calculate_ref_tol_flotsam_jacobian(&
    ichan, pf_interp, pf_components, ssa_val, &
    ref_s, aod, &
    mu_sun, mu_view, saa, vaa, &
    tol_reflectance, &
    n_sca,&
    n_pfc, n_layers, n_aero_layers)

! NOTE regarding allocatable vs automatic arrays
!  automatic arrays have shorter syntax and are place on the stack
!    on some systems stack size can be aquite limited (few kB)
!    rule of thumb: use for small arrays
!  allocatable arrays are placed on the heap -> use for larger arrays
!  the placement can be adjusted by compile options.
!  placement may be compiler dependent

integer, intent(in) :: ichan, n_sca, n_pfc, n_layers, n_aero_layers
real, intent(in) :: pf_interp, pf_components(n_pfc,1)
real, intent(in) :: ssa_val, ref_s, saa, vaa, aod, mu_sun, mu_view
real, intent(inout) :: tol_reflectance
real(FLOTSAM_REAL) ::  d_albedo, d_aod, d_ssa, d_pf, d_pfc
!real(FLOTSAM_REAL), allocatable :: phase_function(:,:)
real(FLOTSAM_REAL) :: radiance, azim, sca, mu_sun_flotsam, mu_view_flotsam
real(FLOTSAM_REAL) :: pfc_array(n_aero_layers, n_pfc), phase_function_flotsam(n_sca,1)
real(FLOTSAM_REAL), dimension(n_aero_layers) :: pf_interp_array, aod_array, ssa_array
real(FLOTSAM_REAL), dimension(n_aero_layers) :: ref_s_array
integer :: iband, iprof, rstatus
integer :: n_pf=1 ! check: what is this exactly?

! i guess this changes if no lambertian?
! --> indeed, might have to take ocean pixels into account here?
!   in that case should use different refl. model, see ECRAD
! need to ask Xavier if we are doing that too, or just use lambertian here
integer(c_int) :: n_albedo_components = 1

! NOTE ichan coming from above is 1-based, but should be 0-based?
!  -> pass I-1 in start_exe...?
!print*, 'calc ref tol for ', ichan
! TODO: steps towards _init_band_profile
!  is aod_pos in aerus-geo the same as aod_index in my codes?
!iband = flotsam_new_band_profile()
!print*, iband

!iprof = flotsam_new_background_profile()

! if several aerosol layers, dim-2 size might be >1...?
!  (should be n_pf instead of 1 in that case...?)
! phase_function = reshape(pf_smooth, (/n_sca,1/))
! can these 2 lines not be done in 1 step...?
!pf_components = reshape(pf_components, (/n_pfc,1/))
pfc_array(1, :) = pf_components(:,1)

!rstatus = flotsam_set_edge_pressure(iprof, n_layers, tol_edge_pressure)

! when do I have to re-init....?
!rstatus = flotsam_init_band_profile(iband, ichan, iprof)

ref_s_array(1) = ref_s
aod_array(1) = aod
ssa_array(1) = ssa_val
pf_interp_array(1) = pf_interp

!print*, n_sca, ssa_array, aod_array, loc, ref_s_array, mu_sun_flotsam, azim
!print*, '++++++++'
!print*, phase_function
!print*, pf_components
!print*, real(pf_components, FLOTSAM_REAL)
!print*, '=========='

radiance = 0.

rstatus = flotsam_reflectance_jacobian(iband, n_albedo_components, &
    ref_s_array, n_layers, loc, aod_array, ssa_array,&
    pf_interp_array, real(pf_components, FLOTSAM_REAL), radiance, d_albedo, d_aod, d_ssa, d_pf, d_pfc)

tol_reflectance = radiance * pi_flotsam / mu_sun
!tol_reflectance = 0.1

! do this after refl. is calculated
rstatus = flotsam_free_band_profile(iband)
rstatus = flotsam_free_background_profile(iprof)

end subroutine calculate_ref_tol_flotsam_jacobian

subroutine calculate_ref_tol_flotsam(&
                ichan, phase_function, pf_components, ssa_val, &
                ref_s, aod, &
                mu_sun, mu_view, saa, vaa, &
                tol_reflectance, &
                n_sca,&
                n_pfc, n_layers, n_aero_layers)

    ! NOTE regarding allocatable vs automatic arrays
    !  automatic arrays have shorter syntax and are place on the stack
    !    on some systems stack size can be aquite limited (few kB)
    !    rule of thumb: use for small arrays
    !  allocatable arrays are placed on the heap -> use for larger arrays
    !  the placement can be adjusted by compile options.
    !  placement may be compiler dependent

    integer, intent(in) :: ichan, n_sca, n_pfc, n_layers, n_aero_layers
    real, intent(in) :: phase_function(n_sca,1), pf_components(n_pfc,1)
    real, intent(in) :: ssa_val, ref_s, saa, vaa, aod, mu_sun, mu_view
    real, intent(inout) :: tol_reflectance
    !real(FLOTSAM_REAL), allocatable :: phase_function(:,:)
    real(FLOTSAM_REAL) :: radiance, azim, sca, mu_sun_flotsam, mu_view_flotsam
    real(FLOTSAM_REAL) :: pfc_array(n_aero_layers, n_pfc), pf_interp, phase_function_flotsam(n_sca,1)
    real(FLOTSAM_REAL), dimension(n_aero_layers) :: pf_interp_array, aod_array, ssa_array
    real(FLOTSAM_REAL), dimension(n_aero_layers) :: ref_s_array
    integer :: iband, iprof, rstatus
    integer :: n_pf=1 ! check: what is this exactly?

    ! i guess this changes if no lambertian?
    ! --> indeed, might have to take ocean pixels into account here?
    !   in that case should use different refl. model, see ECRAD
    ! need to ask Xavier if we are doing that too, or just use lambertian here
    integer(c_int) :: n_albedo_components = 1

    !allocate(phase_function(n_sca, n_pf))
    mu_sun_flotsam = mu_sun
    mu_view_flotsam = mu_view
    phase_function_flotsam = phase_function

    ! NOTE ichan coming from above is 1-based, but should be 0-based?
    !  -> pass I-1 in start_exe...?
    !print*, 'calc ref tol for ', ichan
    ! TODO: steps towards _init_band_profile
    !  is aod_pos in aerus-geo the same as aod_index in my codes?
    iband = flotsam_new_band_profile()
    !print*, iband

    iprof = flotsam_new_background_profile()

    ! if several aerosol layers, dim-2 size might be >1...?
    !  (should be n_pf instead of 1 in that case...?)
    ! phase_function = reshape(pf_smooth, (/n_sca,1/))
    ! can these 2 lines not be done in 1 step...?
    !pf_components = reshape(pf_components, (/n_pfc,1/))
    pfc_array(1, :) = pf_components(:,1)

    rstatus = flotsam_set_edge_pressure(iprof, n_layers, tol_edge_pressure)
    rstatus = flotsam_init_band_profile(iband, ichan, iprof)

    azim = vaa - saa
    sca = flotsam_scattering_angle(&
            real(mu_sun, FLOTSAM_REAL), real(mu_view, FLOTSAM_REAL), azim)
    print*, sca
    rstatus = flotsam_set_geometry(iband, mu_sun_flotsam, mu_view_flotsam, azim)

    pf_interp = flotsam_interp_phase_func(&
            n_sca, real(phase_function, FLOTSAM_REAL), sca)

    ref_s_array(1) = ref_s
    aod_array(1) = aod
    ssa_array(1) = ssa_val
    pf_interp_array(1) = pf_interp

    !print*, n_sca, ssa_array, aod_array, loc, ref_s_array, mu_sun_flotsam, azim
    !print*, '++++++++'
    !print*, phase_function
    !print*, pf_components
    !print*, real(pf_components, FLOTSAM_REAL)
    !print*, '=========='

    radiance = 0.

    rstatus = flotsam_reflectance(iband, n_albedo_components, &
        ref_s_array, n_layers, loc, aod_array, ssa_array,&
        pf_interp_array, real(pf_components, FLOTSAM_REAL), radiance) !, d_albedo, d_aod, d_ssa, d_pf, d_pfc)

    tol_reflectance = radiance * pi_flotsam / mu_sun
    !tol_reflectance = 0.1

    ! do this after refl. is calculated
    rstatus = flotsam_free_band_profile(iband)
    rstatus = flotsam_free_background_profile(iprof)

end subroutine calculate_ref_tol_flotsam


Subroutine calculate_ref_tol_msa(ref_s,p_hg_tilde,ssa_tilde,g_tilde,tau0_tilde,albed_b_in,uv,us,ref_tol)
  Real, Intent(In)     :: ref_s,p_hg_tilde,ssa_tilde,g_tilde,tau0_tilde,albed_b_in,uv,us
  Real, Intent(InOut)    :: ref_tol

  Real :: TS, TV, Salb, trans_coeff, rho_1, R_SS, x1_tilde, R_MS_us, R_MS_uv, R_MS

    x1_tilde = 3.*g_tilde
    ! Aerosol transmittances from Katsev et al.
    TS = exp(-tau0_tilde*(1.-ssa_tilde*(1.-(1.-g_tilde)/2.))/us)
    TV = exp(-tau0_tilde*(1.-ssa_tilde*(1.-(1.-g_tilde)/2.))/uv)
!    ! Spherical albedo from Wiscombe et al.
!    Salb = 2.*ssa_tilde*beta*tau0_tilde
    ! Spherical albedo from Katsev et al.
    Salb = tau0_tilde/(tau0_tilde+4./(3.-x1_tilde))
    trans_coeff = TS*TV/(1.-Salb*albed_b_in)
    ! Aerosol reflectance from Katsev et al.
    !     Single scattering
    rho_1 = 1./(4.*(us+uv))*(1.-exp(-tau0_tilde*(1./us+1./uv)))
    R_SS = ssa_tilde*p_hg_tilde*rho_1
    !     Multiple scattering
    R_MS_us = 1. + 1.5*us + (1.-1.5*us)*exp(-tau0_tilde/us)
    R_MS_uv = 1. + 1.5*uv + (1.-1.5*uv)*exp(-tau0_tilde/uv)
    R_MS = 1. - R_MS_us*R_MS_uv/(4.+(3.-x1_tilde)*tau0_tilde) + ((3.+x1_tilde)*us*uv-2.*(us+uv))*rho_1

    ! Calculating reconstruction error
    ref_tol = trans_coeff * ref_s + R_SS + R_MS

!    print*, 'trans_coeff*ref_s, R_SS+R_MS: ', trans_coeff*ref_s, R_SS+R_MS

End Subroutine calculate_ref_tol_msa


subroutine prepare_phase_functions_flotsam(scattering_angles, pf_lut, n_pfc, &
    & pf_smooth, pf_components, n_sca, n_aod, n_chan) bind(c, name='prepare_phase_functions_flotsam')
    integer(c_int), value, intent(in) :: n_sca, n_aod, n_chan, n_pfc
    real(FLOTSAM_REAL), intent(in) :: scattering_angles(n_sca)
    real, intent(in) :: pf_lut(n_chan, n_sca, n_aod)
    real(FLOTSAM_REAL), allocatable :: pf(:,:)
    real, intent(inout) :: pf_smooth(n_chan, n_sca, n_aod)
    real, intent(inout) :: pf_components(n_chan, n_pfc, n_aod)
    real(FLOTSAM_REAL) :: tmp_pf_smooth(n_sca, n_aod)
    real(FLOTSAM_REAL) :: tmp_pf_components(n_pfc, n_aod)
    integer(c_int) :: n_pf
    integer :: i_chan, i_aod, return_status

    ! execute flotsam_analyse_phase_functions for all channels

    n_pf = n_aod

    !allocate(pf_smooth(n_aero, n_sca, n_pf))
    !allocate(pf_components(n_aero, flotsam_n_phase_function_components(), n_pf))
    allocate(pf(n_sca, n_pf))

    do i_chan = 1, n_chan
            !print*, i_chan
            pf = pf_lut(i_chan,:,:)

            return_status = flotsam_analyse_phase_functions(&
                & n_pf, n_sca, scattering_angles, pf, 4.*pi_flotsam, &
                !& pf_smooth(i_aero, :,:), pf_components(i_aero, :,:) )
                & tmp_pf_smooth, tmp_pf_components )

            pf_components(i_chan,:,:) = tmp_pf_components
            pf_smooth(i_chan,:,:) = tmp_pf_smooth
    end do


end subroutine prepare_phase_functions_flotsam

end module tol_reflectance
