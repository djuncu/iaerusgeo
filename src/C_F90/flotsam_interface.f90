module flotsam_interface
    use iso_c_binding
    use py_ifc
    implicit none
#include <flotsam.inc>
    integer, parameter :: max_iter1=8, max_iter2=8 ! similar results for 3 iterations for all models except dust
    real, parameter :: alpha_lm=2., d_aod=0.01 
    real, parameter :: pi_lm=4.D0*ATAN(1.D0)
    real(FLOTSAM_REAL), parameter :: pi_flotsam=4.D0*ATAN(1.D0)
    ! i guess this changes if no lambertian?
    ! --> indeed, might have to take ocean pixels into account here?
    !   in that case should use different refl. model, see ECRAD
    ! need to ask Xavier if we are doing that too, or just use lambertian here
    integer(c_int) :: n_albedo_components=1, n_layers=1
    ! edge pressure for single layer case
    ! i guess pressure at TOL should not be 0 though?
    ! i think it could be around ~90200 assuming 1km altitude (mcclatchey midsum)
    real(FLOTSAM_REAL), dimension(2) :: tol_edge_pressure = (/90000, 101300/)
    ! IMPORTANT: this index vector is intended for C++, must be 0-based!
    ! in theory depends on number of atm. layers, which is 1 in case of TOL
    integer(c_int), dimension(1) :: loc = (/0/)

    contains

subroutine LM_flotsam(&
    Sa, ssa_pix, uv,us, &
    ! including flotsam preparation for testing
    saa, vaa, sca, &
    refl, refl_out, ref_s, tau_ap, tau0_out,&
    jacoAOD, AOD_max_steps,debug_flag, &
    pf_smooth_pix, pf_components_pix, &
    !iband, iprof, &
    i_channel, &
    n_aero_layers, n_sca, n_aod, n_pfc) !sa, tau_ap)

    implicit none
    !Real, Intent(In) :: sa, tau_ap
    integer :: kk
    Real, Dimension(1:N_AOD_max), Intent(In) :: AOD_max_steps   ! Aerosol parameters
    integer, intent(in) :: i_channel, n_sca, n_aod, n_pfc, n_aero_layers
    !integer(c_int), intent(in) :: iband, iprof
    integer :: iband, iprof
    real :: sca_tmp
    real, intent(in) :: saa, vaa, sca
    Real, dimension(N_Channels, n_aod), intent(in) :: ssa_pix
    Real, dimension(N_channels, n_sca, n_aod), intent(in) :: pf_smooth_pix
    Real, dimension(N_channels, n_pfc, n_aod), intent(in) :: pf_components_pix
    real, intent(in) :: tau_ap ! tau_ap is tau_force in start_exe
    real, intent(in) :: Sa, uv, us, refl,  ref_s
    Logical, Intent(In) :: debug_flag
    real, intent(out) :: refl_out, tau0_out, jacoAOD
    ! 2nd dim should be n_pf instead of 1...?
    real :: tau_apriori
    !Real :: Sy ! Matrices de covariance (ici ce sont des reels en fait)
    !Integer :: max_iter1,max_iter2
    !Real :: gamma_lm
    real :: gamma_lm, Sy
    !real, parameter :: alpha_lm = 2.
    real :: ref_tol_obs, ref_tol_calc!, rad_tol_obs, rad_tol_calc
    Real :: x_plus_1, xhi2, xhi2plus1
    Real :: jac, x, x_apriori
    Integer :: iter1, iter2
    Real :: tau_plus_1
    !Real :: channel

    interface
        subroutine calc_xi(x,x_apriori,jac,ref_tol_obs,ref_tol_calc,Sy,Sa,gamma_lm,x_plus_1)
           Real, Intent(In) :: x !xi a l'iteration precedente
           Real, Intent(In) :: x_apriori !x a priori
           Real, Intent(In) :: jac !valeur de la jacobienne en x (ici de dim 1)
           Real, Intent(In) :: ref_tol_obs !y dans le schema, la reflectance observee
           Real, Intent(In) :: ref_tol_calc !la reflectance calculee pour x
           Real, Intent(In) :: Sy !la variance de l'erreur de mesure (sur ref_tol_obs)
           Real, Intent(In) :: sa ! la variance de l'apriori
           Real, Intent(In) :: gamma_lm !la valeur prise par gamma
           Real, Intent(Out) :: x_plus_1 !la nouvelle valeur de xi
        end subroutine calc_xi
    end interface

    interface
        subroutine calc_xhi2(x,x_apriori,ref_tol_obs,ref_tol_calc,Sa,Sy,xhi2)
           real, intent(in) :: x !xi à l'iteration precedente
           real, intent(in) :: x_apriori !x a priori
           real, intent(in) :: ref_tol_obs !y dans le schema, la reflectance observee
           real, intent(in) :: ref_tol_calc !la reflectance calculee pour x
           real, intent(in) :: Sy !la variance de l'erreur de mesure (sur ref_tol_obs)
           real, intent(in) :: Sa ! la variance de l'apriori
           real, intent(out) :: xhi2 !la valeur du xhi2 à cette iteration
        end subroutine calc_xhi2
    end interface

    gamma_lm = 1.
    Sy = 0.0001 

    !real(FLOTSAM_REAL), parameter :: pi_flotsam=4.D0*ATAN(1.D0)
    call prep_flotsam(us, uv, saa, vaa, i_channel, iband, iprof)

    ! What is AOD_max_steps exactly?
    ! -> seems to be just the range of AOD values.. see aod_parameters.f90

    ! can define them in module scope, if I turn this into a module
    ref_tol_obs = refl
    
    !d=0.01
    tau_apriori = tau_ap
    x_apriori = tau_apriori ! why?
    x = x_apriori ! why?

    call execute_flotsam(AOD_max_steps, ssa_pix, pf_components_pix, &
        pf_smooth_pix, &
        us, sca, ref_s, &
        tau_apriori, i_channel, iband, &
        n_sca, n_pfc, n_aod, n_aero_layers, &
        .True., &
        !n_layers, &
        ref_tol_calc, jac)

    !print*, 'REF TOL CALC / RAD TOL CA ', ref_tol_calc, rad_tol_calc
    !print*, 'INIT FLOTSAM REFL / TAU: ', ref_tol_calc, tau_apriori

    call calc_xhi2(x, x_apriori, ref_tol_obs, ref_tol_calc, Sa, Sy, xhi2)
    IF (debug_flag) print*,"Initialisation de xhi2 =",xhi2

    !-----------------------  Levenberg - Marquardt  -------------------------
    ! --> Calcul de xi_plus_1
    call calc_xi(x, x_apriori, jac, ref_tol_obs, ref_tol_calc, Sy, Sa, gamma_lm, x_plus_1)
    x=x_plus_1
    IF (debug_flag) print*,"Initialisation de x_plus_1 =",x_plus_1 

    iter1 = 0 
    !print*, 'GAMMA / ALPHA: ', gamma_lm, alpha_lm
    do while(iter1 .LT. max_iter1)
        tau_plus_1 = x_plus_1
        
        call execute_flotsam(AOD_max_steps, ssa_pix, pf_components_pix, &
            pf_smooth_pix, &
            us, sca, ref_s, &
            tau_plus_1, i_channel, iband, &
            n_sca, n_pfc, n_aod, n_aero_layers, &
            .True., &
            !n_layers, &
            ref_tol_calc, jac)
        
        !print*, iter1, 'FLOTSAM OBS / REFL / TAU / JAC: ', ref_tol_obs, ref_tol_calc, tau_plus_1, jac
        call calc_xhi2(x_plus_1, x_apriori, ref_tol_obs, ref_tol_calc, Sa, Sy, xhi2plus1)
        
        iter2=0
        !print*, 'OUTER GAMMA / ALPHA: ', gamma_lm, alpha_lm
        do while(iter2 .LT. max_iter2 .and. xhi2plus1 .GE. xhi2)  ! .GE. = >=
        !do while(iter2 .LT. max_iter2)    
            iter2=iter2+1 
            call calc_xi(x, x_apriori, jac,ref_tol_obs, ref_tol_calc, &
                Sy, Sa, gamma_lm, x_plus_1)
            x=x_plus_1
            tau_plus_1 = x_plus_1
            call execute_flotsam(AOD_max_steps, ssa_pix, pf_components_pix, &
                pf_smooth_pix, &
                us, sca, ref_s, &
                tau_plus_1, i_channel, iband,&
                n_sca, n_pfc, n_aod, n_aero_layers, &
                .True., &
                !n_layers, &
                ref_tol_calc, jac)
                !print*, iter1, iter2, ' INNER FLOTSAM REFL / TAU: ', ref_tol_calc, tau_plus_1

            call calc_xhi2(x_plus_1, x_apriori, ref_tol_obs, ref_tol_calc, Sa, Sy, xhi2plus1)
        end do

        iter1=iter1+1
        xhi2=xhi2plus1
!      IF (debug_flag) print*,"xhi2=",xhi2
!
        gamma_lm = gamma_lm / alpha_lm
        call calc_xi(x, x_apriori, jac,ref_tol_obs, ref_tol_calc, Sy, Sa, gamma_lm, x_plus_1)
        x=x_plus_1

    end do 

    tau_plus_1 = x_plus_1
    tau0_out=tau_plus_1

    call execute_flotsam(AOD_max_steps, ssa_pix, pf_components_pix, &
        pf_smooth_pix, &
        us, sca, ref_s, &
        tau_plus_1, i_channel, iband, &
        n_sca, n_pfc, n_aod, n_aero_layers, &
        .False., &
        !n_layers, &
        ref_tol_calc, jac)

    !print*, 'LM DONE'

    call free_flotsam_profiles(iband, iprof)


end subroutine LM_flotsam

subroutine prep_flotsam(us, uv, saa, vaa, i_channel_base1, iband, iprof)
    integer, intent(in) :: i_channel_base1
    real, intent(in) :: us, uv, saa, vaa
    integer, intent(out) :: iband, iprof
    integer(c_int) :: ichan_c, iprof_c, iband_c
    real(FLOTSAM_REAL) :: azim, us_c, uv_c
    integer :: callstat
    
    ichan_c = i_channel_base1 - 1

    ! band profile holds the optical properties of atmospheric gases for a
    ! particular atmospheric profile
    iband_c = flotsam_new_band_profile()

    ! background profile holds the background properties of an atmospheric profile:
    ! pressure, temperature and the mixing ratios of absorbing gases
    iprof_c = flotsam_new_background_profile()

    callstat = flotsam_set_edge_pressure(iprof_c, n_layers, tol_edge_pressure)
    callstat = flotsam_init_band_profile(iband_c, ichan_c, iprof_c)

    azim = vaa - saa
    us_c = us
    uv_c = uv
    callstat = flotsam_set_geometry(iband_c, us_c, uv_c, azim)

    iband = iband_c
    iprof = iprof_c

end subroutine prep_flotsam


subroutine free_flotsam_profiles(iband, iprof)
    integer, intent(in) :: iband, iprof
    integer(c_int) :: iband_c, iprof_c
    integer :: callstat

    iband_c = iband
    iprof_c = iprof
    callstat = flotsam_free_band_profile(iband_c)
    callstat = flotsam_free_background_profile(iprof_c)

end subroutine free_flotsam_profiles

subroutine run_flotsam_daily(ssa_pix, pf_components_pix, &
    pf_smooth_pix, &
    uv, us, saa, vaa, sca, ref_s, &
    aod, AOD_max_steps, ph_val_msa, &
    ssa_tilde, aod_tilde, &
    i_channel, &
    n_aero_layers, n_sca, n_aod, n_pfc, &  
    rho_1, trans_coeff, r_ms)

    Real, Dimension(1:N_AOD_max), Intent(In) :: AOD_max_steps   ! Aerosol parameters
    Real, dimension(N_Channels, n_aod), intent(in) :: ssa_pix
    Real, dimension(N_channels, n_sca, n_aod), intent(in) :: pf_smooth_pix
    Real, dimension(N_channels, n_pfc, n_aod), intent(in) :: pf_components_pix
    integer, intent(in) :: n_sca, n_pfc, n_aod, n_aero_layers, i_channel
    real, intent(in) :: ssa_tilde, aod_tilde !, phase_tilde
    !integer, intent(in) :: n_layers
    integer(c_int) :: iband_c
    real, intent(in) :: aod, us, uv, sca, ref_s, saa, vaa, ph_val_msa
    !real, intent(in) :: ssa_tilde_val
    real, Intent(InOut)  :: r_ms, trans_coeff, rho_1
    real :: rho_tol_0, rho_tol_1, dummy_real, r_ss, ssa_val, pf_interp
    !real :: r_ss_tilde, rho_1_tilde
    real :: pf_cmp_pix_contiguous(n_pfc,1)
    integer :: iprof, iband

    call prep_flotsam(us, uv, saa, vaa, i_channel, iband, iprof)

    ! call 1: get rho_tol_0
    call execute_flotsam(AOD_max_steps, &
        ssa_pix, pf_components_pix, &
        pf_smooth_pix, &
        us, sca, ref_s, &
        aod, i_channel, iband, &
        n_sca, n_pfc, n_aod, n_aero_layers, &
        .False., &
        !n_layers, &
        rho_tol_0, dummy_real)

    !print*, 'RHO TOL 0: ', rho_tol_0

    ! call 2: get rho_tol_1
    call execute_flotsam(AOD_max_steps, &
        ssa_pix, pf_components_pix, &
        pf_smooth_pix, &
        us, sca, 0., &
        aod, i_channel, iband, &
        n_sca, n_pfc, n_aod, n_aero_layers, &
        .False., &
        !n_layers, &
        rho_tol_1, dummy_real)

    !print*, 'RHO TOL 1: ', rho_tol_1

    !rho_1 = 1./(4.*(us+uv))*(1.-exp(-aod*(1./us+1./uv)))
    rho_1 = 1./(4.*(us+uv))*(1.-exp(-aod_tilde*(1./us+1./uv)))
    !print*, 'RHO 1: ', rho_1

    !call get_interp_aer_parms(sca, aod, AOD_max_steps, &
    !    ssa_pix, pf_smooth_pix, pf_components_pix, &
    !    n_sca, n_pfc, n_aod, i_channel, &
    !    pf_interp, ssa_val, pf_cmp_pix_contiguous)

    !r_ss = ssa_val*pf_interp*rho_1
    !print*, 'FLOTSAM SSA: ', ssa_val
    !r_ss = ssa_val*ph_val_msa*rho_1
    r_ss = ssa_tilde*ph_val_msa*rho_1
    !print*, 'R SS: ', r_ss
    !rho_1_tilde = 1./(4.*(us+uv))*(1.-exp(-aod_tilde*(1./us+1./uv)))
    !r_ss_tilde = ssa_tilde*phase_tilde*rho_1_tilde
    !print*, 'R SS ~: ', r_ss_tilde
    
    r_ms = rho_tol_1 - r_ss
    trans_coeff = (rho_tol_0 - rho_tol_1) / ref_s

    !print*, 'rho 1, tc, r_ms: ', rho_1, trans_coeff, r_ms
    call free_flotsam_profiles(iband, iprof)

end subroutine run_flotsam_daily

subroutine get_interp_aer_parms(sca, aod, AOD_max_steps, &
    ssa_pix, pf_smooth_pix, pf_components_pix, &
    n_sca, n_pfc, n_aod, i_channel, &
    pf_interp, ssa_val, pf_cmp_pix_contiguous)
    Real, Dimension(1:N_AOD_max), Intent(In) :: AOD_max_steps   ! Aerosol parameters
    Real, dimension(N_Channels, n_aod), intent(in) :: ssa_pix
    Real, dimension(N_channels, n_sca, n_aod), intent(in) :: pf_smooth_pix
    Real, dimension(N_channels, n_pfc, n_aod), intent(in) :: pf_components_pix
    real, intent(in) :: sca, aod
    integer, intent(in) :: n_sca, n_pfc, n_aod, i_channel
    real, intent(out) :: pf_interp, ssa_val, pf_cmp_pix_contiguous(n_pfc,1)
    real(FLOTSAM_REAL) :: pf_smoo_pix_contiguous(n_sca,1), sca_c
    integer :: i_aod    

    sca_c = sca
    
    ! retrieving position from AOD steps of dynamic properties
    i_aod = MINLOC(ABS(AOD_max_steps-aod), DIM=1)
   
    pf_cmp_pix_contiguous = reshape(&
        pf_components_pix(i_channel,:, i_aod),&
        (/n_pfc,1/))
    pf_smoo_pix_contiguous = reshape(&
        pf_smooth_pix(i_channel,:, i_aod),&
        (/n_sca,1/))

    ! technically I could just use the interpolated phase value as it is 
    ! calculated for MSA.. maybe in the future
    pf_interp = flotsam_interp_phase_func(&
            n_sca, pf_smoo_pix_contiguous, sca_c)

    ssa_val = ssa_pix(i_channel, i_aod)

end subroutine get_interp_aer_parms

subroutine execute_flotsam(AOD_max_steps, ssa_pix, pf_components_pix, &
    pf_smooth_pix, &
    mu_sun, sca, ref_s, &
    aod, i_channel, iband, &
    n_sca, n_pfc, n_aod, n_aero_layers, &
    do_jac, &
    !n_layers, &
    ref_tol_calc, jac)

    ! Does a single execute of FLOTSAM, or a double execute to compute the Jacobian    
    ! -------------------------------------------------------------------------

    Real, Dimension(1:N_AOD_max), Intent(In) :: AOD_max_steps   ! Aerosol parameters
    Real, dimension(N_Channels, n_aod), intent(in) :: ssa_pix
    Real, dimension(N_channels, n_sca, n_aod), intent(in) :: pf_smooth_pix
    Real, dimension(N_channels, n_pfc, n_aod), intent(in) :: pf_components_pix
    integer, intent(in) :: n_sca, n_pfc, n_aod, n_aero_layers, i_channel
    !integer, intent(in) :: n_layers
    integer, intent(in) :: iband
    integer(c_int) :: iband_c
    real, intent(in) :: aod, mu_sun, sca, ref_s
    logical, intent(in) :: do_jac
    real, intent(out) :: ref_tol_calc, jac
    real :: aod_jac, r_jac_1, r_jac_2, jac_fac
    real(FLOTSAM_REAL) :: sca_c, pf_interp
    Real(FLOTSAM_REAL) :: pf_smoo_pix_contiguous(n_sca,1), pf_cmp_pix_contiguous(n_pfc,1)
    real(FLOTSAM_REAL), dimension(n_aero_layers) :: pf_interp_array, aod_array, ssa_array
    real(FLOTSAM_REAL), dimension(n_aero_layers) :: ref_s_array
    real(FLOTSAM_REAL) :: radiance, d_albedo, d_ssa, d_pf, d_pfc, d_od
    !integer(c_int) :: ichan !iprof
    integer :: i_aod, callstat, i_aod_jac

    iband_c = iband
   
    ! if several aerosol layers, dim-2 size might be >1...?
    !  (should be n_pf instead of 1 in that case...?)
    ! phase_function = reshape(pf_smooth, (/n_sca,1/))
   
    sca_c = sca
    !print*, mu_sun_flotsam, mu_view_flotsam, azim   
    
    ! retrieving position from AOD steps of dynamic properties
    i_aod = MINLOC(ABS(AOD_max_steps-aod), DIM=1)
   
    pf_cmp_pix_contiguous = reshape(&
        pf_components_pix(i_channel,:, i_aod),&
        (/n_pfc,1/))
    pf_smoo_pix_contiguous = reshape(&
        pf_smooth_pix(i_channel,:, i_aod),&
        (/n_sca,1/))

    ! technically I could just use the interpolated phase value as it is 
    ! calculated for MSA.. maybe in the future
    pf_interp = flotsam_interp_phase_func(&
            n_sca, pf_smoo_pix_contiguous, sca_c)
    !print*, 'FLOTSAM INIT PHF: ', aod, i_aod, sca_c, pf_interp

    ! wait... what is it now....? need to double check 
    ! if passing only one phase function, the shape does not seem to matter
    ! this is how it is done in ECRAD:
    ! allocate(pf_components(n_pf_components, nlev)) [nlev: number of model levels]
    ! in our case this should be (n_pfc, 1) then, I think...?
    !pfc_array(1, :) = pf_cmp_pix_contiguous(:,1)

    ref_s_array(1) = ref_s
    aod_array(1) = aod
    ssa_array(1) = ssa_pix(i_channel, i_aod)
    pf_interp_array(1) = pf_interp

    !print*, 'FLOTSAM INPUTS: '
    !print*, iband, n_albedo_components, ref_s_array, aod_array, ssa_array, pf_interp_array, n_layers, loc
    !print*, real(pf_cmp_pix_contiguous, FLOTSAM_REAL)
    radiance = 0.
    if (.false.) then 
        ! i don't think we can use flotsam jacobian
        ! if we want to have dynamic aerosol properties...
        callstat = flotsam_reflectance_jacobian(iband, n_albedo_components, &
            ref_s_array, n_layers, loc, aod_array, ssa_array,&
            pf_interp_array, &
            pf_cmp_pix_contiguous, &
            !pfc_array, &
            radiance, d_albedo, d_od, d_ssa, d_pf, d_pfc)
    else
        !print*, 'FLOTSAM INPUTS: ', i_aod, sca_flotsam, pf_interp, ssa_array(1), ref_s, aod
        callstat = flotsam_reflectance(iband, n_albedo_components, &
            ref_s_array, n_layers, loc, aod_array, ssa_array,&
            pf_interp_array, &
            pf_cmp_pix_contiguous, &
            !pfc_array, &
            radiance)


        ref_tol_calc = radiance * pi_lm / mu_sun
        !print*, 'INIT FLOTSAM REFL / TAU: ', ref_tol_calc, aod
        !rad_tol_calc = radiance
        r_jac_1 = ref_tol_calc

        if (do_jac) then
            if (i_aod .eq. N_AOD_max) then
                i_aod_jac = i_aod - 1
                aod_jac = aod - d_aod
                jac_fac = -1.   
            else
                i_aod_jac = i_aod + 1
                aod_jac = aod + d_aod
                jac_fac = 1.
            endif

            ! need to change this later..
            pf_cmp_pix_contiguous = reshape(&
                pf_components_pix(i_channel,:, i_aod_jac),&
                (/n_pfc,1/))
            pf_smoo_pix_contiguous = reshape(&
                pf_smooth_pix(i_channel,:, i_aod_jac),&
                (/n_sca,1/))

            pf_interp = flotsam_interp_phase_func(&
                    n_sca, pf_smoo_pix_contiguous, sca_c)

            aod_array(1) = aod_jac
            ssa_array(1) = ssa_pix(i_channel, i_aod_jac)
            pf_interp_array(1) = pf_interp

            callstat = flotsam_reflectance(iband, n_albedo_components, &
                ref_s_array, n_layers, loc, aod_array, ssa_array,&
                pf_interp_array, &
                pf_cmp_pix_contiguous, &
                !pfc_array, &
                radiance)
            r_jac_2 = radiance * pi_lm / mu_sun

            jac = jac_fac * (r_jac_2 - r_jac_1) / d_aod
            !print*, 'JACO :', rad_jac_1, rad_jac_2, d_aod, jac_fac, jac
        else 
            ! return dummy result for jacobian
            jac = 0.
        endif 

    endif

end subroutine execute_flotsam

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

subroutine init_flotsam_channel()
    integer :: ichan

    ! Create channel (wavelength-dependent information)
    ! -------------------------------------------------------------------
    ! use ..._vacuum if rayleigh scattering is OFF
    ichan = flotsam_new_channel_vacuum()

end subroutine init_flotsam_channel

! xi is phase function here
! (function was just for testing..)
subroutine calc_flotsam_xi(saa, vaa, us, uv, xi)
    real, intent(in) :: saa, vaa, us, uv
    real(FLOTSAM_REAL) :: azim, us_flotsam, uv_flotsam, xi_flotsam
    real, intent(out) :: xi

    azim = vaa - saa
    us_flotsam = us
    uv_flotsam = uv
    xi_flotsam = flotsam_scattering_angle(us_flotsam, uv_flotsam, azim)
    xi = xi_flotsam

end subroutine calc_flotsam_xi
 

end module flotsam_interface

