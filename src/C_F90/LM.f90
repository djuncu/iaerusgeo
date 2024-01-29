!-----------------------------------------------------------------!
!                                                                 !
!    Estimation optimale par la méthode de Levenberg Marquardt    !
!                                                                 !
!-----------------------------------------------------------------!
Subroutine LM (Sa, g_tilde,ssa,ssa_tilde,eta,uv,us,phf,refl,refl_out,ref_s,albed_b_in,tau_ap,tau0_out,jacoAOD,AOD_max_steps,debug_flag)

   Use py_ifc

   Implicit None 

   Real, Dimension(1:N_AOD_max), Intent(In) :: g_tilde, ssa, ssa_tilde, eta, phf, AOD_max_steps   ! Aerosol parameters
   Real, Intent(In)                         :: refl,ref_s,albed_b_in ! Surface reflectances and albedo
   Real, Intent(In)                         :: Sa, tau_ap
   Real, Intent(In)                         :: uv, us ! Geometry
   Logical, Intent(In)                      :: debug_flag
   Real, Intent(Out)                        :: tau0_out         ! AOD
   Real, Intent(Out)                        :: jacoAOD         ! AOD solution for which fitting error is lower than maximum error
   Real, Intent(Out)                        :: refl_out
   Real :: tau_apriori, tau_apriori_tilde ! AOD value a priori
   Real :: Sy ! Matrices de covariance (ici ce sont des reels en fait)
   Integer ::max_iter1,max_iter2
   Real :: gamma_lm,alpha_lm
   Real :: x_plus_1, ref_tol_calc,xhi2,xhi2plus1
   Real :: d,jac,ref_tol,x,x_apriori
   Integer :: iter1, iter2, kk
   Real :: tau_plus_1, phf_tilde
   Real :: channel,ref_tol_obs

   !----------------------------------------------------------------------------
   !         --------------->  RTM   <----------------
   !----------------------------------------------------------------------------

   ! interface of instantaneous retrieval
   Interface
      Subroutine calculate_ref_tol(ref_s,phf_tilde,ssa_tilde_kk,g_tilde_kk,tau0_tilde,albed_b_in,uv,us,ref_tol)
      Implicit None
	   Real, Intent(In)     :: ref_s,phf_tilde,ssa_tilde_kk,g_tilde_kk,tau0_tilde,albed_b_in,uv,us
	   Real, Intent(InOut)    :: ref_tol
      End Subroutine calculate_ref_tol
   End Interface

   ref_tol_obs=refl
   channel=1
   max_iter1=8 ! similar results for 3 iterations for all models except dust
   max_iter2=8
   gamma_lm=1
   alpha_lm=2
   Sy=0.0001
   !Sa=0.1
   d=0.01
   !tau_apriori=0.15
   tau_apriori=tau_ap

   ! retrieving position from AOD steps of dynamic properties
   kk = MINLOC(ABS(AOD_max_steps-tau_apriori), DIM=1)
   call calculate_tau_tilde(tau_apriori, ssa(kk), eta(kk), tau_apriori_tilde)
   x_apriori=tau_apriori_tilde
   x=x_apriori

   !----------------------------------------------------------------------------
   !         --------------->  levenberg-marquardt   <----------------
   !----------------------------------------------------------------------------

   !============== Calcul variables tilde =================================
   !call calculate_tilde_variables(4, angle, tau_apriori, 1, eta_out, g_tilde, ssa_out, ssa_tilde, p_tilde)
   phf_tilde = phf(kk) / (1.-eta(kk)) ! truncation
   !print*, 'MSA INIT PHF / SSA: ', phf(kk), ssa_tilde(kk)

   !======================================  calcul jacobienne  ==================================
   call calculate_jac_ref_tol(ref_s,phf,ssa,ssa_tilde,g_tilde,eta,AOD_max_steps,tau_apriori,albed_b_in,uv,us,kk,d,Jac)
   jac=Jac

   !-----------------------  Calcul du transfert radiatif   -------------------------
   !print*, 'MSA INPUTS: ', phf(kk)
   ! is this regular ssa...?
   !print*, ssa_tilde(kk) / (1. - eta(kk) * (1-ssa_tilde(kk)))
   call calculate_ref_tol(ref_s,phf_tilde,ssa_tilde(kk),g_tilde(kk),x_apriori,albed_b_in,uv,us,ref_tol)
   ref_tol_calc = ref_tol
   !print*, 'INIT MSA REFL / TAU: ', ref_tol_calc, tau_apriori

   call calc_xhi2(x, x_apriori, ref_tol_obs, ref_tol_calc,Sa, Sy, xhi2)
   IF (debug_flag) print*,"Initialisation de xhi2 =",xhi2

   !-----------------------  Schéma de Levenberg - Marquardt  -------------------------
   ! --> Calcul de xi_plus_1
   call calc_xi(x, x_apriori, jac, ref_tol_obs, ref_tol_calc, Sy, Sa, gamma_lm, x_plus_1)
   x=x_plus_1
   IF (debug_flag) print*,"Initialisation de x_plus_1 =",x_plus_1
  
   !-----------------------  Fonction de cout a minimiser  -------------------------
  
   !----------------------------------------------------------------------------
   !           --------------->  Partie itération   <----------------
   !----------------------------------------------------------------------------

   iter1 = 0
   iter2 = 0

   do while(iter1 .LT. max_iter1)
      call calculate_tau(x_plus_1,ssa(kk), eta(kk), tau_plus_1)
      kk = MINLOC(ABS(AOD_max_steps-tau_plus_1), DIM=1)
      phf_tilde = phf(kk) / (1.-eta(kk)) ! truncation
      call calculate_jac_ref_tol(ref_s,phf,ssa,ssa_tilde,g_tilde,eta,AOD_max_steps,tau_plus_1,albed_b_in,uv,us,kk,d,Jac)
      jac=Jac
      call calculate_ref_tol(ref_s,phf_tilde,ssa_tilde(kk),g_tilde(kk),x_plus_1,albed_b_in,uv,us,ref_tol)
      ref_tol_calc = ref_tol
      !print*, iter1, 'MSA OBS / REFL / TAU / JAC: ', ref_tol_obs, ref_tol_calc, tau_plus_1, jac
      call calc_xhi2(x_plus_1, x_apriori, ref_tol_obs, ref_tol_calc, Sa, Sy, xhi2plus1)
      iter2=0
        IF (debug_flag) THEN
         write(*,*) "boucle iter 1"
         write(*,*) "iter1=", iter1
         write(*,*) "Xhi2plus1=", xhi2plus1
         write(*,*) "Gamma=", gamma_lm
         write(*,*) "Tau_est=",x
         write(*,*) "-------------"
        ENDIF

         do while(iter2 .LT. max_iter2 .and. xhi2plus1 .GE. xhi2)    ! .GE. = >=

           IF (debug_flag) THEN 
            write(*,*) "boucle iter 2"
            write(*,*) "iter2=", iter2
            write(*,*) "Xhi2plus=", xhi2plus1
            write(*,*) "Gamma=", gamma_lm
            write(*,*) "Tau_est=",x
            write(*,*) "-------------"
           ENDIF
            gamma_lm=gamma_lm*alpha_lm
            iter2=iter2+1   
            ! xi_plus_1 ne fait pas baisser la fonction de cout donc calcul d'un nouveau xi_plus_1 avec le nouveau gamma_lm 
            call calc_xi(x, x_apriori,jac, ref_tol_obs, ref_tol_calc, Sy, Sa, gamma_lm, x_plus_1)
            x=x_plus_1
            ! calcul de la nouvelle reflectance estimee 
            call calculate_tau(x_plus_1,ssa(kk), eta(kk), tau_plus_1)
            kk = MINLOC(ABS(AOD_max_steps-tau_plus_1), DIM=1)
            phf_tilde = phf(kk) / (1.-eta(kk)) ! truncation
            call calculate_ref_tol(ref_s,phf_tilde,ssa_tilde(kk),g_tilde(kk),x_plus_1,albed_b_in,uv,us,ref_tol_calc)
            IF (debug_flag) print*,'ref_tol_calc : ',ref_tol_calc,ref_tol_obs
            !print*, iter1, iter2, 'MSA REFL / TAU: ', ref_tol_calc, tau_plus_1
            call calculate_jac_ref_tol(ref_s,phf,ssa,ssa_tilde,g_tilde,eta,AOD_max_steps,tau_plus_1,albed_b_in,uv,us,kk,d,Jac)
            call calc_xhi2(x_plus_1, x_apriori, ref_tol_obs, ref_tol_calc, Sa, Sy, xhi2plus1)
            xhi2plus1=xhi2plus1
         end do

      iter1=iter1+1
      !x = x_plus_1      
      xhi2=xhi2plus1
      IF (debug_flag) print*,"xhi2=",xhi2

      gamma_lm = gamma_lm / alpha_lm
      call calc_xi(x, x_apriori, jac,ref_tol_obs, ref_tol_calc, Sy, Sa, gamma_lm, x_plus_1)
      x=x_plus_1
      !call calculate_tau(x_plus_1,ssa(kk), eta(kk), tau_plus_1)
      !kk = MINLOC(ABS(AOD_max_steps-tau_plus_1), DIM=1)
      !phf_tilde = phf(kk) / (1.-eta(kk)) ! truncation
      !call calculate_jac_ref_tol(ref_s,phf,ssa,ssa_tilde,g_tilde,eta,AOD_max_steps,x_plus_1,albed_b_in,uv,us,kk,d,Jac)
      !call calculate_ref_tol(ref_s,phf_tilde,ssa_tilde(kk),g_tilde(kk),x_plus_1,albed_b_in,uv,us,ref_tol)
      !ref_tol_calc=ref_tol
    
   end do 
   call calculate_tau(x_plus_1,ssa(kk), eta(kk), tau_plus_1)
   tau0_out=tau_plus_1
   kk = MINLOC(ABS(AOD_max_steps-tau_plus_1), DIM=1)
   phf_tilde = phf(kk) / (1.-eta(kk)) ! truncation
   call calculate_jac_ref_tol(ref_s,phf,ssa,ssa_tilde,g_tilde,eta,AOD_max_steps,tau_plus_1,albed_b_in,uv,us,kk,d,Jac)
   call calculate_ref_tol(ref_s,phf_tilde,ssa_tilde(kk),g_tilde(kk),x_plus_1,albed_b_in,uv,us,ref_tol)
   ref_tol_calc=ref_tol
   refl_out=ref_tol_calc
   jacoAOD=jac

   return

End Subroutine LM

!------------------------------------------------------------------------------------------------
!                       --------------->  SUBROUTINES   <----------------
!------------------------------------------------------------------------------------------------
!   !============== Calcul variables tilde =================================
!   !call calculate_tilde_variables(4, angle, tau_apriori, 1, eta_out, g_tilde, ssa_out, ssa_tilde, p_tilde)
!   phf_tilde = phf(kk) / (1.-eta(kk)) ! truncation!

!   !======================================  calcul jacobienne  ==================================
!   call calculate_jac_ref_tol(ref_s,phf,ssa,ssa_tilde,g_tilde,eta,AOD_max_steps,tau_apriori,albed_b_in,uv,us,kk,d,Jac)
!   jac=Jac
!
!   !-----------------------  Calcul du transfert radiatif   -------------------------
!   call calculate_ref_tol(ref_s,phf_tilde,ssa_tilde(kk),g_tilde(kk),x_apriori,albed_b_in,uv,us,ref_tol)
!   ref_tol_calc = ref_tol

!end subroutine calc_msa_jac_ref_tol


Subroutine calculate_jac_ref_tol(ref_s,phf,ssa,ssa_tilde,g_tilde,eta,AOD_max_steps,tau0_out,albed_b_in,uv,us,kk,d,Jac)

   Use py_ifc

   Implicit none 

   ! declaration des variables
   Real, Dimension(1:N_AOD_max), Intent(In) :: g_tilde, ssa, ssa_tilde, eta, phf, AOD_max_steps   ! Aerosol parameters
   Real, Intent(In) :: ref_s, tau0_out, albed_b_in, uv, us, d
   Integer, Intent(In) :: kk
   Real, intent(Out) :: Jac

   Real :: phf_tilde, tau0_tilde, ref_tol
   Real :: c1,c2
   Integer :: kkk

   ! jacobian is calculated at the AOD step of aerosol properties (aka. kk and kk+1, which is 0.01 for the moment)

   kkk = kk

   ! respecting AOD limits of LUTs
   IF (kk .GE. N_AOD_max) kkk = kkk-1
   IF (kk .LE. 0) kkk = kkk+1

   ! truncated phase function from Katsev et al.
   phf_tilde = phf(kkk) / (1.-eta(kkk)) ! truncation
   ! calculation of tilde AOD
   tau0_tilde = (1.-ssa(kkk)*eta(kkk))*tau0_out
   ! retrieval for current slot
   call calculate_ref_tol(ref_s,phf_tilde,ssa_tilde(kkk),g_tilde(kkk),tau0_tilde,albed_b_in,uv,us,ref_tol)
   c1=ref_tol

   ! truncated phase function from Katsev et al.
   phf_tilde = phf(kkk+1) / (1.-eta(kkk+1)) ! truncation
   ! calculation of tilde AOD
   tau0_tilde = (1.-ssa(kkk+1)*eta(kkk+1))*(tau0_out+d)
   ! retrieval for current slot
   call calculate_ref_tol(ref_s,phf_tilde,ssa_tilde(kkk+1),g_tilde(kkk+1),tau0_tilde,albed_b_in,uv,us,ref_tol)
   c2=ref_tol      

   Jac=(c2-c1)/d

   return

end subroutine calculate_jac_ref_tol

Subroutine calculate_tau(tau_tilde,w,eta,tau)

   Implicit None
   Real, Intent(In) :: tau_tilde
   Real, Intent(In) :: w, eta
   Real, Intent(Out) :: tau
   tau=tau_tilde/(1-w*eta)
   return

End Subroutine

Subroutine calculate_tau_tilde(tau,w,eta,tau_tilde)

   Implicit None
   Real, Intent(In) :: tau
   Real, Intent(In) :: w, eta
   Real, Intent(Out) :: tau_tilde
   tau_tilde=tau*(1-w*eta)
   return

End Subroutine

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

   x_plus_1=x_apriori+(1./(jac*(1./Sy)*jac+(1./Sa)+gamma_lm/Sa))*(jac*(1./Sy)* &
      & ((ref_tol_obs-ref_tol_calc)+jac*(x-x_apriori))+gamma_lm*(1./Sa)*(x-x_apriori))

   return

end subroutine calc_xi

subroutine calc_xhi2(x,x_apriori,ref_tol_obs,ref_tol_calc,Sa,Sy,xhi2)

   real, intent(in) :: x !xi à l'iteration precedente
   real, intent(in) :: x_apriori !x a priori
   real, intent(in) :: ref_tol_obs !y dans le schema, la reflectance observee
   real, intent(in) :: ref_tol_calc !la reflectance calculee pour x
   real, intent(in) :: Sy !la variance de l'erreur de mesure (sur ref_tol_obs)
   real, intent(in) :: Sa ! la variance de l'apriori
   real, intent(out) :: xhi2 !la valeur du xhi2 à cette iteration

   xhi2=((ref_tol_obs-ref_tol_calc)*(1/Sy)*(ref_tol_obs-ref_tol_calc))+(x-x_apriori)*(1/Sa)*(x-x_apriori)

   return

end subroutine calc_xhi2

! JUST TESTING HERE FOR NOW
subroutine prepare_ref_tol_flotsam_local(&
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
   Use iso_c_binding

   implicit none
#include <flotsam.inc>
   
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
    ! edge pressure for single layer case
    real(FLOTSAM_REAL), dimension(2) :: tol_edge_pressure = (/0, 101300/)
    ! IMPORTANT: this index vector is intended for C++, must be 0-based!
    ! in theory depends on number of atm. layers, which is 1 in case of TOL
    integer(c_int), dimension(1) :: loc = (/0/)
   
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
   !print*, sca
   sca_aerus = sca
   
   rstatus = flotsam_set_geometry(iband, mu_sun_flotsam, mu_view_flotsam, azim)
   
   !pf_interp = flotsam_interp_phase_func(&
   !    n_sca, real(phase_function, FLOTSAM_REAL), sca)
   !print*, pf_interp

end subroutine prepare_ref_tol_flotsam_local
