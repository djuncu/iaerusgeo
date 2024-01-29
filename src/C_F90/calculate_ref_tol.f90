Subroutine calculate_ref_tol(ref_s,p_hg_tilde,ssa_tilde,g_tilde,tau0_tilde,albed_b_in,uv,us,ref_tol)
  Implicit None
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

End Subroutine calculate_ref_tol
