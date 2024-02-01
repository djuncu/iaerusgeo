Module brdfmodels

  ! size of the look-up-table for hemispherical integrals I1 and I2
  Integer, Parameter :: n_table = 18
  ! step width of solar zenith angle (in degrees) in the look-up-table
  Real, Parameter :: theta_step = 5.

  ! look-up-table for hemispherical integrals I1 and I2
  Real, Dimension(0:n_table-1,1:3) :: hint
  Real, Dimension(0:n_table-1,1:3) :: hint_land
  Real, Dimension(0:n_table-1,1:3) :: hint_ocean

  ! bi-hemispherical integrals of the model kernel functions
  Real, Dimension(0:3)             :: bihi_calc
  Real, Dimension(0:3)             :: bihi_calc_land
  Real, Dimension(0:3)             :: bihi_calc_ocean

Contains

  Function brdfmodel(theta_obs, theta_sun, phi_obs, phi_sun, phi_del, wind_speed, phiw, I, model, ocean_flag, instantaneous)
    Implicit None

    Real, Intent(In)     :: theta_obs, theta_sun, phi_obs, phi_sun, phi_del, wind_speed, phiw
    Integer, Intent(In)  :: model, I
    Logical, Intent(In)  :: ocean_flag, instantaneous ! flag for land/water mask
    Real, Dimension(0:2) :: brdfmodel

    If (ocean_flag) Then
       brdfmodel = brdfmodel_ocean(theta_obs, theta_sun, phi_obs, phi_sun, phi_del, wind_speed, phiw, I, model)
    Else
       brdfmodel = brdfmodel_land(theta_obs, phi_del, theta_sun, model, instantaneous)
    Endif

  End Function brdfmodel

  Function brdfmodel_land(theta_obs, phi_del, theta_sun, model_land, instantaneous)
    Implicit None
  
    Real, Intent(In)     :: theta_obs, phi_del, theta_sun
    Integer, Intent(In)  :: model_land
    Real, Dimension(0:2) :: brdfmodel_land
    Logical, Intent(In)  :: instantaneous

    Select Case (model_land)
     Case (1)
       brdfmodel_land = roujean(theta_obs, phi_del, theta_sun)
     Case (2)
       brdfmodel_land = rtls(theta_obs, phi_del, theta_sun)
     Case (3)
       if (instantaneous) then
         brdfmodel_land = lirosshotspot(theta_obs, phi_del, theta_sun, 1.5)
       Else
         brdfmodel_land = lirosshotspot(theta_obs, phi_del, theta_sun, 1.5)
       Endif
     Case Default
       print *, 'wrong BRDF model for land'
       stop 
    End Select

  End Function brdfmodel_land

  Function brdfmodel_ocean(theta_obs, theta_sun, phi_obs, phi_sun, phi_del, wind_speed, phiw, I, model_ocean)
    Implicit None
  
    Real, Intent(In)     :: theta_obs, theta_sun, phi_obs, phi_sun, phi_del, wind_speed, phiw
    Integer, Intent(In)  :: model_ocean, I
    Real, Dimension(0:2) :: brdfmodel_ocean
    
    Select Case (model_ocean)
     Case (1)
       brdfmodel_ocean = coxmunk(theta_obs, theta_sun, phi_obs, phi_sun, phi_del, wind_speed, phiw, I)
     Case Default
       print *, 'wrong BRDF model for ocean'
       stop 
    End Select

  End Function brdfmodel_ocean

  Function brdfmodel_aerosol(theta_obs, theta_sun, phi_obs, phi_sun, phi_del, wind_speed, phiw, I, model, ocean_flag, tau0, g_, ssa_, trans_coeff, phFunc, instantaneous)
    Implicit None

    Real, Intent(In)      :: theta_obs, theta_sun, phi_obs, phi_sun, phi_del, wind_speed, phiw, tau0, g_, ssa_, trans_coeff, phFunc
    Integer, Intent(In)  :: model, I
    Logical, Intent(In)  :: ocean_flag, instantaneous ! flag for land/water mask
    Real, Dimension(0:3) :: brdfmodel_aerosol

    If (ocean_flag) Then
      brdfmodel_aerosol = brdfmodel_ocean_aerosol(theta_obs, theta_sun, phi_obs, phi_sun, phi_del, wind_speed, phiw, I, model, tau0, g_, ssa_, trans_coeff, phFunc)
    Else
      brdfmodel_aerosol = brdfmodel_land_aerosol(theta_obs, phi_del, theta_sun, model, tau0, g_, ssa_, trans_coeff, phFunc, instantaneous)
    Endif

  End Function brdfmodel_aerosol

  Function brdfmodel_land_aerosol(theta_obs, phi_del, theta_sun, model_land, tau0, g_, ssa_, trans_coeff, phFunc, instantaneous)
    Implicit None
  
    Real, Intent(In)     :: theta_obs, phi_del, theta_sun, tau0, g_, ssa_, trans_coeff, phFunc
    Integer, Intent(In)  :: model_land
    Real, Dimension(0:3) :: brdfmodel_land_aerosol
    Logical, Intent(In)  :: instantaneous

    Select Case (model_land)
     Case (1)
       brdfmodel_land_aerosol = roujean_aerosol(theta_obs, phi_del, theta_sun, tau0, g_, ssa_, trans_coeff, phFunc)
     Case (2)
       brdfmodel_land_aerosol = rtls_aerosol(theta_obs, phi_del, theta_sun, tau0, g_, ssa_, trans_coeff, phFunc)
     Case (3)
       if (instantaneous) then
         brdfmodel_land_aerosol = lirosshotspot_aerosol(theta_obs, phi_del, theta_sun, tau0, g_, ssa_, trans_coeff, phFunc, 1.5)
       Else
         brdfmodel_land_aerosol = lirosshotspot_aerosol(theta_obs, phi_del, theta_sun, tau0, g_, ssa_, trans_coeff, phFunc, 1.5)
       Endif
     Case Default
       print *, 'wrong BRDF model for land'  
       stop 
    End Select

  End Function brdfmodel_land_aerosol

  Function brdfmodel_ocean_aerosol(theta_obs, theta_sun, phi_obs, phi_sun, phi_del, wind_speed, phiw, I, model_ocean, tau0, g_, ssa_, trans_coeff, phFunc)
    Implicit None
  
    Real, Intent(In)     :: theta_obs, theta_sun, phi_obs, phi_sun, phi_del, wind_speed, phiw, tau0, g_, ssa_, trans_coeff, phFunc
    Integer, Intent(In)  :: model_ocean, I
    Real, Dimension(0:3) :: brdfmodel_ocean_aerosol

    Select Case (model_ocean)
     Case (1)
       brdfmodel_ocean_aerosol = coxmunk_aerosol(theta_obs, theta_sun, phi_obs, phi_sun, phi_del, wind_speed, phiw, I, tau0, g_, ssa_, trans_coeff, phFunc)
     Case Default
       print *, 'wrong BRDF model for ocean'
       stop 
    End Select

  End Function brdfmodel_ocean_aerosol

  ! BRDF kernel model from Roujean, Leroy, and Deschamps (1992)
  Function roujean(theta_obs, phi_del, theta_sun)
    Implicit None
    
    Real, Intent(In)     :: theta_obs, phi_del, theta_sun
    Real, Dimension(0:2) :: roujean
    
    Real, Parameter :: pi           = 3.141592653589793238462643383279502884197
    Real, Parameter :: four_threepi = 4./(3.*pi)
    Real, Parameter :: pi_half      = pi/2.
    Real, Parameter :: one_third    = 1./3.
    Real, Parameter :: one_twopi    = 1./(2.*pi)
    Real, Parameter :: one_pi       = 1./pi
    
    Real :: phi, cos_tobs, cos_tsun, cos_phi, phaseAng, tan_tobs, tan_tsun

    phi = phi_del
    If (phi .lt. 0.) phi = phi + 2.*pi
    If (phi .gt. pi) phi = 2.*pi - phi
  
    ! isotropic term
    roujean(0) = 1.0
    
    ! geometric kernel
    cos_phi  = Cos(phi)
    tan_tobs = Tan(theta_obs)
    tan_tsun = Tan(theta_sun)
    roujean(1) = one_twopi * ( (pi-phi)*cos_phi+Sin(phi) ) *&
         & tan_tobs*tan_tsun -&
         & one_pi * ( tan_tobs + tan_tsun +&
         &   Sqrt( tan_tobs**2 + tan_tsun**2 - 2*tan_tobs*tan_tsun*cos_phi ) )
    
    ! volume-scattering kernel
    cos_tobs = Cos(theta_obs)
    cos_tsun = Cos(theta_sun)
    phaseAng = Acos( cos_tobs*cos_tsun +&
         & Sin(theta_obs)*Sin(theta_sun)*cos_phi )
    roujean(2) = four_threepi * ( (pi_half-phaseAng)*Cos(phaseAng)+Sin(phaseAng) ) /&
         & ( cos_tobs+cos_tsun) - one_third
    
  End Function roujean

  ! BRDF kernel model from Lucht et al. (2000)
 Function rtls(theta_obs, phi_del, theta_sun)
    Implicit None

    Real, Intent(In)     :: theta_obs, phi_del, theta_sun
    Real, Dimension(0:2) :: rtls
  
    Real, Parameter :: pi           = 3.141592653589793238462643383279502884197
    Real, Parameter :: pi_over_two  = pi/2.
    Real, Parameter :: pi_over_four = pi/4.
    Real, Parameter :: one_two      = 1./2.
    Real, Parameter :: one_pi       = 1./pi
  
    Real :: phi, cos_phi, phaseAng, cos_phaseAng, mm, t, cos_t, sin_t
    Real :: cos_tobs, cos_tsun, tan_tobs, tan_tsun
    Real :: big_O
   
    phi = phi_del
    If (phi .lt. 0.) phi = phi + 2.*pi
    If (phi .gt. pi) phi = 2.*pi - phi

    cos_tobs = Cos(theta_obs)
    cos_tsun = Cos(theta_sun)
    cos_phi  = Cos(phi)

    cos_phaseAng = cos_tobs*cos_tsun + Sin(theta_obs)*Sin(theta_sun)*cos_phi
    phaseAng = Acos(cos_phaseAng)
    mm = 1./cos_tobs+1./cos_tsun

    ! isotropic term
    rtls(0) = 1.0

    ! geometric kernel
    tan_tobs = Tan(theta_obs)
    tan_tsun = Tan(theta_sun)
    cos_t = 2./mm *&
         & Sqrt( tan_tobs**2 + tan_tsun**2 - 2*tan_tobs*tan_tsun*cos_phi + &
         &       (tan_tobs*tan_tsun*Sin(phi))**2 )
    cos_t = Minval( (/cos_t, 1./) )
    t = Acos( cos_t )
    sin_t = sin(t)
    big_O = one_pi * mm * (t - sin_t * cos_t)
    rtls(1) = big_O - ((1/cos_tobs)+(1/cos_tsun))  + &
                & one_two * ( 1 + cos_phaseAng ) / (cos_tobs*cos_tsun)
 
   ! volume-scattering kernel
   
    rtls(2) = (((pi_over_two - phaseAng) * cos_phaseAng + Sin(phaseAng)) &
                          & / (cos_tobs+cos_tsun))  &
                          & - pi_over_four

  End Function rtls

  ! BRDF kernel model with hot spot from Maignan et al. (2004)
  Function lirosshotspot(theta_obs, phi_del, theta_sun, xi_0)
    Implicit None
  
    Real, Intent(In)     :: theta_obs, phi_del, theta_sun, xi_0
    Real, Dimension(0:2) :: lirosshotspot
    
    Real, Parameter :: pi           = 3.141592653589793238462643383279502884197
    Real, Parameter :: four_threepi = 4./(3.*pi)
    Real, Parameter :: pi_half      = pi/2.
    Real, Parameter :: one_third    = 1./3.

    Real :: phaseAng0
    Real :: phi, cos_phi, phaseAng, cos_phaseAng, mm, t, cos_t, sin_t
    Real :: cos_tobs, cos_tsun, tan_tobs, tan_tsun
  
    phaseAng0 = xi_0/180.*pi

    ! isotropic term
    lirosshotspot(0) = 1.0
  
    ! geometric kernel
    phi = phi_del
    If (phi .lt. 0.) phi = phi + 2.*pi
    If (phi .gt. pi) phi = 2.*pi - phi 
    cos_phi  = Cos(phi)
    cos_tobs = Cos(theta_obs)
    cos_tsun = Cos(theta_sun)
    tan_tobs = Tan(theta_obs)
    tan_tsun = Tan(theta_sun)
    mm = 1./cos_tobs+1./cos_tsun
    cos_t = 2./mm *&
         & Sqrt( tan_tobs**2 + tan_tsun**2 - 2*tan_tobs*tan_tsun*cos_phi +&
         &       (tan_tobs*tan_tsun*Sin(phi))**2 )
    cos_t = Minval( (/cos_t, 1./) )
    t = Acos( cos_t )
    sin_t = sin(t)
    cos_phaseAng = cos_tobs*cos_tsun + Sin(theta_obs)*Sin(theta_sun)*cos_phi
    phaseAng = Acos(cos_phaseAng)
    lirosshotspot(1) = mm/pi * (t - sin_t*cos_t - pi) +&
         & (1.+cos_phaseAng) / (2*cos_tobs*cos_tsun)

    ! volume-scattering kernel
    lirosshotspot(2) = four_threepi * ( (pi_half-phaseAng)*cos_phaseAng+Sin(phaseAng) ) *&
         & (1.+1./(1.+phaseAng/phaseAng0)) / (cos_tobs+cos_tsun) - one_third
  
  End Function lirosshotspot

  ! BRDF kernel model for ocean surfaces from Cox and Munk (1954)
  ! 1st kernel: Lambertian term (water-leaving reflectance + whitecaps + foam)
  ! 2nd kernel: Fresnel reflection (Sun glint)
  ! 3rd kernel: none
  Function coxmunk(theta_obs, theta_sun, phi_obs, phi_sun, phi_del, wind_speed, phiw, I)
    Implicit None
  
    Real, Intent(In)     :: theta_obs, theta_sun, phi_obs, phi_sun, phi_del, wind_speed, phiw
    Integer, Intent(In)  :: I
    Real, Dimension(0:2) :: coxmunk
    
    Real, Parameter :: pi           = 3.141592653589793238462643383279502884197
  
    Real :: phi, cos_phi
    Real :: mu_n1, cos_theta, sigmac2, sigmau2, sig1
    ! Variables for slope distribution
    Real :: slope_dist, zx, zy, zu, zc, sin0, sinv, cphi0, cphiv, sphi0, sphiv
    ! Variables for shadowing
    Real :: shadow, S1, s3, s2, mu_ir, mu_sq, cot, t1, t2, ShadI, ShadR
    ! Variables for Fresnel reflection coefficient
    Real :: Fresnel_c, mu0, muv, sin_tv, sin_ti, zz, cos_tet, sin_tet, tet0v, st_min=0.001*pi/180., sin_eps=0.00000001, alf, mu_alf
    Complex :: r1, r2, r_par, r_par1, r_par2, r_perp, r_perp1, r_perp2
    Real :: r_par_sq, r_perp_sq
    Complex :: m ! relative complex refractive index of the ocean

    ! 1st kernel
    coxmunk(0) = 1.0

    ! From Table 3 - Sayer et al. 2010  (indexes are cubic-spline-interpolated from AATSR wavelenghts to SEVIRI ones!)
    Select Case (I)
     Case (1) ! from 660nm to 635nm
       m = cmplx(1.33861699,0.0)
     Case (2) ! from 870nm to 810nm
       m = cmplx(1.33495937,0.0)
     Case (3) ! no interpolation, as 1600 is the maximum wavelength in Sayer et al. 2010
       m = cmplx(1.323,0.0)
     Case Default
       print *, 'Stopping because of wrong channel for ocean refractive index?'
       stop
    End Select

    ! wind-related inputs (to become real inputs in the future)
!    phiw = pi/2.0 ! wind direction (in radians); coming from East to represent trade winds (-30 to 30°N)
!    wind_speed = 7.0! in (m/s); 7.0m/s is the annual average over SEVIRI disk from ERA-Interim reanalysis here: http://science.globalwindatlas.info/datasets.html

    ! determining wind-related parameters
    sig1 = (0.003 + 0.00512*wind_speed)/2.0 ! input of GRASP function
    sigmac2 = wind_speed*0.00192+0.003
    sigmau2 = wind_speed*0.00316

    ! 2nd kernel (routine mod_brm.f90 from GRASP open)
    phi = phi_del
    If (phi .lt. 0.) phi = phi + 2.*pi
    If (phi .gt. pi) phi = 2.*pi - phi
    cos_phi  = Cos(phi)
    muv = Cos(theta_obs)
    mu0 = Cos(theta_sun)
    sinv = Sin(theta_obs)
    sin0 = Sin(theta_sun)
    cos_theta = -mu0*muv - sin0*sinv*cos_phi
    mu_n1 = (muv + mu0)/(2.0*sin(0.5*acos(cos_theta)))

    ! Calculation of slope distribution (Cox and Munck, 1954)
!    cphi0 =-cos(phi_sun)
!    cphiv = cos(phi_obs)
!    sphi0 =-sin(phi_sun)
!    sphiv = sin(phi_obs)
!    zx    = (sinv*cphiv-sin0*cphi0)/(mu0+muv)
!    zy    = (sinv*sphiv-sin0*sphi0)/(mu0+muv)
    zx    = (-sinv*sin(phi))/(mu0+muv)
    zy    = (sin0+sinv*cos_phi)/(mu0+muv)
    zu    = cos(phiw)*zx+sin(phiw)*zy
    zc    =-sin(phiw)*zx+cos(phiw)*zy
    slope_dist = exp(-0.5*(zc*zc/sigmac2+zu*zu/sigmau2))
    slope_dist = slope_dist/(2.0*pi*sqrt(sigmac2*sigmau2)*mu_n1**3)

    ! Calculation of Fresnel reflection coefficient for Gaussian surface
    sin_tv=sqrt(1.0-muv*muv)
    sin_ti=sqrt(1.0-mu0*mu0)
    zz=muv*mu0+sin_ti*sin_tv*cos_phi
    cos_tet=-zz
    sin_tet=sqrt(1.0-cos_tet*cos_tet)
    If (abs(sin_tv) .lt. sin_eps) then
       tet0v=acos(muv)
       tet0v=abs(tet0v-st_min)
       muv=cos(tet0v)
       sin_tv=sqrt(1.0-muv*muv)
       zz=muv*mu0+sin_ti*sin_tv*cos_phi
    end if
    If (abs(sin_ti) .lt. sin_eps) then
       tet0v=acos(mu0)
       tet0v=abs(tet0v-st_min)
       mu0=cos(tet0v)
       sin_ti=sqrt(1.0-mu0*mu0)
       zz=muv*mu0+sin_ti*sin_tv*cos_phi
    end if
    if (abs(sin_tet) .lt. sin_eps) then
       tet0v=acos(mu0)
       tet0v=abs(tet0v-st_min)
       muv=cos(tet0v)
       sin_tv=sqrt(1.0-muv*muv)
       zz=muv*mu0+sin_ti*sin_tv*cos_phi
    end if
    alf=acos(zz)
    mu_alf=cos(alf*0.5)
    r1=m*m*mu_alf
    r2=csqrt(m*m-1.+mu_alf*mu_alf)
    r_par1=r1-r2
    r_par2=r1+r2
    r_perp1=mu_alf-r2
    r_perp2=mu_alf+r2
    r_par=r_par1/r_par2
    r_perp=r_perp1/r_perp2
    r_par_sq=r_par*conjg(r_par)
    r_perp_sq=r_perp*conjg(r_perp)
    Fresnel_c=(r_par_sq+r_perp_sq)/2.0

    ! Shadowing function (Mischenko and Travis, 1997 ?)
    S1=SQRT(2.0*sig1/pi)
    S3=1.0/(SQRT(2.0*sig1))
    S2=S3*S3
    mu_ir=mu0
    mu_sq=mu_ir*mu_ir
    cot=mu_ir/SQRT(1.0-mu_sq)
    T1=EXP(-cot*cot*S2)
    T2=ERFC(cot*S3)
    ShadI=0.5*(S1*T1/cot-T2)
    mu_ir=muv
    mu_sq=mu_ir*mu_ir
    cot=mu_ir/SQRT(1.0-mu_sq)
    T1=EXP(-cot*cot*S2)
    T2=ERFC(cot*S3)
    ShadR=0.5*(S1*T1/cot-T2)
    shadow=1.0/(1.0+ShadI+ShadR)

    coxmunk(1) = pi*slope_dist*Fresnel_c*shadow/4.0/muv/mu0/mu_n1

    ! 3rd kernel
    coxmunk(2) = 0.0
  
  End Function coxmunk

  ! BRDF kernel model from Roujean, Leroy, and Deschamps (1992) plus an aerosol kernel using a double Henyey-Greenstein phase function
  Function roujean_aerosol(theta_obs, phi_del, theta_sun, tau0, g_, ssa_, trans_coeff, phFunc)
    Implicit None
    
    Real, Intent(In)  :: theta_obs, phi_del, theta_sun
    Real, Intent(In)  :: g_, ssa_, trans_coeff, tau0, phFunc

    Real, Dimension(0:3) :: roujean_aerosol
    
    Real, Parameter :: pi           = 3.141592653589793238462643383279502884197
    Real, Parameter :: four_threepi = 4./(3.*pi)
    Real, Parameter :: pi_half      = pi/2.
    Real, Parameter :: one_third    = 1./3.
    Real, Parameter :: one_twopi    = 1./(2.*pi)
    Real, Parameter :: one_pi       = 1./pi
    
    Real :: phi, cos_tobs, cos_tsun, cos_phi, phaseAng, tan_tobs, tan_tsun, p, exp_x, nu

    phi = phi_del
    If (phi .lt. 0.) phi = phi + 2.*pi
    If (phi .gt. pi) phi = 2.*pi - phi

    ! isotropic term
    roujean_aerosol(0) = 1.0
    
    ! geometric kernel
    cos_phi  = Cos(phi)
    tan_tobs = Tan(theta_obs)
    tan_tsun = Tan(theta_sun)
    roujean_aerosol(1) =1.*( one_twopi * ( (pi-phi)*cos_phi+Sin(phi) ) *&
         & tan_tobs*tan_tsun -&
         & one_pi * ( tan_tobs + tan_tsun +&
         &   Sqrt( tan_tobs**2 + tan_tsun**2 - 2*tan_tobs*tan_tsun*cos_phi ) ))
    
    ! volume-scattering kernel
    cos_tobs = Cos(theta_obs)
    cos_tsun = Cos(theta_sun)
    phaseAng = Acos( cos_tobs*cos_tsun +&
         & Sin(theta_obs)*Sin(theta_sun)*cos_phi )
    roujean_aerosol(2) = four_threepi * ( (pi_half-phaseAng)*Cos(phaseAng)+Sin(phaseAng) ) /&
         & ( cos_tobs+cos_tsun) - one_third
    
    !aerosol kernel
!    p=(1.-g_**2)/(1.+g_**2-2.*g_*cos(pi-phaseAng))**1.5 + (1.-g_**2)*f_2HG*(3.*cos(pi-phaseAng)**2.0-1.0)/(2.0*(1.+g_**2.)**1.5)
    p=phFunc
    !p=p/(1.-eta) ! truncation
    nu=1./cos_tobs+1./cos_tsun
    exp_x=(840.-60.*(tau0*nu)+20.*(tau0*nu)**2-(tau0*nu)**3)/(840.+360.*(tau0*nu)+60.*(tau0*nu)**2+4*(tau0*nu)**3)
!    roujean_aerosol(3) =ssa_/4./cos_tsun/cos_tobs*p*exp_x*(1.+(7.-tau0)*tau0/5.)
    roujean_aerosol(3) =ssa_/4./cos_tsun/cos_tobs*p*exp_x

    if (abs(roujean_aerosol(3)) .le. 1e-5) then
            roujean_aerosol(3) = 0.
    endif
    
    roujean_aerosol(3) = roujean_aerosol(3) / trans_coeff
    
  End Function roujean_aerosol 

  ! BRDF kernel model from Lutch et al. (2000) plus an aerosol kernel using a double Henyey-Greenstein phase function
  Function rtls_aerosol(theta_obs, phi_del, theta_sun, tau0, g_, ssa_, trans_coeff, phFunc)
    Implicit None
  
    Real, Intent(In)     :: theta_obs, phi_del, theta_sun
    Real, Dimension(0:3) :: rtls_aerosol
    Real, Intent(In)  :: g_, ssa_, trans_coeff, tau0, phFunc

    Real, Parameter :: pi           = 3.141592653589793238462643383279502884197
    Real, Parameter :: one_two      = 1./2.
    Real, Parameter :: one_pi       = 1./pi
    Real, Parameter :: pi_over_four = pi/4.
    Real, Parameter :: pi_over_two = pi/2.

    Real :: big_O
    Real :: phi, cos_phi, phaseAng, cos_phaseAng, mm, t, cos_t, sin_t
    Real :: cos_tobs, cos_tsun, tan_tobs, tan_tsun
    Real :: p, exp_x, nu
    
    phi = phi_del
    If (phi .lt. 0.) phi = phi + 2.*pi
    If (phi .gt. pi) phi = 2.*pi - phi

    cos_tobs = Cos(theta_obs)
    cos_tsun = Cos(theta_sun)
    cos_phi  = Cos(phi)

    cos_phaseAng = cos_tobs*cos_tsun + Sin(theta_obs)*Sin(theta_sun)*cos_phi
    phaseAng = Acos(cos_phaseAng)
    mm = 1./cos_tobs+1./cos_tsun

    ! isotropic term
    rtls_aerosol(0) = 1.0

    ! geometric kernel
    tan_tobs = Tan(theta_obs)
    tan_tsun = Tan(theta_sun)
    cos_t = 2./mm *&
         & Sqrt( tan_tobs**2 + tan_tsun**2 - 2*tan_tobs*tan_tsun*cos_phi + &
         &       (tan_tobs*tan_tsun*Sin(phi))**2 )
    cos_t = Minval( (/cos_t, 1./) )
    t = Acos( cos_t )
    sin_t = sin(t)
    big_O = one_pi * mm * (t - sin_t * cos_t)
    rtls_aerosol(1) = big_O - ((1/cos_tobs)+(1/cos_tsun))  + &
                & one_two * ( 1 + cos_phaseAng ) / (cos_tobs*cos_tsun)
 
   ! volume-scattering kernel
   
    rtls_aerosol(2) = (((pi_over_two - phaseAng) * cos_phaseAng + Sin(phaseAng)) &
                          & / (cos_tobs+cos_tsun))  &
                          & - pi_over_four
  
   ! aerosol kernel
!    p=(1.-g_**2)/(1.+g_**2-2.*g_*cos(pi-phaseAng))**1.5 + (1.-g_**2)*f_2HG*(3.*cos(pi-phaseAng)**2.0-1.0)/(2.0*(1.+g_**2.)**1.5)
    p=phFunc
    !p=p/(1.-eta) ! truncation
    nu=1./cos_tobs+1./cos_tsun
    exp_x=(840.-60.*(tau0*nu)+20.*(tau0*nu)**2-(tau0*nu)**3)/(840.+360.*(tau0*nu)+60.*(tau0*nu)**2+4*(tau0*nu)**3)
!    rtls_aerosol(3) =ssa_/4./cos_tsun/cos_tobs*p*exp_x*(1.+(7.-tau0)*tau0/5.)
    rtls_aerosol(3) =ssa_/4./cos_tsun/cos_tobs*p*exp_x

    if (abs(rtls_aerosol(3)) .le. 1e-5) then
            rtls_aerosol(3) = 0.
    endif
  
    rtls_aerosol(3) = rtls_aerosol(3) / trans_coeff
    
  End Function rtls_aerosol

  ! BRDF kernel model with hot spot from Maignan et al. (2004) plus an aerosol kernel using a double Henyey-Greenstein phase function
  Function lirosshotspot_aerosol(theta_obs, phi_del, theta_sun, tau0, g_, ssa_, trans_coeff, phFunc, xi_0)
    Implicit None
  
    Real, Intent(In)     :: theta_obs, phi_del, theta_sun, xi_0
    Real, Dimension(0:3) :: lirosshotspot_aerosol
    Real, Intent(In)  :: g_, ssa_, trans_coeff, tau0, phFunc
    
    Real, Parameter :: pi           = 3.141592653589793238462643383279502884197
    Real, Parameter :: four_threepi = 4./(3.*pi)
    Real, Parameter :: pi_half      = pi/2.
    Real, Parameter :: one_third    = 1./3.

    Real :: phaseAng0
    Real :: phi, cos_phi, phaseAng, cos_phaseAng, mm, t, cos_t, sin_t
    Real :: cos_tobs, cos_tsun, tan_tobs, tan_tsun

    Real :: p, exp_x, nu

    phaseAng0 = xi_0/180.*pi

    ! isotropic term
    lirosshotspot_aerosol(0) = 1.0
  
    ! geometric kernel
    phi = phi_del
    If (phi .lt. 0.) phi = phi + 2.*pi
    If (phi .gt. pi) phi = 2.*pi - phi 
    cos_phi  = Cos(phi)
    cos_tobs = Cos(theta_obs)
    cos_tsun = Cos(theta_sun)
    tan_tobs = Tan(theta_obs)
    tan_tsun = Tan(theta_sun)
    mm = 1./cos_tobs+1./cos_tsun
    cos_t = 2./mm *&
         & Sqrt( tan_tobs**2 + tan_tsun**2 - 2*tan_tobs*tan_tsun*cos_phi +&
         &       (tan_tobs*tan_tsun*Sin(phi))**2 )
    cos_t = Minval( (/cos_t, 1./) )
    t = Acos( cos_t )
    sin_t = sin(t)
    cos_phaseAng = cos_tobs*cos_tsun + Sin(theta_obs)*Sin(theta_sun)*cos_phi
    phaseAng = Acos(cos_phaseAng)
    lirosshotspot_aerosol(1) = mm/pi * (t - sin_t*cos_t - pi) +&
         & (1.+cos_phaseAng) / (2*cos_tobs*cos_tsun)

    ! volume-scattering kernel
    lirosshotspot_aerosol(2) = four_threepi * ( (pi_half-phaseAng)*cos_phaseAng+Sin(phaseAng) ) *&
         & (1.+1./(1.+phaseAng/phaseAng0)) / (cos_tobs+cos_tsun) - one_third

    !aerosol kernel
!    p=(1.-g_**2)/(1.+g_**2-2.*g_*cos(pi-phaseAng))**1.5 + (1.-g_**2)*f_2HG*(3.*cos(pi-phaseAng)**2.0-1.0)/(2.0*(1.+g_**2.)**1.5)
    p=phFunc
    !p=p/(1.-eta) ! truncation
    nu=1./cos_tobs+1./cos_tsun
    exp_x=(840.-60.*(tau0*nu)+20.*(tau0*nu)**2-(tau0*nu)**3)/(840.+360.*(tau0*nu)+60.*(tau0*nu)**2+4*(tau0*nu)**3)
!    lirosshotspot_aerosol(3) =ssa_/4./cos_tsun/cos_tobs*p*exp_x*(1.+(7.-tau0)*tau0/5.)
    lirosshotspot_aerosol(3) =ssa_/4./cos_tsun/cos_tobs*p*exp_x

    if (abs(lirosshotspot_aerosol(3)) .le. 1e-5) then
            lirosshotspot_aerosol(3) = 0.
    endif
    
    lirosshotspot_aerosol(3) = lirosshotspot_aerosol(3) / trans_coeff
  
  End Function lirosshotspot_aerosol

  ! BRDF kernel model for ocean surfaces from Cox and Munk (1954) plus the aerosol kernel
  ! 1st kernel: Lambertian term (water-leaving reflectance + whitecaps + foam)
  ! 2nd kernel: Fresnel reflection (Sun glint)
  ! 3rd kernel: none
  ! 4th kernel: aerosols
  Function coxmunk_aerosol(theta_obs, theta_sun, phi_obs, phi_sun, phi_del, wind_speed, phiw, I, tau0, g_, ssa_, trans_coeff, phFunc)
    Implicit None
  
    Real, Intent(In)     :: theta_obs, theta_sun, phi_obs, phi_sun, phi_del, wind_speed, phiw
    Integer, Intent(In)  :: I
    Real, Dimension(0:3) :: coxmunk_aerosol
    Real, Intent(In)  :: g_, ssa_, trans_coeff, tau0, phFunc
    
    Real, Parameter :: pi           = 3.141592653589793238462643383279502884197
  
    Real :: phi, cos_phi
    Real :: mu_n1, cos_theta, sigmac2, sigmau2, sig1
    ! Variables for slope distribution
    Real :: slope_dist, zx, zy, zu, zc, sin0, sinv, cphi0, cphiv, sphi0, sphiv
    ! Variables for shadowing
    Real :: shadow, S1, s3, s2, mu_ir, mu_sq, cot, t1, t2, ShadI, ShadR
    ! Variables for Fresnel reflection coefficient
    Real :: Fresnel_c, mu0, muv, sin_tv, sin_ti, zz, cos_tet, sin_tet, tet0v, st_min=0.001*pi/180., sin_eps=0.00000001, alf, mu_alf
    Complex :: r1, r2, r_par, r_par1, r_par2, r_perp, r_perp1, r_perp2
    Real :: r_par_sq, r_perp_sq
    Complex :: m ! relative complex refractive index of the ocean
    ! Variables for aerosol kernel
    Real :: p, exp_x, nu

    ! 1st kernel
    coxmunk_aerosol(0) = 1.0
  
    ! From Table 3 - Sayer et al. 2010  (indexes are cubic-spline-interpolated from AATSR wavelenghts to SEVIRI ones!)
    Select Case (I)
     Case (1) ! from 660nm to 635nm
       m = cmplx(1.33861699,0.0)
     Case (2) ! from 870nm to 810nm
       m = cmplx(1.33495937,0.0)
     Case (3) ! no interpolation, as 1600 is the maximum wavelength in Sayer et al. 2010
       m = cmplx(1.323,0.0)
     Case Default
       print *, 'Stopping because of wrong channel for ocean refractive index?'
       stop
    End Select

    ! wind-related inputs (to become real inputs in the future)
!    phiw = pi/4.0 ! wind direction (in radians); coming from East to represent trade winds (-30 to 30°N)
!    wind_speed = 7.0! in (m/s); 7.0m/s is the annual average over SEVIRI disk from ERA-Interim reanalysis here: http://science.globalwindatlas.info/datasets.html

    ! determining wind-related parameters
    sig1 = (0.003 + 0.00512*wind_speed)/2.0 ! input of GRASP function
    sigmac2 = wind_speed*0.00192+0.003
    sigmau2 = wind_speed*0.00316

    ! 2nd kernel (routine mod_brm.f90 from GRASP open)
    phi = phi_del
    If (phi .lt. 0.) phi = phi + 2.*pi
    If (phi .gt. pi) phi = 2.*pi - phi
    cos_phi  = Cos(phi)
    muv = Cos(theta_obs)
    mu0 = Cos(theta_sun)
    sinv = Sin(theta_obs)
    sin0 = Sin(theta_sun)
    cos_theta = -mu0*muv - sin0*sinv*cos_phi
    mu_n1 = (muv + mu0)/(2.0*sin(0.5*acos(cos_theta)))

    ! Calculation of slope distribution (Cox and Munck, 1954)
!    cphi0 =-cos(phi_sun)
!    cphiv = cos(phi_obs)
!    sphi0 =-sin(phi_sun)
!    sphiv = sin(phi_obs)
!    zx    = (sinv*cphiv-sin0*cphi0)/(mu0+muv)
!    zy    = (sinv*sphiv-sin0*sphi0)/(mu0+muv)
    zx    = (-sinv*sin(phi))/(mu0+muv)
    zy    = (sin0+sinv*cos_phi)/(mu0+muv)
    zu    = cos(phiw)*zx+sin(phiw)*zy
    zc    =-sin(phiw)*zx+cos(phiw)*zy
    slope_dist = exp(-0.5*(zc*zc/sigmac2+zu*zu/sigmau2))
    slope_dist = slope_dist/(2.0*pi*sqrt(sigmac2*sigmau2)*mu_n1**3)

!    print*, 'slope_dist', slope_dist
!    print*, 'exp(-0.5*(zc*zc/sigmac2+zu*zu/sigmau2))', exp(-0.5*(zc*zc/sigmac2+zu*zu/sigmau2))
!    print*, '-0.5*(zc*zc/sigmac2+zu*zu/sigmau2)', -0.5*(zc*zc/sigmac2+zu*zu/sigmau2)
!    print*, 'zc*zc/sigmac2', zc*zc/sigmac2
!    print*, 'zu*zu/sigmau2', zu*zu/sigmau2
!    print*, 'zc, zu, zx, zy', zc, zu, zx, zy
!    print*, 'sigmac2, sigmau2', sigmac2, sigmau2
!    print*, 'theta_obs, theta_sun, phi_obs, phi_sun, phi', theta_obs*180./pi, theta_sun*180./pi, phi_obs*180./pi, phi_sun*180./pi, phi*180./pi
!    print*, 'muv, mu0, sinv, sin0, cos_phi, cos_theta, mu_n1', muv, mu0, sinv, sin0, cos_phi, cos_theta, mu_n1
!!    print*, 'cphi0, cphiv, sphi0, sphiv', cphi0, cphiv, sphi0, sphiv 
!    print*, '(2.0*pi*sqrt(sigmac2*sigmau2)*mu_n1**3)', (2.0*pi*sqrt(sigmac2*sigmau2)*mu_n1**3)
!    print*, 'mu_n1', mu_n1

    ! Calculation of Fresnel reflection coefficient for Gaussian surface
    sin_tv=sqrt(1.0-muv*muv)
    sin_ti=sqrt(1.0-mu0*mu0)
    zz=muv*mu0+sin_ti*sin_tv*cos_phi
    cos_tet=-zz
    sin_tet=sqrt(1.0-cos_tet*cos_tet)
    If (abs(sin_tv) .lt. sin_eps) then
       tet0v=acos(muv)
       tet0v=abs(tet0v-st_min)
       muv=cos(tet0v)
       sin_tv=sqrt(1.0-muv*muv)
       zz=muv*mu0+sin_ti*sin_tv*cos_phi
    end if
    If (abs(sin_ti) .lt. sin_eps) then
       tet0v=acos(mu0)
       tet0v=abs(tet0v-st_min)
       mu0=cos(tet0v)
       sin_ti=sqrt(1.0-mu0*mu0)
       zz=muv*mu0+sin_ti*sin_tv*cos_phi
    end if
    if (abs(sin_tet) .lt. sin_eps) then
       tet0v=acos(mu0)
       tet0v=abs(tet0v-st_min)
       muv=cos(tet0v)
       sin_tv=sqrt(1.0-muv*muv)
       zz=muv*mu0+sin_ti*sin_tv*cos_phi
    end if
    alf=acos(zz)
    mu_alf=cos(alf*0.5)
    r1=m*m*mu_alf
    r2=csqrt(m*m-1.+mu_alf*mu_alf)
    r_par1=r1-r2
    r_par2=r1+r2
    r_perp1=mu_alf-r2
    r_perp2=mu_alf+r2
    r_par=r_par1/r_par2
    r_perp=r_perp1/r_perp2
    r_par_sq=r_par*conjg(r_par)
    r_perp_sq=r_perp*conjg(r_perp)
    Fresnel_c=(r_par_sq+r_perp_sq)/2.0

    ! Shadowing function (Mischenko and Travis, 1997 ?)
    S1=SQRT(2.0*sig1/pi)
    S3=1.0/(SQRT(2.0*sig1))
    S2=S3*S3
    mu_ir=mu0
    mu_sq=mu_ir*mu_ir
    cot=mu_ir/SQRT(1.0-mu_sq)
    T1=EXP(-cot*cot*S2)
    T2=ERFC(cot*S3)
    ShadI=0.5*(S1*T1/cot-T2)
    mu_ir=muv
    mu_sq=mu_ir*mu_ir
    cot=mu_ir/SQRT(1.0-mu_sq)
    T1=EXP(-cot*cot*S2)
    T2=ERFC(cot*S3)
    ShadR=0.5*(S1*T1/cot-T2)
    shadow=1.0/(1.0+ShadI+ShadR)

    coxmunk_aerosol(1) = pi*slope_dist*Fresnel_c*shadow/4.0/muv/mu0/mu_n1

    ! 3rd kernel
    coxmunk_aerosol(2) = 0.0

    ! 4th (and aerosol) kernel
    p=phFunc
    !p=p/(1.-eta) ! truncation
    nu=1./muv+1./mu0
    exp_x=(840.-60.*(tau0*nu)+20.*(tau0*nu)**2-(tau0*nu)**3)/(840.+360.*(tau0*nu)+60.*(tau0*nu)**2+4*(tau0*nu)**3)
    coxmunk_aerosol(3) =ssa_/4./mu0/muv*p*exp_x
    if (abs(coxmunk_aerosol(3)) .le. 1e-5) coxmunk_aerosol(3) = 0.
    coxmunk_aerosol(3) = coxmunk_aerosol(3) / trans_coeff
  
  End Function coxmunk_aerosol

!-------------------------------------------------------

  ! select look-up-table
  Subroutine brdfinit_land(model_land)
    Implicit None

    Integer, Intent(In) :: model_land

    ! Roujean, Leroy, and Deschamps (1992)
    Real, Dimension(0:n_table-1,1:3), Parameter :: hint_roujean = Reshape( &
         (/ -0.997910  , -0.998980  , -1.00197   , -1.00702 ,   -1.01438    , &
         -1.02443   , -1.03773   , -1.05501   , -1.07742 ,   -1.10665    , &
         -1.14526   , -1.19740   , -1.27008   , -1.37595 ,   -1.54059    , &
         -1.82419   , -2.40820   , -4.20369   , &
         -0.00894619, -0.00837790, -0.00665391, -0.00371872,  0.000524714, &
         0.00621877,  0.0135606 ,  0.0228129 ,  0.0343240 ,  0.0485505  , &
         0.0661051 ,  0.0878086 ,  0.114795  ,  0.148698  ,  0.191944   , &
         0.248471  ,  0.325351  ,  0.438371, &
         0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. , 0., 0., 0., 0./), (/ n_table, 3 /) )
    Real, Parameter, Dimension(0:3) :: bihi_calc_roujean = (/ 1.,-1.28159,8.02838E-02 , 0./)
    
    ! "LiRossHotspot" Maignan et al. (2004)
    Real, Dimension(0:n_table-1,1:3), Parameter :: hint_hotspot = Reshape( &
         (/ -1.2872    , -1.2883    , -1.29142   , -1.2966  ,   -1.30384    , &
         -1.31307   , -1.3243    , -1.33744   , -1.35237 ,   -1.3689     , &
         -1.38686   , -1.40582   , -1.4253    , -1.44471 ,   -1.46328    , &
         -1.48025   , -1.49538   , -1.51218   , &
         0.0052371 ,  0.00581059,  0.00754731,  0.0105049 ,  0.0147809  , &
         0.0205190 ,  0.0279187 ,  0.0372467 ,  0.0488525 ,  0.0632023  , &
         0.0809112 ,  0.102814  ,  0.130061  ,  0.164302  ,  0.208015   , &
         0.265199  ,  0.343066  ,  0.457749, &
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. , 0., 0., 0., 0. /), (/ n_table, 3 /) )
    Real, Parameter, Dimension(0:3) :: bihi_calc_hotspot = (/ 1.,-1.37760,9.52950E-02, 0. /)

   ! RTLS
    Real, Dimension(0:n_table-1, 1:3), Parameter :: hint_rtls=Reshape( &
         (/-1.2889 , -1.2899 , -1.293 , -1.2981 , -1.3053 , &
          -1.3145  , -1.3256 , -1.3387, -1.3535 , -1.3698 , &
          -1.3875  , -1.4062 , -1.4253, -1.4441, -1.4618 , &
          -1.4773  , -1.4895 , -1.4973, &
	  -0.02107921, -0.01973968, -0.01567785, -0.00876165, 0.00123677, &
           0.01465346, 0.03195199, 0.05375256, 0.08087403, 0.11439660, &
	   0.15575623, 0.20689142, 0.27048166, 0.35035791, 0.45226682, &
           0.58545999, 0.76661237, 1.03292777, &
           0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. , 0., 0., 0., 0. /), (/ n_table, 3 /) )
    Real, Parameter, Dimension(0:3) :: bihi_calc_rtls = (/ 1.,-1.37751,1.89186E-01, 0. /)
     
    Select Case (model_land)
     Case (1)
       print *, 'BRDF model for land is Roujean'  
       hint_land = hint_roujean
       bihi_calc_land = bihi_calc_roujean
     Case (2)
       print *, 'BRDF model for land is RTLS'  
       hint_land = hint_rtls
       bihi_calc_land = bihi_calc_rtls
     Case (3)
       print *, 'BRDF model for land is Maignan'  
       hint_land = hint_hotspot
       bihi_calc_land = bihi_calc_hotspot
     Case Default
       print *, 'Wrong BRDF model for land'
       stop 
    End Select

  End Subroutine brdfinit_land

  Subroutine brdfinit_ocean(model_ocean)
    Implicit None

    Integer, Intent(In) :: model_ocean

    ! Cox and Munk (195X?)
    Real, Dimension(0:n_table-1,1:3), Parameter :: hint_coxmunk = Reshape( & ! -> TBD!!!
         (/ -0.997910  , -0.998980  , -1.00197   , -1.00702 ,   -1.01438 , & ! -> TBD!!!
         -1.02443   , -1.03773   , -1.05501   , -1.07742 ,   -1.10665    , & ! -> TBD!!!
         -1.14526   , -1.19740   , -1.27008   , -1.37595 ,   -1.54059    , & ! -> TBD!!!
         -1.82419   , -2.40820   , -4.20369   , &                            ! -> TBD!!!
         0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. , 0., 0., 0., 0., &
         0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. , 0., 0., 0., 0./), (/ n_table, 3 /) )
    Real, Parameter, Dimension(0:3) :: bihi_calc_coxmunk = (/ 1.,-1.28159, 0. , 0./) ! -> TBD!!!
    
    Select Case (model_ocean)
     Case (1)
       print *, 'BRDF model for ocean is Cox and Munk'  
       hint_ocean = hint_coxmunk
       bihi_calc_ocean = bihi_calc_coxmunk
     Case Default
       print *, 'Wrong BRDF model for ocean'
       stop 
    End Select

  End Subroutine brdfinit_ocean

End Module brdfmodels
