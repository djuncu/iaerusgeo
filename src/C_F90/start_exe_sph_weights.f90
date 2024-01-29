Subroutine start_exe(ProcessingStatus_p)
   Use nr, Only: svdvar
   Use py_ifc
   Use algoconf
   Use start_mod
   Use brdfmodels
   Use aerosol_parameters

   Implicit None

   Integer(Kind=taokind) tmptaufrc

   Integer, Intent(Out) :: ProcessingStatus_p

   External reportLog, getStopStatus, stopping

   ! number of lines to be read
   Integer :: Lines

   ! character variable definitions
   Character(LEN=4) :: linestr
   Character(LEN=7) :: status

   ! variables for date_and_time subroutine
   Integer :: day_of_year ! day of the year

   Real, Dimension(:, :), Allocatable :: F3_res

   ! days since the last execution of the algorithm
   Real    :: days_last_in
   ! relative increase of a priori covariance matrix
   Real    :: delta

   Real    :: alpha_intrp ! auxiliary variable for linear interpolation
   Real    :: t_sat_lim ! satellite zenith angle limit (in radians)
   Real    :: t_sol_lim ! solar zenith angle limit (in radians)
   Real    :: t_sat_wlim ! limit for weighting equation (in radians)
   Real    :: t_sol_wlim ! limit for weighting equation (in radians)
   Real    :: x_lim ! limit for scattering angle (in radians)
   Real    :: lat ! geogr. latitude of the processed pixel
   Real    :: lon ! geogr. longitude of the processed pixel
   Real    :: theta_ref_dh ! reference angle for dir.-hem. albedo
   Real    :: theta_sol_midi ! solar zenith angle at local solar noon
   Integer :: age ! "age" of last observation used
   Integer :: N_valid_obs ! number of valid observations for DAILY run
   Integer :: N_valid_snow ! number of valid observations with snow flag set
   Integer :: N_snow_limit ! minimum number for setting the snow flag
   Integer :: model ! model for land or ocean, depending on the pixel
   Real    :: albed ! auxiliary variable for albedo
   Real    :: sigma ! auxiliary variable for albedo error
   Real    :: variance ! auxiliary variable for variance
   Logical :: observations_d = .false. ! .true. if a sufficient number of observations are available for DAILY
   Logical :: observations_i = .false. ! .true. if a sufficient number of observations are available for INST
   Logical :: previous ! .true. if previous estimate is available
   Logical :: snow ! .true. if pixel is declared as snow covered
   Logical :: convergence_OK, limit_AOD, persistent_low_AOD, persistent_high_AOD

   Real, Dimension(0:MM)         :: k_reg
   Real, Dimension(0:MM)         :: sig_k_reg

   Integer(Kind=maskind) :: lw_mask ! auxiliary variable for land/water mask

   Real    :: eta_sph=100./6371.
   Integer :: N, NB, NN, X, Y, YY, I, J, IHF, T, S, KF ! various loop indices
   Real    :: tau, tau_0, albed_b_in, albed_b_in_, albed_d_in, refl_mean, tau_1, tau_badValue = -1.0
   Real    :: delta_tau=0.01, albed_force, albed_force_land=0.15, albed_force_ocean=0.00, tau_force, Sa_val=0.05, bright_alb = 0.20, k0_ocean_max=1.0
   Integer :: aod_pos, aod_pos_j
   Real    :: tau_0_tilde, tau_0_tilde_j, x1_tilde
   Real, Dimension(1:N_AOD_max) :: phFunc
   Real, Dimension(1:MaxNHistFiles) :: tau_in_hist
   Real    :: TS, TV, Salb, trans_coeff, rho_1, R_MS, R_MS_us, R_MS_uv, us, uv, xi
   Integer :: cpt_tau, cpt_tau_conv, cpt_tau_max, counter

   Integer :: aer_type_1 = -1, aer_type_2 = -1, aer_type_3 = -1, aer_type_prev_1 = -1
   Logical :: flag_aer_type_2, flag_aer_type_3
   Real    :: weight_mod_1, weight_mod_2, weight_mod_3

   Logical :: freezeDaily = .false. ! update of surface BRDF is disabled
   Logical :: smooth_flag = .true. ! spatiotemporal smoothing
   Logical :: debug_flag, ocean_flag, land_flag, bright_surface, tau_in_hist_flag
   Logical :: flag_failedInvExit, flag_invClimato, flag_invFailed, flag_kinIsClimato = .false.
   Real    :: tau_in_hist_mean, tau_in_hist_std, tau_in_hist_flag_val=0.03
   Logical :: writeInversion_flag = .false.

   Logical :: read_tmp_file, write_tmp_file, IsError
   Integer :: ipct = 0, BlockRest = 0

   Integer                                                   :: halfBoxSize = (boxSize+1)/2
   Integer, Dimension(0:boxSize*boxSize-1)                   :: X_X, Y_Y ! nearest pixels
   Integer                                                   :: cptX, cptX_final, NS, i_smooth
   Real, Dimension(0:boxSize*boxSize*N_scenes-1, 0:MM)       :: k_in_near, k_in_near_ ! a priori information for parameters
   Real, Dimension(0:boxSize*boxSize*N_scenes-1)             :: uv_near, us_near, ref_tol_near
   Real, Dimension(0:boxSize*boxSize*N_scenes-1)             :: albed_b_in_near, albed_b_in_near_, ref_s_near, ref_s_near_
   Real, Dimension(0:boxSize*boxSize*N_scenes-1)             :: weight_near, xi_near
   Real, Dimension(0:boxSize*boxSize*N_scenes-1,1:N_AOD_max) :: phFunc_near
   Integer, Dimension(0:boxSize*boxSize*N_scenes-1)          :: time_near
   Real                                                      :: theta_sol_near, theta_sat_near, phi_del_near, phi_sol_near, phi_sat_near, sunglint_near
   Real                                                      :: ref_tol_near_mean, ref_tol_near_std, weight_near_sum

   Real :: ref_s, ref_s_, w_ref_s, w_ref_s_, ref_tol, ref_tol_j, ref_tol_out, fit_error, fit_error_force, jacobian
   Real, Dimension(1:N_AOD_max) :: P_tilde, P_tilde_j

   Real :: jacoAOD, jacoAOD_thrs=0.2, conf_Meas = -1

   ! interface of solar zenith angle subroutine
   Interface
      Subroutine solzenith(day_of_year_ref, latitude, theta_sol_midi)
         Implicit None
         Integer, Intent(In) :: day_of_year_ref
         Real, Intent(In)    :: latitude
         Real, Intent(Out)   :: theta_sol_midi
      End Subroutine solzenith
   End Interface

   ! interface for Julian Date subroutine
   Interface
      Function julday(mm, id, iyyy)
         Integer, Intent(In) :: mm, id, iyyy
         Integer  :: julday
      End Function julday
   End Interface

   ! interface of instantaneous retrieval
   Interface
      Subroutine LM(Sa, g_tilde, ssa, ssa_tilde, eta, uv, us, phFunc, refl_, refl_out_, ref_s, albed_b_in, tau_ap, tau0, jacoAOD, AOD_dynamic_steps, debug_flag)
         Use py_ifc
         Implicit None
         Real, Dimension(1:N_AOD_max), Intent(In)     :: g_tilde, ssa, ssa_tilde, eta, phFunc ! Aerosol parameters
         Real, Intent(In)                             :: refl_, ref_s ! Reflectances
         Real, Intent(In)                             :: Sa, tau_ap
         Real, Intent(Out)                            :: refl_out_
         Real, Intent(Out)                            :: tau0 ! AOD
         Real, Intent(Out)                            :: jacoAOD ! Jacobian
         Real, Intent(In)                             :: albed_b_in
         Real, Intent(In)                             :: uv, us
         Logical, Intent(In)                          :: debug_flag
         Real, Dimension(1:N_AOD_max), Intent(In)     :: AOD_dynamic_steps
      End Subroutine LM
   End Interface

   ! initialise getStopStatus counter
   Call System_Clock(Count_Old, Count_Rate, Count_Max)
   Count_Wait = Count_Rate*YFreqInSecs
   If (Count_Rate .le. 0) Then
      Count_Wait = 1
      LogInfoMsg = 'Warning: System Clock does not exist! No call to getStopStatus!'
      Call reportLog(LogInfoMsg, LogStatus)
   End If

   ! define integral tables according to BRDF model used
   Call brdfinit_land(model_land)
   Call brdfinit_ocean(model_ocean)

   ! transform angle to radians
   t_sat_lim = theta_sat_limit*rad
   t_sol_lim = theta_sol_limit*rad
   t_sat_wlim = theta_sat_wlimit*rad
   t_sol_wlim = theta_sol_wlimit*rad
   x_lim = xi_limit*rad

   ! initialise snow flag limit to one valid snowy input observation
   N_snow_limit = 1

   ! calculate the day of the year from the calender date
   day_of_year = julday(month, day_of_month, year)-julday(1, 1, year)+1
   ! calculate the "age" of the previous estimate
   days_last_in = Abs(julday(month, day_of_month, year)-julday(month_in, day_of_month_in, year_in))

   ! full-disk or AERONET processing?
   If (MSGpixY .EQ. 9) Then
      fullDisk = .false.
   Else
      fullDisk = .true.
   End If
   print *, ' Full-disk (T) or AERONET (F): ', fullDisk

   if (LinesRest.gt.0) BlockRest = 1

   Missing_ANG = Missing_SAA(1)
   Missing_COV = Missing_COV_in(1,1)
   Missing_PAR = Missing_PAR_in(1,0)
   Missing_QUA_REF = -1

   ! memory allocation for input variables
   Allocate (                                                               &
        & reflectance(1:MSGpixX, 0:LinesBlock+1, 1:N_channels, 1:N_scenes), &
        & lwcs_mask  (1:MSGpixX, 0:LinesBlock+1, 1:N_channels, 1:N_scenes), &
        & zenith_sol (1:MSGpixX, 0:LinesBlock+1, 1:N_scenes),               &
        & azimuth_sol(1:MSGpixX, 0:LinesBlock+1, 1:N_scenes),               &
        & zenith_sat (1:MSGpixX, 0:LinesBlock+1, 1:N_scenes),               &
        & azimuth_sat(1:MSGpixX, 0:LinesBlock+1, 1:N_scenes),               &
        & wind_speed (1:MSGpixX, 0:LinesBlock+1, 1:N_scenes),               &
        & wind_dir   (1:MSGpixX, 0:LinesBlock+1, 1:N_scenes),               &
        & latitude   (1:MSGpixX, 0:LinesBlock+1),                           &
        & longitude  (1:MSGpixX, 0:LinesBlock+1),                           &
!        & aer_mod    (1:MSGpixX, 0:LinesBlock+1),                           &
        & aer_mod_1  (1:MSGpixX, 0:LinesBlock+1),                           &
        & aer_mod_2  (1:MSGpixX, 0:LinesBlock+1),                           &
        & aer_mod_3  (1:MSGpixX, 0:LinesBlock+1),                           &
        & tot_aod    (1:MSGpixX, 0:LinesBlock+1),                           &
        & bdw_mod_2  (1:MSGpixX, 0:LinesBlock+1),                           &
        & bdw_mod_3  (1:MSGpixX, 0:LinesBlock+1),                           &
        & coast_pix  (1:MSGpixX, 0:LinesBlock+1),                           &
        & STAT=astat)
   If (IsError(astat .eq. 0, 'Error when allocating input variables!', ProcessingStatus_p)) Return

   ! if required allocate input variables for the previous estimate
   If (recursion .and. .not. startseries) Then
      Allocate (                                                                               &
           & brdf_in      (1:MSGpixX, 0:LinesBlock+1, 1:NHistFiles, 1:N_channels, 0:MM      ), &
           & brdf_in_     (1:MSGpixX, 0:LinesBlock+1, 1:NHistFiles, 1:N_channels, 0:MM-1    ), &
           & covariance_in(1:MSGpixX, 0:LinesBlock+1, 1:NHistFiles, 1:N_channels, 0:MM, 0:MM), &
           & quality_in   (1:MSGpixX, 0:LinesBlock+1, 1:NHistFiles, 1:N_channels            ), &
           & age_obs_in   (1:MSGpixX, 0:LinesBlock+1, 1:NHistFiles, 1:N_channels            ), &
           & STAT=astat)
      If (IsError(astat .eq. 0, 'Error when allocating input variables!', ProcessingStatus_p)) Return
   End If

   ! memory allocation for time series variables
   Allocate (mask_cma(1:N_scenes),                    &
      & ang_ok(1:N_scenes), processed(1:N_scenes),    &
      & cloudy(1:N_scenes), snowy(1:N_scenes),        &
      & bad_cma(1:N_scenes), valid(1:N_scenes),       &
      & valid_near(0:boxSize*boxSize*N_scenes-1),     &
      & refl(1:N_scenes), sigrefl(1:N_scenes),        &
      & sigrefl_(1:N_scenes), wi_angular(1:N_scenes), &
      & wi_angular_(1:N_scenes),                      &
      & scat_ang(1:N_scenes), sunglint(1:N_scenes),   &
      & theta_sat(1:N_scenes), theta_sol(1:N_scenes), &
      & t_sat_rel(1:N_scenes), t_sol_rel(1:N_scenes), &
      & phi_sat(1:N_scenes), phi_sol(1:N_scenes),     &
      & phi_del(1:N_scenes),                          &
      & wspeed(1:N_scenes), wdir(1:N_scenes),         &
      & w_ok(1:N_scenes),                             &
      & A(1:N_scenes, 0:MM),                          &
      & AT(0:MM, 1:N_scenes),                         &
      & AT_(0:MM-1, 1:N_scenes),                      &
      & A_(1:N_scenes, 0:MM-1),                       &
      & b(1:N_scenes), b_(1:N_scenes),                &
      & F3_res(0, 0:MM),                              &
      & STAT=astat)
   If (IsError(astat .eq. 0, 'Error when allocating time series variables!', ProcessingStatus_p)) Return

   ! memory allocation for output variables
   If (recursion) Then
      Allocate ( &
           & brdf            (1:MSGpixX, 0:LinesBlock+1, 1:N_channels, 0:MM),       &
           & brdf1           (1:MSGpixX, 0:LinesBlock+1, 1:N_channels, 0:MM-1),     &
           & aodins          (1:MSGpixX, 0:LinesBlock+1, 1:N_channels),             &
           & xiins           (1:MSGpixX, 0:LinesBlock+1),                           &
           & refins          (1:MSGpixX, 0:LinesBlock+1),                           &
           & jacoAODins      (1:MSGpixX, 0:LinesBlock+1),                           &
           & covariance      (1:MSGpixX, 0:LinesBlock+1, 1:N_channels, 0:MM, 0:MM), &
           & age_obs         (1:MSGpixX, 0:LinesBlock+1, 1:N_channels),             &
           & age_obs_al      (1:MSGpixX, 0:LinesBlock+1, 1:N_channels+1),           &
           & CMins           (1:MSGpixX, 0:LinesBlock+1, 1:N_channels),             &
           & quality         (1:MSGpixX, 0:LinesBlock+1, 1:N_channels),             &
           & quality_al      (1:MSGpixX, 0:LinesBlock+1, 1:N_channels+1),           &
           & albedo_sdh      (1:MSGpixX, 0:LinesBlock+1, 1:N_channels),             &
           & sigma_albedo_sdh(1:MSGpixX, 0:LinesBlock+1, 1:N_channels),             &
           & albedo_sbh      (1:MSGpixX, 0:LinesBlock+1, 1:N_channels),             &
           & sigma_albedo_sbh(1:MSGpixX, 0:LinesBlock+1, 1:N_channels),             &
           & albedo_bdh      (1:MSGpixX, 0:LinesBlock+1),                           &
           & sigma_albedo_bdh(1:MSGpixX, 0:LinesBlock+1),                           &
           & albedo_bbh      (1:MSGpixX, 0:LinesBlock+1),                           &
           & sigma_albedo_bbh(1:MSGpixX, 0:LinesBlock+1),                           &
           & albedo_vdh      (1:MSGpixX, 0:LinesBlock+1),                           &
           & sigma_albedo_vdh(1:MSGpixX, 0:LinesBlock+1),                           &
           & albedo_ndh      (1:MSGpixX, 0:LinesBlock+1),                           &
           & sigma_albedo_ndh(1:MSGpixX, 0:LinesBlock+1),                           &
           & STAT=astat)
      If (IsError(astat .eq. 0, 'Error when allocating output variables!', ProcessingStatus_p)) Return
   End If

   ! divide image into successive blocks to be processed
   Do NB = 1, N_Blocks

      ! calculate the record length and the record number
      If (NB .le. N_Blocks_F) Then
         Lines = LinesBlock
      Else
         Lines = LinesRest
      End If

      ! read latitude and longitude
      If (.not. read_tmp_file(Trim(YFILEINP(1)), latkind, MSGpixX, NB, LinesBlock, N_Blocks, &
         & ProcessingStatus_p, latitude (:, 0:LinesBlock+1), "", Missing_LAT)) Return
      If (.not. read_tmp_file(Trim(YFILEINP(2)), latkind, MSGpixX, NB, LinesBlock, N_Blocks, &
         & ProcessingStatus_p, longitude(:, 0:LinesBlock+1), "", Missing_LON)) Return

      ! read aerosol models, total AOD and boundary weights
      If (.not. read_tmp_file(Trim(YFILEINP(3)), aemkind, MSGpixX, NB, LinesBlock, N_Blocks, &
         & ProcessingStatus_p, aer_mod_1(:, 0:LinesBlock+1), "", Missing_AEM)) Return
      If (.not. read_tmp_file(Trim(YFILEINP(4)), aemkind, MSGpixX, NB, LinesBlock, N_Blocks, &
         & ProcessingStatus_p, aer_mod_2(:, 0:LinesBlock+1), "", Missing_AEM)) Return
      If (.not. read_tmp_file(Trim(YFILEINP(5)), aemkind, MSGpixX, NB, LinesBlock, N_Blocks, &
         & ProcessingStatus_p, aer_mod_3(:, 0:LinesBlock+1), "", Missing_AEM)) Return
      If (.not. read_tmp_file(Trim(YFILEINP(6)), taokind, MSGpixX, NB, LinesBlock, N_Blocks, &
         & ProcessingStatus_p, tot_aod  (:, 0:LinesBlock+1), "", Missing_TAO)) Return
      If (.not. read_tmp_file(Trim(YFILEINP(7)), bdwkind, MSGpixX, NB, LinesBlock, N_Blocks, &
         & ProcessingStatus_p, bdw_mod_2  (:, 0:LinesBlock+1), "", Missing_BDW)) Return
      If (.not. read_tmp_file(Trim(YFILEINP(8)), bdwkind, MSGpixX, NB, LinesBlock, N_Blocks, &
         & ProcessingStatus_p, bdw_mod_3  (:, 0:LinesBlock+1), "", Missing_BDW)) Return

      ! read coast pixels
      If (.not. read_tmp_file(Trim(YFILEINP(9)), cstkind, MSGpixX, NB, LinesBlock, N_Blocks, &
         & ProcessingStatus_p, coast_pix(:, 0:LinesBlock+1), "", Missing_CST)) Return

      ! read previous estimate
      If (recursion .and. .not. startseries) Then

         ! loop over previous esimates
         Do IHF = 1, NHistFiles

            ! loop over channels
            Do I = 1, N_channels

               KF = prv_off + (N_channels*(IHF-1) + (I-1))*2

               ! loop over model parameter indices
               Do J = 0, MM
                  ! read parameter
                  If (.not. read_tmp_file(Trim(YFILEINP(KF)), parkind, MSGpixX, NB, LinesBlock, N_Blocks, &
                       & ProcessingStatus_p, brdf_in(:, 0:LinesBlock+1, IHF, I, J), "K"//Char(48+J)//"in", Missing_PAR)) Return
               End Do
               Do J = 0, MM-1
                  ! read parameter
                  If (.not. read_tmp_file(Trim(YFILEINP(KF)), parkind, MSGpixX, NB, LinesBlock, N_Blocks, &
                       & ProcessingStatus_p, brdf_in_(:, 0:LinesBlock+1, IHF, I, J), "R"//Char(48+J)//"in", Missing_PAR)) Return
               End Do

               ! read processing flag
               If (.not. read_tmp_file(Trim(YFILEINP(KF)), quakind, MSGpixX, NB, LinesBlock, N_Blocks, &
                    & ProcessingStatus_p, quality_in(:, 0:LinesBlock+1, IHF, I), "K012.qua_in", Missing_QUA)) Return

               ! read age of information
               If (.not. read_tmp_file(Trim(YFILEINP(KF)), agekind, MSGpixX, NB, LinesBlock, N_Blocks, &
                    & ProcessingStatus_p, age_obs_in(:, 0:LinesBlock+1, IHF, I), "K012.age_in", Missing_AGE)) Return

               ! loop over covariance matrix indices
               J = 0
               Do N = 0, MM
                  Do NN = N, MM
                     J = J+1
                     ! read covariance matrix element
                     If (.not. read_tmp_file(Trim(YFILEINP(KF+1)), covkind, MSGpixX, NB, LinesBlock, N_Blocks, &
                          & ProcessingStatus_p, covariance_in(:, 0:LinesBlock+1, IHF, I, N, NN), "E"//Char(48+J)//"in",     &
                          & Missing_COV)) Return
                  End Do
               End Do
            End Do
         End Do
      End If

      ! loop over scenes
      Do N = 1, N_scenes
         ! loop over channels
         Do I = 1, N_channels !correctif du 2006116
            ! read reflectance
            If (.not. read_tmp_file(Trim(YFILEIN((N-1)*N_files_in_scene+I)), refkind, MSGpixX, NB, LinesBlock, N_Blocks, &
               & ProcessingStatus_p, reflectance(:, 0:LinesBlock+1, I, N), "c"//Char(48+I), Missing_REF(1))) Return

            ! read quality file including land/water/cloud/snow information
            If (.not. read_tmp_file(Trim(YFILEIN((N-1)*N_files_in_scene+I)), maskind, MSGpixX, NB, LinesBlock, N_Blocks, &
               & ProcessingStatus_p, lwcs_mask(:, 0:LinesBlock+1, I, N), "c"//Char(48+I)//".qua", Missing_QUA_REF)) Return
         End Do

         ! read angles
         If (.not. read_tmp_file(Trim(YFILEIN((N-1)*N_files_in_scene+N_channels+1)), angkind, MSGpixX, NB, LinesBlock, N_Blocks, &
            & ProcessingStatus_p, azimuth_sol(:, 0:LinesBlock+1, N), "", Missing_ANG)) Return
         If (.not. read_tmp_file(Trim(YFILEIN((N-1)*N_files_in_scene+N_channels+2)), angkind, MSGpixX, NB, LinesBlock, N_Blocks, &
            & ProcessingStatus_p, zenith_sol (:, 0:LinesBlock+1, N), "", Missing_ANG)) Return
         If (.not. read_tmp_file(Trim(YFILEIN((N-1)*N_files_in_scene+N_channels+3)), angkind, MSGpixX, NB, LinesBlock, N_Blocks, &
            & ProcessingStatus_p, azimuth_sat(:, 0:LinesBlock+1, N), "", Missing_ANG)) Return
         If (.not. read_tmp_file(Trim(YFILEIN((N-1)*N_files_in_scene+N_channels+4)), angkind, MSGpixX, NB, LinesBlock, N_Blocks, &
            & ProcessingStatus_p, zenith_sat (:, 0:LinesBlock+1, N), "", Missing_ANG)) Return

         ! read ECMWF climatos : wind speed (m/s) and direction (radians from north)
         If (.not. read_tmp_file(Trim(YFILEIN((N-1)*N_files_in_scene+N_channels+5)), spdkind, MSGpixX, NB, LinesBlock, N_Blocks, &
            & ProcessingStatus_p, wind_speed (:, 0:LinesBlock+1, N), "", Missing_ANG)) Return
         If (.not. read_tmp_file(Trim(YFILEIN((N-1)*N_files_in_scene+N_channels+6)), angkind, MSGpixX, NB, LinesBlock, N_Blocks, &
            & ProcessingStatus_p, wind_dir   (:, 0:LinesBlock+1, N), "", Missing_ANG)) Return

      End Do

      ! initialize the output variables
      If (recursion) Then
         brdf = MissingValue_kcov
         brdf1 = MissingValue_kcov
         aodins = MissingValue_kcov
         xiins = MissingValue_kcov
         refins = MissingValue_kcov
         jacoAODins = MissingValue_kcov
         covariance = MissingValue_kcov
         quality = MissingValue
         quality_al = MissingValue
         age_obs = MissingValue
         age_obs_al = MissingValue
         CMins = MissingValue
         albedo_sdh = MissingValue
         sigma_albedo_sdh = MissingValue
         albedo_sbh = MissingValue
         sigma_albedo_sbh = MissingValue
         albedo_bdh = MissingValue
         sigma_albedo_bdh = MissingValue
         albedo_bbh = MissingValue
         sigma_albedo_bbh = MissingValue
         albedo_vdh = MissingValue
         sigma_albedo_vdh = MissingValue
         albedo_ndh = MissingValue
         sigma_albedo_ndh = MissingValue
      End If

      If (instantaneous) then
         print *, ' Instantaneous processing of ', Lines, ' lines'
      Else
         print *, ' Daily processing of ', Lines, ' lines'
      End If

      If (writeInversion_flag) Then
         print *, ''
         print * ,'START - Debugging info for improvement of inversion method'
         print * , '*_*_*_*_*_*_*_*_*_*'
         print *, 'X, Y, I, uv, us, xi, albed_b_in, ref_s, ref_tol, ref_tol_out, tau_0'
      EndIf

      !print*, 'Lines: ', Lines
      !print*, 'MSGpixX: ', MSGpixX
      Do Y = 1, Lines

         YY = (NB-1)*LinesBlock + Y

         If (YY*100 .ge. ipct*MSGpixY .and. Lines .gt. 1) Then
            If (.not. writeInversion_flag) Write (*, FMT='(A,"  Processed : ",I3,"%",TL4)', Advance='no') ACHAR(13), ipct
            ipct = ipct+1
         End If

         Call System_Clock(Count_New, Count_Rate, Count_Max)
         If ((Count_New-Count_Old) .gt. Count_Wait .or. (Count_New-Count_Old) .lt. 0) Then
            Count_Old = Count_New
            Call getStopStatus(StopStatus)
            If (StopStatus .eq. STOP_EXECUTION) Then
               LogInfoMsg = 'Algorithm stopped by ICARE system!'
               Call reportLog(LogInfoMsg, LogStatus)
               Call stopping
               ProcessingStatus_p = UNABLE_TO_PROCESS
               Return
            Endif
         End If

         Do X = 1, MSGpixX

            ! Are coast pixels processed?
            if (.not. COAST_OK .and. coast_pix(X, Y) .eq. 1) Then
               If (instantaneous) Then
                  quality(X, Y, :) = 5
                  CMins(X, Y, :) = 10 ! special case
               EndIf
               ! no processing
               cycle
            EndIf

            ! in AERONET model we process only the central pixel (#5), the rest being used for spatiotemporal smoothing only
            If ((instantaneous) .and. (.not. fullDisk) .and. (Y .ne. 5)) cycle

            ! propagate land/water mask information from the last input scene
            lw_mask = iand(lwcs_mask(X, Y, 1, 1), B'00000011')

            ocean_flag = .false.
            If ((lw_mask .eq. MLW_OCEAN) .or. (lw_mask .eq. MLW_WATER)) ocean_flag = .true.
            ! Are ocean pixels processed?
            If (.not. OCEAN_OK .and. ocean_flag) cycle

            land_flag = .false.
            If (lw_mask .eq. MLW_LAND) land_flag = .true.
            ! Are land pixels processed?
            If (.not. LAND_OK .and. land_flag) cycle

            ! Reading a priori AOD from CAMS-based climatology
            tau_force = (tot_aod(X, Y)+Offset_TAO)/Scale_TAO
            ! tau_force .eq. 0.0 for pixels out of the disk
            if (tau_force .eq. 0.0) cycle
            ! checking validity of values
            if (tau_force .lt. 0.0) then
               print*, 'Negative CAMS AOD?', X, Y, YY, tot_aod(X, Y), Offset_TAO, Scale_TAO, tau_force
               STOP
            endif

            ! Reading a priori aerosol model from CAMS-based climatology
            ! Model 0 is Sayer's 'maritime - Atlantic'
            ! Model 1 is Lyapustin's 'continental - GSFC'
            ! Model 2 is Lyapustin's 'arid areas - western USA'
            ! Model 4 is Lyapustin's 'continental - Europe'
            ! Model 6 is Lyapustin's 'desert - Solar Village'
            ! Model 7 is Lyapustin's 'biomass-burning - Mongu'
            ! Model 8 is Lyapustin's 'industrial - India'

            ! Primary aerosol model
            aer_type_1 = aer_mod_1(X, Y) ! "Aerosol model 1" in MSG+0000.3km.aer.h5
            ! Models 3 and 5 should not be in SEVIRI disk but...
            if (aer_type_1 .eq. 3) aer_type_1  = 1
            if (aer_type_1 .eq. 5) aer_type_1  = 8
            if (aer_type_1 .eq. -1) cycle ! pixels out of the disk
            ! checking validity of values
            if ((aer_type_1 .lt. 0) .or. (aer_type_1 .gt. 8)) then
               print*, 'Unknown aerosol type 1!', X, YY, aer_type_1
               STOP
            endif

            ! Secondary aerosol model
            aer_type_2 = aer_mod_2(X, Y) ! "Aerosol model 2" in MSG+0000.3km.aer.h5
            flag_aer_type_2 = .false.
!            if (aer_type_2 .ne. 255) Then
            if (aer_type_2 .ne. -1) Then
               flag_aer_type_2 = .true.
               ! Models 3 and 5 should not be in SEVIRI disk but...
               if (aer_type_2 .eq. 3) aer_type_2  = 1
               if (aer_type_2 .eq. 5) aer_type_2  = 8
               if ((aer_type_2 .lt. 0) .or. (aer_type_2 .gt. 8)) then
                  print*, 'Unknown aerosol type 2!', X, YY, aer_type_2
                  STOP
               endif
               weight_mod_2 = (bdw_mod_2(X, Y)+Offset_BDW)/Scale_BDW ! "Boundary weight 2" in MSG+0000.3km.aer.h5
               ! Needed values are Offset_BDW=0, Scale_BDW=10000
               if ((weight_mod_2 .lt. 0.0) .or. (weight_mod_2 .gt. 1.0)) then
                  print*, 'Weight for model 2 out of [0-1] limits!', X, YY, aer_type_2, weight_mod_2
                  STOP
               endif
            EndIf

            ! Tertiary aerosol model
            aer_type_3 = aer_mod_3(X, Y) ! "Aerosol model 3" in MSG+0000.3km.aer.h5
            flag_aer_type_3 = .false.
!            if (aer_type_3 .ne. 255) Then
            if (aer_type_3 .ne. -1) Then
               flag_aer_type_3 = .true.
               ! Models 3 and 5 should not be in SEVIRI disk but...
               if (aer_type_3 .eq. 3) aer_type_3  = 1
               if (aer_type_3 .eq. 5) aer_type_3  = 8
               if ((aer_type_3 .lt. 0) .or. (aer_type_3 .gt. 8)) then
                  print*, 'Unknown aerosol type 3!', X, YY, aer_type_3
                  STOP
               endif
               weight_mod_3 = (bdw_mod_3(X, Y)+Offset_BDW)/Scale_BDW ! "Boundary weight 3" in MSG+0000.3km.aer.h5
               ! Needed values are  Offset_BDW=0, Scale_BDW=10000
               if ((weight_mod_3 .lt. 0.0) .or. (weight_mod_3 .gt. 1.0)) then
                  print*, 'Weight for model 3 out of [0-1] limits!', X, YY, aer_type_3, weight_mod_3
                  STOP
               endif
            EndIf

            If (flag_aer_type_3 .and. .not. flag_aer_type_2) Then
               print*, 'Absence of secondary model but tertiary model exists!'
               STOP
            endif

            ! Three possibilities with aerosol models
            If (flag_aer_type_3) Then ! Three models exist
               weight_mod_1 = 1.0 - weight_mod_2 - weight_mod_3
               Call aerosol_parameters_init(aer_type_1,weight_mod_1,aer_type_2,weight_mod_2,aer_type_3,weight_mod_3)
               aer_type_prev_1 = -1
            Else If (flag_aer_type_2) Then ! Two models exist
               weight_mod_1 = 1.0 - weight_mod_2
!               Call aerosol_parameters_init(aer_type_1,weight_mod_1,aer_type_2,weight_mod_2,255,0.0)
               Call aerosol_parameters_init(aer_type_1,weight_mod_1,aer_type_2,weight_mod_2,-1,0.0)
               aer_type_prev_1 = -1
            Else ! Mono-mode case
!               If (aer_type_1 .ne. aer_type_prev_1) Call aerosol_parameters_init(aer_type_1,1.0,255,0.0,255,0.0)
               If (aer_type_1 .ne. aer_type_prev_1) Call aerosol_parameters_init(aer_type_1,1.0,-1,0.0,-1,0.0)
               aer_type_prev_1 = aer_type_1
            EndIf

            AOD_max = AOD_max_steps(N_AOD_max) ! 3.0 in the current version
            AOD_min = AOD_max_steps(1) ! 0.0 in the current version

            If (recursion) Then
               ! different use of quality() for DAILY and INST modes
               If (instantaneous) Then
                  quality(X, Y, :) = 0
               Else
                  quality(X, Y, :) = lw_mask
               EndIf
               quality_al(X, Y, :) = lw_mask
            End If

            If (ocean_flag) Then
               model = model_ocean
               hint = hint_ocean
               bihi_calc = bihi_calc_ocean
               albed_force = albed_force_ocean
               sig_k_reg = sig_k_reg_ocean
               k_reg = k_reg_ocean
            Else If (land_flag) Then
               model = model_land
               hint = hint_land
               bihi_calc = bihi_calc_land
               albed_force = albed_force_land
               sig_k_reg = sig_k_reg_land
               k_reg = k_reg_land
            Else
               ! Space pixels (not to be processed?)
               cycle
            EndIf

            If ((latitude(X, Y) .ne. Missing_LAT) .and. (longitude(X, Y) .ne. Missing_LON)) Then

               ! determine solar angle at local solar noon
               lat = (latitude(X, Y)+Offset_LAT)/Scale_LAT
               lon = (longitude(X, Y)+Offset_LON)/Scale_LON
               Call solzenith(day_of_year, lat, theta_sol_midi)
               theta_ref_dh = Minval((/theta_sol_midi, theta_ref_dh_limit/))
               ! linear interpolation in look-up-table
               T = Int(theta_ref_dh/theta_step)
               alpha_intrp = theta_ref_dh/theta_step-T
               If (T .ge. n_table-1) Then
                  T = n_table-2
                  alpha_intrp = 1.
               End If
               dihi_calc(0) = 1.
               dihi_calc(1) = (1.-alpha_intrp)*hint(T, 1)+alpha_intrp*hint(T+1, 1)
               dihi_calc(2) = (1.-alpha_intrp)*hint(T, 2)+alpha_intrp*hint(T+1, 2)

               ! determine the valid observations in the time series
               ! extract time series and apply the scale factors
               ang_ok = zenith_sat(X, Y, :) .ne. Missing_VZA .and. &
                  &   zenith_sol(X, Y, :) .ne. Missing_SZA .and.   &
                  &   azimuth_sat(X, Y, :) .ne. Missing_VAA .and.  &
                  &   azimuth_sol(X, Y, :) .ne. Missing_SAA

               Where (ang_ok)
                  theta_sat = (zenith_sat(X, Y, :)+Offset_VZA)/Scale_VZA*rad
                  t_sat_rel = theta_sat/t_sat_wlim
                  Where (t_sat_rel .gt. 0.999) t_sat_rel = 0.999
                  theta_sol = (zenith_sol(X, Y, :)+Offset_SZA)/Scale_SZA*rad
                  t_sol_rel = theta_sol/t_sol_wlim
                  Where (t_sol_rel .gt. 0.999) t_sol_rel = 0.999
! int16 input azimuth angles are in [-180,180] because fortran can't deal with unsigned integers; need to convert in [0,360]
!                  phi_sat = (azimuth_sat(X, Y, :)+Offset_VAA)/Scale_VAA*rad
!                  Where (phi_sat .gt. two_pi) phi_sat = two_pi
!                  phi_sol = (azimuth_sol(X, Y, :)+Offset_SAA)/Scale_SAA*rad
!                  Where (phi_sol .gt. two_pi) phi_sol = two_pi
                  phi_sat = modulo((azimuth_sat(X,Y,:) + Offset_VAA)/Scale_VAA * rad - two_pi, two_pi)
                  phi_sol = modulo((azimuth_sol(X,Y,:) + Offset_SAA)/Scale_SAA * rad - two_pi, two_pi)
                  phi_del = phi_sat-phi_sol
                  Where (phi_del .lt. 0.) phi_del = phi_del+two_pi
                  Where (phi_del .gt. pi) phi_del = two_pi-phi_del
                  ! calculation of scattering angle - 0deg: fwd scatt. 180deg: bck scatt.
                  scat_ang = pi-Acos(Cos(theta_sat)*Cos(theta_sol)+Sin(theta_sat)*Sin(theta_sol)*Cos(phi_del))
               End Where

               ! calculate angular weights via scattering angle
               if (land_flag) Then
                  Where (ang_ok)
!                     wi_angular = scat_ang/pi ! for retrieving AOD (Inv_1)
!                     wi_angular_ = (pi-scat_ang)/pi ! for retrieving BRDF (Inv_2)
                     wi_angular = (pi-xi_trunc)/(10.0*(pi-scat_ang)) ! for retrieving AOD (Inv_1)
                     wi_angular_ = (pi-xi_trunc)/(10.0*(scat_ang-xi_trunc)) ! for retrieving BRDF (Inv_2)                    
                  End Where
               Else If (ocean_flag) Then! no weights over ocean to avoid band-shaped artifact
                  Where (ang_ok)
                     wi_angular = 1.0 ! for retrieving AOD
                     wi_angular_ = 1.0 ! not used for ocean surfaces as Inv_1 == Inv_2
                  End Where
               Else
                  print*, 'Pixel is not land, nor ocean... what is it then? Stopping...'
                  STOP
               EndIf

               ang_ok = (ang_ok .and. theta_sol .le. t_sol_lim .and. theta_sat .le. t_sat_lim)

               ! observations with scattering angle lower than 'xi_trunc' are not considered due to the truncation of the phase function
               ang_ok = (ang_ok .and. scat_ang .gt. xi_trunc)

               ! removing sunglint-affected pixels for ocean only.
               if (ocean_flag) Then
                  Where (ang_ok)
                     ! calculation of Sun glint angle
                     sunglint = Acos((Cos(theta_sol)*Cos(theta_sat))-(Sin(theta_sol)*Sin(theta_sat)*Cos(phi_del)))
                  End Where
                  ang_ok = (ang_ok .and. sunglint .gt. x_lim)
                  ! second condition to consider the bigger sunglint hole for high SZA (not implemented for DAILY multi-geo study)
                  ! ang_ok = (ang_ok .and. sunglint .gt. x_lim .and. Cos(sunglint) .lt. 0.9*Cos(theta_sol))
               End if

               valids = .false.
               snow = .false.
               albed_d = 0.
               albed_b = 0.
               sigma_d = 0.
               sigma_b = 0.

               ! reading wind data
               wspeed = (wind_speed(X, Y, :)+Offset_SWND)/Scale_SWND
               wdir   = (wind_dir  (X, Y, :)+Offset_DWND)/Scale_DWND

               ! checking that wind data exist
               w_ok   = ((wind_speed(X, Y, :) .ne. Missing_SWND) .and. (wind_dir(X, Y, :) .ne. Missing_DWND))

               ! loop over channels
               Do I = 1, N_channels

                  debug_flag = .false.
!                  if ((I .eq. 1) .and. (Y .eq. 5)) debug_flag = .true.
                  if ((I .eq. 1) .and. (Y .eq. 5)) debug_flag = .false.

                  ! check for correctly processed slots
                  processed = BTest(lwcs_mask(X, Y, I, :), BIT_PROC) .and. reflectance(X, Y, I, :) .ne. Missing_Ref
                  ! check cloud mask information
                  mask_cma = iand(lwcs_mask(X, Y, I, :), B'00011100')
                  cloudy = mask_cma .eq. MCL_CLOUD .or. mask_cma .eq. MCL_CONTAM
                  snowy = mask_cma .eq. MCL_SNOW .or. mask_cma .eq. MCL_SNOW_X
                  bad_cma = mask_cma .eq. MCL_CLEAR_X .or. mask_cma .eq. MCL_SNOW_X
                  ! define valid scenes
                  valid = ang_ok .and. processed .and. .not. cloudy .and. .not. snowy
                  ! valid scenes need wind data over ocean
                  if (ocean_flag) valid = valid .and. w_ok
                  ! eliminate observations with bad CMa quality
                  If (instantaneous) Then
                    If (bad_CMa_elim_inst) valid = valid .and. .not. bad_cma
                  Else
                    If (bad_CMa_elim_daily) valid = valid .and. .not. bad_cma
                  EndIf
                  ! eliminate potentially affected observations before or after cloudy slots for DAILY
                  !If (.not. instantaneous) Then
                  If (N_slot_elim .gt. 0) Then
                     Do N = 1, N_scenes
                        If (cloudy(N)) valid(Maxval((/N-N_slot_elim, 1/)):Minval((/N+N_slot_elim, N_scenes/))) = .false.
                     End Do
                  End If
                  !End If
                  ! calculate the number of valid observations
                  N_valid_obs = Count(valid)
                  ! DAILY runs only if there are N_obs_limit observations at least for pixel X,Y during one full day
                  if (.not. instantaneous) observations_d = (N_valid_obs .ge. N_obs_limit)
                  ! INST runs only if pixel X,Y is available at t0...
                  if (instantaneous) Then
                     if (smooth_flag) Then
                        ! checking validity of slots over last 30-min
                        observations_i = valid(1) .and. (N_valid_obs .ge. NINT(N_scenes/3.0))

                        if (.not. observations_i) Then
                           if (.not. ang_ok(1) .and. .not. ang_ok(2) .and. .not. ang_ok(3)) quality(X, Y, I) = 10
                           if (.not. processed(1) .and. .not. processed(2) .and. .not. processed(3)) quality(X, Y, I) = 12
                           if (ocean_flag .and. (.not. w_ok(1) .and. .not. w_ok(2) .and. .not. w_ok(3))) quality(X, Y, I) = 14
                           if (snowy(1) .and. snowy(2) .and. snowy(3)) quality(X, Y, I) = 20
                           if (cloudy(1) .and. cloudy(2) .and. cloudy(3)) quality(X, Y, I) = 30
                           if (bad_CMa_elim_inst .and. (bad_cma(1) .and. bad_cma(2) .and. bad_cma(3))) quality(X, Y, I) = 32
                        Else
                           ! ... and if there is at least 1/3 of observations in the considered time frame (today: 3 slots out of 9)
                           !observations_i = (N_valid_obs .ge. NINT(N_scenes/3.0))
                           !if (.not. observations_i) quality(X, Y, I) = 40
                           if (N_valid_obs .lt. NINT(N_scenes/3.0)) quality(X, Y, I) = 110

                           if (.not. bad_CMa_elim_inst .and. (bad_cma(1) .or. bad_cma(2) .or. bad_cma(3))) quality(X, Y, I) = 112
                        EndIf
                     Else
                        ! checking validity of current slot only
                        observations_i = valid(1)
                        if (.not. ang_ok(1)) quality(X, Y, I) = 10
                        if (.not. processed(1)) quality(X, Y, I) = 12
                        if (ocean_flag .and. .not. w_ok(1)) quality(X, Y, I) = 14
                        if (snowy(1)) quality(X, Y, I) = 20
                        if (cloudy(1)) quality(X, Y, I) = 30
                        if (bad_CMa_elim_inst .and. bad_cma(1)) quality(X, Y, I) = 32
                     EndIf
                  End If

                  ! calculate number of valid snowy observations
                  N_valid_snow = Count(snowy .and. valid)
                  If (.not. snow_flag_one) N_snow_limit = N_valid_obs/2+1

                  ! regularisation
                  ! land :
                  !    k_reg_land      =  0.0, 0.03,    0.2,  0.0
                  !    sig_k_reg_land  = 10.0, 0.05,    0.5, 50.0
                  ! ocean (k0 is the water reflectance, k1 is the fraction of Fresnel reflection, k2 does not exist, and k3 is the AOD) :
                  !    k_reg_ocean     = 0.01,  1.0,    0.0,  0.0
                  !    sig_k_reg_ocean =  0.1,  1.0, 0.0001, 50.0
                  CkI_reg = 0.
                  Do N = 0, MM
                     CkI_reg(N, N) = sig_k_reg(N)**(-2)
                  End Do
                  k_reg_(0)=k_reg(0)
                  k_reg_(1)=k_reg(1)
                  k_reg_(2)=k_reg(2)

                  ! if applicable generate result based on previous estimate as a priori information
                  If (recursion) Then

                     If (.not. startseries) Then
                        previous = BTest(quality_in(X, Y, 1, I), BIT_MSG) .and. .not. BTest(quality_in(X, Y, 1, I), BIT_FAILS)
                     Else
                        previous = .false.
                     End If

                     ! if previous estimate is available read parameter vector and covariance matrix
                     If (previous) Then

                        ! tau_force becomes tau_in if stable during the last NHistFiles (3) days
                        tau_in_hist = 0.
                        tau_in_hist_flag = .true.
                        Do NS = 1, NHistFiles
                           tau_in_hist(NS)=brdf_in(X, Y, NS, I, 3) / Scale_PAR_in(I, 3)
                           if ((tau_in_hist(NS) .le. AOD_min) .or. (tau_in_hist(NS) .ge. AOD_max)) tau_in_hist_flag = .false.
                        EndDo
                        if (tau_in_hist_flag) then
                           tau_in_hist_mean=SUM(tau_in_hist)/NHistFiles
                           tau_in_hist_std = 0.
                           DO NS = 1, NHistFiles
                              tau_in_hist_std = tau_in_hist_std + ABS(tau_in_hist(NS)-tau_in_hist_mean)**2
                           End Do
                           tau_in_hist_std = SQRT(tau_in_hist_std/NHistFiles)

                           if (tau_in_hist_std .le. tau_in_hist_flag_val) tau_force = tau_in_hist_mean
                        end if

                        ! use previous result as a priori information
                        k_in = brdf_in(X, Y, 1, I, :)/Scale_PAR_in(I, :)

                        flag_kinIsClimato = ((k_in(0) .eq. k_reg(0)) .and. (k_in(1) .eq. k_reg(1)) .and. (k_in(2) .eq. k_reg(2)))

                        k_in_(0:2) = brdf_in_(X, Y, 1, I, :)/Scale_PAR_in(I,0:MM-1)

                        !tau_in = k_in(3) not used anymore but definition is left for clarity

                        Ck_in = 0
                        J = 0
                        Do N = 0, MM
                           Do NN = N, MM
                              J = J+1
                              Ck_in(N, NN) = (covariance_in(X, Y, 1, I, N, NN)/Scale_Cov_in(I, J))**2
                              If (N .ne. NN) Then
                                 Ck_in(N, NN) = 0.
                                 Ck_in(NN, N) = Ck_in(N, NN)
                              End If
                           End Do
                        End Do

                        !If (land_flag) then
                        Do NN = 0, 2
                           If (Ck_in(NN, NN) .le. (0.0015)**2 .and. Ck_in(NN, NN) .gt. (0.00001)**2) then
                              if (debug_flag) print *, 'dispersive factor over land&ocean to constrain the error of priori information'
                              Do N = 0, 2
                                 Ck_in(N, N) = Ck_in(N, N)*(0.0015)**2/Ck_in(NN, NN)
                              End Do
                           EndIF
                        EndDo
                        !EndIf
                        ! assure positive variances
                        Do N = 0, MM
                           Ck_in(N, N) = Maxval((/Ck_in(N, N), epsilon_var/))
                        End Do

                        ! add "temporal uncertainty"
                        Do N = 0, MM
                           Do NN = 0, MM
                              if (N .eq. 3 .or. NN .eq. 3) then
                                 Ck_in(N, NN) = Ck_in(N, NN)*(50.)**days_last_in
                              else
                                 if (N .ne. 0 .and. NN .ne. 0) then
                                    ! add "temporal uncertainty"
                                    delta = 2.**(2./(timescale*12.))-1. ! timescale*12. = 60
                                    Ck_in(N, NN) = Ck_in(N, NN)*(1.+delta)**days_last_in
                                 else
                                    if (ocean_flag) delta = 2.**(2./(timescale*12.))-1. ! timescale*12 = 60
                                    if (land_flag) delta = 2.**(2./(timescale*2.))-1. ! timescale*12 = 10
                                    Ck_in(N, NN) = Ck_in(N, NN)*(1.+delta)**days_last_in
                                 endif
                              endif
                           End Do
                        End Do
                        ! assure positive variances
                        Do N = 0, MM
                           Ck_in(N, N) = Maxval((/Ck_in(N, N), epsilon_var/))
                        End Do

                        ! calculate the inverse
                        CkI_in = invers(Ck_in)
                        if (CkI_in(1, 1) .eq. -111.) then
                           tau = -1.
                           CkI_in = 0
                           Do N = 0, 2
                              CkI_in(N, N) = 50000
                           End Do
                          !!goto 301
                        endif
                     Else
                        ! no previous estimate available
                        k_in = 0.
                        CkI_in = 0.
                        flag_kinIsClimato = .false. ! to avoid horizontal-shaped artefact for pixels without a priori
                     End If !!end if previous

                     !-*******************************************************************
                     ! if observations are available calculate matrix A and vector b
                     If (observations_i .or. observations_d) Then
                        ! calculate reflectance value and error estimate
                        Where (valid) refl = (reflectance(X, Y, I, :)+Offset_Ref)/Scale_Ref
                        Where (valid) refl = refl*cos(theta_sol)*(sqrt(cos(theta_sol)**2+2.0*eta_sph+eta_sph**2)-cos(theta_sol))/eta_sph
                        Where (valid) refl = refl*cos(theta_sat)*(sqrt(cos(theta_sat)**2+2.0*eta_sph+eta_sph**2)-cos(theta_sat))/eta_sph
!                       Where ( valid ) sigrefl = (sig_nadir_a(I)+sig_nadir_b(I)*refl)
!                       Where ( valid .and. sigrefl .lt. sigrefl_min) sigrefl = sigrefl_min
!                       Where ( valid .and. sigrefl .gt. sigrefl_max) sigrefl = sigrefl_max
!                       Where ( valid ) sigrefl = sigrefl !* wi_angular -> this was found to have a negative impact on BRDF retrieval
                        Where (valid) sigrefl = wi_angular
                        Where (valid .and. bad_cma) sigrefl = sigrefl*bad_CMa_factor
                        Where (valid) sigrefl_ = wi_angular_
                        Where (valid .and. bad_cma) sigrefl_ = sigrefl_*bad_CMa_factor
                     Endif !!end if observations

                     !************************************************************************
                     ! perform calculation if observations and/or previous estimate are available
                     If ((observations_i .or. observations_d) .or. previous) Then

                        if (instantaneous) then
!!$ BEG inst_estimates
                           ! INST inversion is made if there are enough observations
                           If (observations_i) Then
                              ! INST inversion is made if a priori is available
                              If (previous) Then
                                 ! A priori must be different from climatology over land (not over ocean since reflectance is forced)
                                 If (.not. (flag_kinIsClimato .and. land_flag)) Then

                                    ! if AP is older than 25 days...
                                    If (age_obs_in(X, Y, 1, I) .gt. 25) quality(X, Y, I) = 114

                                    If (smooth_flag) Then ! spatiotemporal smoothing
                                       X_X(:) = -1
                                       Y_Y(:) = -1
                                       if (fullDisk) Then ! full disk
                                          ! nearest pixels in boxSize x boxSize box
                                          cptX = 0
                                          Do NS = 1, boxSize
                                             Do S = 1, boxSize
                                                If ( ((X+NS-halfBoxSize) .ge. 1) .and. ((Y+S-halfBoxSize) .ge. 1) .and. &
                                                   & ((X+NS-halfBoxSize) .le. MSGpixX) .and. ((Y+S-halfBoxSize) .le. LinesBlock)  ) Then
                                                   ! keeping pixels which are 0.25 degrees away at most (0.25deg -> 28km at the equator)
                                                   If ( (((latitude(X+NS-halfBoxSize,Y+S-halfBoxSize)+Offset_LAT)/Scale_LAT)+0.25 .ge. lat) .and. &
                                                      & (((latitude(X+NS-halfBoxSize,Y+S-halfBoxSize)+Offset_LAT)/Scale_LAT)-0.25 .le. lat) .and. &
                                                      & (((longitude(X+NS-halfBoxSize,Y+S-halfBoxSize)+Offset_LON)/Scale_LON)+0.25 .ge. lon) .and. &
                                                      & (((longitude(X+NS-halfBoxSize,Y+S-halfBoxSize)+Offset_LON)/Scale_LON)-0.25 .le. lon) ) Then
                                                      If (NS .ne. halfBoxSize .or. S .ne. halfBoxSize) Then
                                                         X_X(cptX) = X+NS-halfBoxSize
                                                         Y_Y(cptX) = Y+S-halfBoxSize
                                                         cptX = cptX+1
                                                      EndIf
                                                   EndIf
                                                EndIf
                                             EndDo
                                          EndDo
                                       else ! AERONET subset
                                          Do NS = 0, boxSize*boxSize-1
                                             X_X(NS) = X ! station position
                                             Y_Y(NS) = NS+1 ! pixel position
                                          EndDo
                                          cptX = boxSize*boxSize ! this may not be always good. Now it is since AERONET subsets correspond to 3x3 boxes
                                       endif
                                       cptX_final = cptX * N_scenes ! number of pixels in spatiotemporal superpixel

                                       !print*, 'X_X: ', X_X
                                       !print*, 'Y_Y: ', Y_Y
                                       !print*, 'cptX: ', cptX
                                       !print*, 'cptX_final: ', cptX_final

                                       ! initialisation for spatiotemporal average
                                       uv_near(:) = -1.
                                       us_near(:) = -1.
                                       phFunc_near(:,:) = -1.
                                       ref_tol_near(:) = -1.
                                       ref_s_near(:) = -1.
                                       ref_s_near_(:) = -1.
                                       albed_b_in_near(:) = -1.
                                       albed_b_in_near_(:) = -1.
                                       weight_near(:) = -1.
                                       xi_near(:) = -1.
                                       valid_near(:) = .false.

                                       ! Check valid obs                                       
                                       i_smooth = -1
                                       Do N = 1, N_scenes
                                          Do NS = 0, cptX-1 ! spatiotemporal loop
                                             
                                             i_smooth = (N-1)*cptX+NS ! index

                                             ! check for correctly processed slots
                                             !print*, 'X_X(NS), Y_Y(NS), I: ', X_X(NS), Y_Y(NS), I
                                             processed = BTest( lwcs_mask(X_X(NS),Y_Y(NS),I,N), BIT_PROC ) .and. &
                                                & reflectance(X_X(NS),Y_Y(NS),I,N) .ne. Missing_Ref

                                             ! check cloud mask information
                                             mask_cma(N) = iand(lwcs_mask(X_X(NS),Y_Y(NS),I,N),B'00011100')
                                             cloudy  (N) = mask_cma(N) .eq. MCL_CLOUD   .or. mask_cma(N) .eq. MCL_CONTAM
                                             snowy   (N) = mask_cma(N) .eq. MCL_SNOW    .or. mask_cma(N) .eq. MCL_SNOW_X
                                             bad_cma (N) = mask_cma(N) .eq. MCL_CLEAR_X .or. mask_cma(N) .eq. MCL_SNOW_X
                                             ! determine if scene is valid
                                             valid_near(i_smooth) = ang_ok(N) .and. processed(N) .and. .not. cloudy(N) .and. .not. snowy(N)
                                             ! eliminate observations without brdf_in_ or brdf_in
                                             valid_near(i_smooth) = valid_near(i_smooth) .and. (brdf_in_(X_X(NS),Y_Y(NS),1,I,0) .ne. MissingValue_kcov)
                                             valid_near(i_smooth) = valid_near(i_smooth) .and. (brdf_in (X_X(NS),Y_Y(NS),1,I,0) .ne. MissingValue_kcov)
                                             ! eliminate observations with bad CMa quality
                                             If (bad_CMa_elim_inst) valid_near(i_smooth) = valid_near(i_smooth) .and. .not. bad_cma(N)
                                             ! eliminate observations not belonging to the same surface type LAND/OCEAN/WATER
                                             if (lw_mask .eq. MLW_LAND) Then
                                                valid_near(i_smooth) = valid_near(i_smooth) .and. &
                                                   & (iand(lwcs_mask(X_X(NS), Y_Y(NS), 1, 1), B'00000011') .eq. MLW_LAND)
                                             EndIf
                                             if (lw_mask .eq. MLW_OCEAN) Then
                                                valid_near(i_smooth) = valid_near(i_smooth) .and. &
                                                   & (iand(lwcs_mask(X_X(NS), Y_Y(NS), 1, 1), B'00000011') .eq. MLW_OCEAN)
                                             EndIf
                                             if (lw_mask .eq. MLW_WATER) Then
                                                valid_near(i_smooth) = valid_near(i_smooth) .and. &
                                                   & (iand(lwcs_mask(X_X(NS), Y_Y(NS), 1, 1), B'00000011') .eq. MLW_WATER)
                                             EndIf

                                             If (valid_near(i_smooth)) Then
                                                ! angles
!                                                phi_sat_near   = ( azimuth_sat(X_X(NS),Y_Y(NS),N)+Offset_VAA(N) ) / Scale_VAA(N) * rad
!                                                If (phi_sat_near .gt. two_pi) phi_sat_near = two_pi
!                                                phi_sol_near   = ( azimuth_sol(X_X(NS),Y_Y(NS),N)+Offset_SAA(N) ) / Scale_SAA(N) * rad
!                                                If (phi_sol_near .gt. two_pi) phi_sol_near = two_pi
                                                phi_sat_near   = modulo((azimuth_sat(X_X(NS),Y_Y(NS),N)+Offset_VAA(N) ) / Scale_VAA(N) * rad - two_pi, two_pi)
                                                phi_sol_near   = modulo((azimuth_sol(X_X(NS),Y_Y(NS),N)+Offset_SAA(N) ) / Scale_SAA(N) * rad - two_pi, two_pi)
                                                phi_del_near   = phi_sat_near - phi_sol_near
                                                If (phi_del_near .lt. 0.) phi_del_near = phi_del_near + two_pi
                                                If (phi_del_near .gt. pi) phi_del_near = two_pi - phi_del_near
                                                theta_sat_near = ( zenith_sat(X_X(NS),Y_Y(NS),N)+Offset_VZA(N) ) / Scale_VZA(N) * rad
                                                theta_sol_near = ( zenith_sol(X_X(NS),Y_Y(NS),N)+Offset_SZA(N) ) / Scale_SZA(N) * rad
                                                uv_near(i_smooth) = Cos(theta_sat_near)
                                                us_near(i_smooth) = Cos(theta_sol_near)
                                                xi_near(i_smooth) = pi-Acos(Cos(theta_sat_near)*Cos(theta_sol_near)+Sin(theta_sat_near)*Sin(theta_sol_near)*Cos(phi_del_near))

                                                ! checking that xi is lower than 'xi_trunc'
                                                if (xi_near(i_smooth) .le. xi_trunc) Then
                                                   valid_near(i_smooth) = .false.
                                                   cptX_final = cptX_final - 1
                                                   go to 134
                                                EndIf

                                                ! removing sunglint-affected ocean pixels
                                                if (ocean_flag) Then
                                                   ! calculation of Sun glint angle
                                                   sunglint_near = Acos((Cos(theta_sol_near)*Cos(theta_sat_near))-(Sin(theta_sol_near)*Sin(theta_sat_near)*Cos(phi_del_near)))
                                                   if (sunglint_near .le. x_lim) Then
                                                   !.and. (Cos(sunglint_near) .ge. 0.9*Cos(theta_sol_near))) ! second condition
                                                      valid_near(i_smooth) = .false.
                                                      cptX_final = cptX_final - 1
                                                      go to 134
                                                   EndIf
                                                End if

                                                ! removing coast pixels
                                                if (.not. COAST_OK .and. coast_pix(X_X(NS),Y_Y(NS)) .eq. 1) then
                                                   valid_near(i_smooth) = .false.
                                                   cptX_final = cptX_final - 1
                                                   go to 134
                                                endif

                                                ! phase function
                                                phFunc_near(i_smooth,:) = phaseFunc_value(xi_near(i_smooth),I)

                                                ! ref TOL
                                                ref_tol_near(i_smooth) = (reflectance(X_X(NS),Y_Y(NS),I,N)+Offset_Ref(N)) / Scale_Ref(N)! where does 'Scale_Ref' come from?

                                                ref_tol_near(i_smooth) = &
                                 & ref_tol_near(i_smooth)*us_near(i_smooth)*(sqrt(us_near(i_smooth)**2+2.0*eta_sph+eta_sph**2)-us_near(i_smooth))/eta_sph
                                                ref_tol_near(i_smooth) = &
                                 & ref_tol_near(i_smooth)*uv_near(i_smooth)*(sqrt(uv_near(i_smooth)**2+2.0*eta_sph+eta_sph**2)-uv_near(i_smooth))/eta_sph

                                                ! surface
                                                k_in_near (i_smooth,0:2) = brdf_in (X_X(NS),Y_Y(NS),1,I,0:MM-1) / Scale_PAR_in(I,0:MM-1)
                                                k_in_near_(i_smooth,0:2) = brdf_in_(X_X(NS),Y_Y(NS),1,I,:     ) / Scale_PAR_in(I,0:MM-1)
                                                k_in_near (i_smooth,3) = 0. ! security just in case
                                                k_in_near_(i_smooth,3) = 0. ! security just in case
                                                ! k0 and albed_b_in are to 'albed_force_ocean' for ocean. Coast are an exception (for SUBSET only, as not processed in FULLDISK)
                                                If (ocean_flag .and. coast_pix(X, Y) .ne. 1) k_in_near(i_smooth,0) = albed_force_ocean
                                                If (ocean_flag .and. coast_pix(X, Y) .ne. 1) k_in_near_(i_smooth,0) = albed_force_ocean
                                                albed_b_in_near(i_smooth) = Dot_Product(k_in_near(i_smooth,:), bihi_calc)
                                                If (ocean_flag .and. coast_pix(X, Y) .ne. 1) albed_b_in_near(i_smooth) = albed_force_ocean
                                                albed_b_in_near(i_smooth) = Minval( (/albed_b_in_near(i_smooth), alb_max/) )
                                                albed_b_in_near(i_smooth) = Maxval( (/albed_b_in_near(i_smooth), alb_min/) )
                                                albed_b_in_near_(i_smooth) = Dot_Product(k_in_near_(i_smooth,:), bihi_calc)
                                                If (ocean_flag .and. coast_pix(X, Y) .ne. 1) albed_b_in_near_(i_smooth) = albed_force_ocean
                                                albed_b_in_near_(i_smooth) = Minval( (/albed_b_in_near_(i_smooth), alb_max/) )
                                                albed_b_in_near_(i_smooth) = Maxval( (/albed_b_in_near_(i_smooth), alb_min/) )
                                                ref_s_near(i_smooth) = &
                                                   & dot_product(brdfmodel(theta_sat_near, theta_sol_near, phi_sat_near, phi_sol_near, phi_del_near, &
                                                   & wspeed(1), wdir(1), I, model, ocean_flag),k_in_near(i_smooth,0:2))
                                                ref_s_near_(i_smooth) = &
                                                   & dot_product(brdfmodel(theta_sat_near, theta_sol_near, phi_sat_near, phi_sol_near, phi_del_near, &
                                                   & wspeed(1), wdir(1), I, model, ocean_flag),k_in_near_(i_smooth,0:2))

                                                ! age of obs
                                                time_near(i_smooth) = N_scenes - N + 1 ! reverted for definition of weights

!                                               print*, 'Debug ref_s_near:'
!                                               print*, ref_s_near(i_smooth), albed_b_in_near(i_smooth)
!                                               print*, brdf_in_(X_X(NS),Y_Y(NS),1,I,:)
!                                               print*, brdf_in (X_X(NS),Y_Y(NS),1,I,:)
!                                               print*, X_X(NS),Y_Y(NS),I
!                                               print*, Missing_Ref
                                             Else ! not valid
                                                cptX_final = cptX_final - 1
                                             EndIf
134                                          continue
                                          EndDo ! end loop boxSize x boxSize
                                       EndDo ! end loop N_scenes

                                       ! Calculating average reflectance
                                       ref_tol_near_mean = 0.
                                       i_smooth = -1
                                       Do N = 1, N_scenes
                                          Do NS = 0, cptX-1 ! spatiotemporal loop
                                             i_smooth = (N-1)*cptX+NS ! index
                                             If (valid_near(i_smooth)) ref_tol_near_mean =  ref_tol_near_mean + ref_tol_near(i_smooth)
                                          EndDo
                                       EndDo
                                       ref_tol_near_mean = ref_tol_near_mean / cptX_final

                                       ! Calculating standard deviation
                                       ref_tol_near_std = 0.
                                       i_smooth = -1
                                         Do N = 1, N_scenes
                                          Do NS = 0, cptX-1 ! spatiotemporal loop
                                             i_smooth = (N-1)*cptX+NS ! index
                                             If (valid_near(i_smooth)) ref_tol_near_std = ref_tol_near_std + ABS(ref_tol_near(i_smooth) - ref_tol_near_mean)**2
                                          EndDo
                                       EndDo
                                       ref_tol_near_std = SQRT(ref_tol_near_std / cptX_final)

                                       ! Procedure to filter out bright pixels (ref > mean + std)
                                       i_smooth = -1
                                       Do N = 1, N_scenes
                                          Do NS = 0, cptX-1 ! spatiotemporal loop
                                             i_smooth = (N-1)*cptX+NS ! index
                                             If (valid_near(i_smooth)) Then
                                                If (ref_tol_near(i_smooth) .gt. ref_tol_near_mean+ref_tol_near_std) Then
                                                   valid_near(i_smooth) = .false.
                                                   ref_tol_near(i_smooth) = -1 ! for the weighting later
                                                EndIf
                                             End If
                                          EndDo
                                       EndDo

                                       ! variables needed for INST retrieval
                                       uv = 0.
                                       us = 0.
                                       ref_s = 0.
                                       ref_s_ = 0.
                                       ref_tol = 0.
                                       phFunc(:) = 0.
                                       albed_b_in = 0.
                                       albed_b_in_ = 0.
                                       xi = 0.
                                       weight_near_sum = 0.

                                       ! Spatiotemporal average using ref_tol as weight
                                       i_smooth = -1
                                       Do N = 1, N_scenes
                                          Do NS = 0, cptX-1 ! spatiotemporal loop
                                             i_smooth = (N-1)*cptX+NS ! index

                                             If (valid_near(i_smooth)) Then
                                                ! spatial weights are greater for low reflectances (2 for MINVALUE and 0 for MAXVALUE)
                                                weight_near(i_smooth) = &
                                                   & 2.0*(MAXVAL(ref_tol_near,MASK=valid_near)-ref_tol_near(i_smooth))/ &
                                                   & (MAXVAL(ref_tol_near,MASK=valid_near)-MINVAL(ref_tol_near,MASK=valid_near))
                                                ! temporal weights are greater if close in time (x2 for current slot, x1.7 for 1h ago and x1 for 2h ago)
                                                weight_near(i_smooth) = &
                                                   & weight_near(i_smooth) * (1.0+((time_near(i_smooth)-1)/(N_scenes-1))**0.5)
                                               ! print*, X, Y, I, N, N_scenes, weight_near(i_smooth), time_near(i_smooth)

                                                ! weighted mean
                                                uv = uv + uv_near(i_smooth)*weight_near(i_smooth)
                                                us = us + us_near(i_smooth)*weight_near(i_smooth)
                                                ref_s = ref_s + ref_s_near(i_smooth)*weight_near(i_smooth)
                                                ref_s_ = ref_s_ + ref_s_near_(i_smooth)*weight_near(i_smooth)
                                                ref_tol = ref_tol + ref_tol_near(i_smooth)*weight_near(i_smooth)
                                                phFunc = phFunc + phFunc_near(i_smooth,:)*weight_near(i_smooth)
                                                albed_b_in = albed_b_in + albed_b_in_near(i_smooth)*weight_near(i_smooth)
                                                albed_b_in_ = albed_b_in_ + albed_b_in_near_(i_smooth)*weight_near(i_smooth)
                                                xi = xi + xi_near(i_smooth)*weight_near(i_smooth)

                                                ! sum of weights
                                                weight_near_sum = weight_near_sum + weight_near(i_smooth)
                                             End If

                                          End Do
                                       End Do
                                       uv = uv / weight_near_sum
                                       us = us / weight_near_sum
                                       ref_s = ref_s / weight_near_sum
                                       ref_s_ = ref_s_ / weight_near_sum
                                       ref_tol = ref_tol / weight_near_sum
                                       phFunc = phFunc / weight_near_sum
                                       albed_b_in = albed_b_in / weight_near_sum
                                       albed_b_in_ = albed_b_in_ / weight_near_sum
                                       xi = xi / weight_near_sum

                                    Else ! .not. smooth_flag

                                       uv = Cos(theta_sat(1))
                                       us = Cos(theta_sol(1))
                                       If (ocean_flag .and. coast_pix(X, Y) .ne. 1) k_in(0) = albed_force_ocean
                                       If (ocean_flag .and. coast_pix(X, Y) .ne. 1) k_in_(0) = albed_force_ocean
                                       ref_s_ = dot_product(brdfmodel(theta_sat(1), theta_sol(1), phi_sat(1), phi_sol(1), phi_del(1), &
                                          & wspeed(1), wdir(1), I, model, ocean_flag), k_in_(0:2))
                                       ref_s = dot_product(brdfmodel(theta_sat(1), theta_sol(1), phi_sat(1), phi_sol(1), phi_del(1), &
                                          & wspeed(1), wdir(1), I, model, ocean_flag), k_in(0:2))
                                       ref_tol = refl(1)
                                       phFunc = phaseFunc_value(scat_ang(1),I)
                                       albed_b_in_ = Dot_Product(k_in_, bihi_calc)
                                       If (ocean_flag .and. coast_pix(X, Y) .ne. 1) albed_b_in_ = albed_force_ocean
                                       albed_b_in_ = Minval((/albed_b_in_, alb_max/))
                                       albed_b_in_ = Maxval((/albed_b_in_, alb_min/))
                                       albed_b_in = Dot_Product(k_in, bihi_calc)
                                       If (ocean_flag .and. coast_pix(X, Y) .ne. 1) albed_b_in = albed_force_ocean
                                       albed_b_in = Minval((/albed_b_in, alb_max/))
                                       albed_b_in = Maxval((/albed_b_in, alb_min/))
                                       xi = scat_ang(1)

                                    EndIf ! If (smooth_flag)

                                    if (debug_flag) Then
                                       !print *, 'uv, us: ', uv, us
                                       !print *, 'ref_s, albed_b_in: ', ref_s, albed_b_in
                                       !print *, 'phFunc, ref_tol: ', phFunc, ref_tol
                                       !print *, 'ref_tol: ', ref_tol
                                       If (smooth_flag) Then
                                          print *, ''
                                          print *, '***Spatiotemporal smoothing ENABLED***'
                                          print *, ''
                                          !print *, 'ref_tol_near: ', ref_tol_near
                                          !print *, '   ref_tol_near_mean, ref_tol_near_std: ', ref_tol_near_mean, ref_tol_near_std, ref_tol_near_mean+ref_tol_near_std
                                          !print *, 'valid_near: ', valid_near
                                          !print *, 'weight_near (%): ', 100*weight_near/weight_near_sum
                                          !print *, 'phFunc_near : ', phFunc_near
                                          !print *, 'ref_s_near : ', ref_s_near
                                          !print *, 'albed_b_in_near : ', albed_b_in_near
                                          !print *, 'k_in_near : ', k_in_near
                                       else
                                          print *, ''
                                          print *, '***Spatiotemporal smoothing DISABLED***'
                                          print *, ''
                                       EndIf
                                    endif

                                    ! Starting instantaneous retrieval
                                    tau_0 = tau_badValue

                                    ! Bright surface?
                                    bright_surface = .False.
                                    if (albed_b_in .ge. bright_alb) bright_surface = .True.

                                    ! using double BRDF model for ref_s
                                    w_ref_s_ = pi/(pi-xi)-1.0
                                    w_ref_s = pi/xi-1.0
                                    ref_s = (w_ref_s*ref_s + w_ref_s_*ref_s_)/(w_ref_s+w_ref_s_)
                                    albed_b_in = (w_ref_s*albed_b_in + w_ref_s_*albed_b_in_)/(w_ref_s+w_ref_s_)

                                    ! security check on surface reflectance
                                    if (ref_s .lt. 0.) Then
                                       ref_s = 0.
                                       quality(X, Y, I) = 124
                                    EndIf
                                    if (ref_s .gt. 1.) Then
                                       ref_s = 1.
                                       quality(X, Y, I) = 126
                                    EndIf

                                    if (debug_flag) Then
                                       print *, '***START INST DEBUG***'
                                       print *, 'X, Y, N_valid_obs: ', X, Y, N_valid_obs
                                       print *, 'k_in: ', k_in
                                       print *, 'ref_tol, ref_s: ', ref_tol, ref_s
                                    endif

                                    ! Calling instantaneous retrieval method
                                    call LM(Sa_val**(1.+ref_s), g_tilde(:,I), ssa(:,I), ssa_tilde(:,I), eta(:,I), &
                                       & uv, us, phFunc, ref_tol, ref_tol_out, ref_s, albed_b_in, &
                                       & tau_force, tau_0, jacoAOD, AOD_max_steps, debug_flag)

                                    ! if AOD is slightly negative we invert again with a higher weight for a priori
                                    if (tau_0 .le. 0.0) Then
                                       call LM(Sa_val**(1.+ref_s)/10.0, g_tilde(:,I), ssa(:,I), ssa_tilde(:,I), eta(:,I), &
                                          & uv, us, phFunc, ref_tol, ref_tol_out, ref_s, albed_b_in, &
                                          & tau_force, tau_0, jacoAOD, AOD_max_steps, debug_flag)

                                       ! higher AP weight if it did not work either
                                       if (tau_0 .le. 0.0) Then
                                          call LM(Sa_val**(1.+ref_s)/50.0, g_tilde(:,I), ssa(:,I), ssa_tilde(:,I), eta(:,I), &
                                             & uv, us, phFunc, ref_tol, ref_tol_out, ref_s, albed_b_in, &
                                             & tau_force, tau_0, jacoAOD, AOD_max_steps, debug_flag)

                                          ! try again with lower surface reflectance if it did not work either
                                          if (tau_0 .le. 0.0) Then
                                             call LM(Sa_val**(1.+ref_s)/50.0, g_tilde(:,I), ssa(:,I), ssa_tilde(:,I), eta(:,I), &
                                                & uv, us, phFunc, ref_tol, ref_tol_out, ref_s*0.85, albed_b_in, &
                                                & tau_force, tau_0, jacoAOD, AOD_max_steps, debug_flag)
                                             if (tau_0 .gt. AOD_min .and. tau_0 .lt. AOD_max) quality(X, Y, I) = 120
                                          Else
                                             if (tau_0 .gt. AOD_min .and. tau_0 .lt. AOD_max) quality(X, Y, I) = 118
                                          EndIf
                                       Else
                                          if (tau_0 .gt. AOD_min .and. tau_0 .lt. AOD_max) quality(X, Y, I) = 116
                                       EndIf
                                    endif

                                    ! repeat inversion if AOD over ocean is negative due to negative values of ref_TOL
                                    if (ocean_flag .and. (ref_tol .gt. -0.1) .and. (ref_tol .le. 0.0) .and. (tau_0 .le. 0.0)) Then
                                       call LM(Sa_val**(1.+ref_s)/10.0, g_tilde(:,I), ssa(:,I), ssa_tilde(:,I), eta(:,I), &
                                          & uv, us, phFunc, 0.0, ref_tol_out, ref_s, albed_b_in, &
                                          & tau_force, tau_0, jacoAOD, AOD_max_steps, debug_flag)

                                       if (tau_0 .gt. AOD_min .and. tau_0 .lt. AOD_max) quality(X, Y, I) = 122
                                    endif

                                    ! checking the robustness of the solution and removing extrema (AOD=AOD_min and AOD=AOD_max)
                                    if (tau_0 .ge. AOD_max .or. tau_0 .le. AOD_min) Then
                                       tau_0 = tau_badValue
                                       quality(X, Y, I) = 75
                                    Else
                                       ! checking that special cases are not met (no_smoothing, bad_cma, etc.)
                                       if (quality(X, Y, I) .eq. 0) quality(X, Y, I) = 100
                                    End if

                                    If (writeInversion_flag) print *, X, Y, I, uv, us, xi, albed_b_in, ref_s, ref_tol, ref_tol_out, tau_0

                                    if (debug_flag) Then
                                       print *, 'Final AOD solution is: ', tau_0
                                       print *, '***END INST DEBUG***'
                                       print *, ' '
                                    endif

                                    ! Building confidence parameter based on Jacobian and surface reflectance
                                    if (tau_0 .eq. tau_badValue) Then
                                       conf_Meas = 0 ! bad retrieval
                                    else
                                       if (ABS(jacoAOD) .gt. jacoAOD_thrs) then
                                          conf_Meas = 6 ! best retrieval
                                       Else If (ABS(jacoAOD) .gt. jacoAOD_thrs/4.0) then
                                          conf_Meas = 4 ! good retrieval
                                       Else
                                          conf_Meas = 2 ! bad retrieval
                                       EndIf
                                       ! bright reflectance makes quality decrease
                                       if (bright_surface) conf_Meas = conf_Meas - 1
                                    EndIf

                                    ! result of inversion
                                    k = 0.0
                                    k(0:2) = k_in_(0:2)
                                    !k(0:2) = k_in(0:2) -> could be a mix but this is never used!
                                    k(3) = tau_0
                                    ! check range and store brdf parameters
                                    Do S = 0, MM
                                       k(S) = Minval((/k(S), par_max/))
                                       k(S) = Maxval((/k(S), par_min/))
                                    End Do
                                    aodins(X, Y, I) = NInt(k(MM)*scale_par, Kind=parkind)
                                    xiins(X, Y) = NInt(xi*scale_par, Kind=parkind)
                                    If (I .eq. 1) then
                                       refins(X, Y) = NInt(ref_s*scale_par, Kind=parkind)
                                       jacoAODins(X, Y) = NInt(jacoAOD*scale_par, Kind=parkind)
                                    EndIf
                                    CMins(X, Y, I) = NInt(conf_Meas, Kind=agekind)
                                    ! store age of information
                                    age_obs(X, Y, I) = 0 ! but there is no notion of age for the instantaneous AOD values!

                                 Else
                                    quality(X, Y, I) = 50
                                    If (writeInversion_flag) print *, X, Y, I, 'BAD'
                                 EndIf ! not climato
                              Else
                                 quality(X, Y, I) = 40
                                 If (writeInversion_flag) print *, X, Y, I, 'BAD'
                              EndIf ! previous
                           Else
                              If (writeInversion_flag) print *, X, Y, I, 'BAD'
                           End If ! observations_i
!!$ END inst_estimates
                        else

                           quality(X, Y, I) = quality(X, Y, I)+QUA_MSG

                           ! calculate bi-hemispherical albedo
                           If (previous) Then
                              albed_b_in = Dot_Product(k_in, bihi_calc)
                              albed_b_in_ = Dot_Product(k_in_, bihi_calc)
                           else
                              albed_b_in = albed_force ! default value
                              albed_b_in_ = albed_force ! default value
                           endIf
                           albed_b_in = Minval((/albed_b_in, alb_max/))
                           albed_b_in = Maxval((/albed_b_in, alb_min/))
                           albed_b_in_ = Minval((/albed_b_in_, alb_max/))
                           albed_b_in_ = Maxval((/albed_b_in_, alb_min/))

                           if (debug_flag) print *, '***START DAILY DEBUG***'
                           if (debug_flag) print *, 'X, Y, I, N_valid_obs: ', X, Y, I, N_valid_obs
                           If (observations_d) Then
!!$ BEG daily_estimates
                              if (debug_flag) print *, 'Enough observations'

                              !no processing of AOD over snow
                              If (N_valid_snow .lt. N_snow_limit) Then
                                 if (debug_flag) print *, 'No snow'
                                 Do N = 0, 2
                                    CkI_in(N, N) = Maxval((/CkI_in(N, N), epsilon_var/))
                                 End Do

                                 ! calculating average TOL reflectance
                                 counter = 0
                                 refl_mean = 0.
                                 Do N = 1, N_scenes
                                    If (valid(N)) then
                                       refl_mean = refl_mean + refl(N)
                                       counter = counter + 1
                                    EndIf
                                 EndDo
                                 refl_mean = refl_mean / counter

                                 ! flags for inversion loop
                                 convergence_OK = .false.
                                 limit_AOD = .false.

                                 CkI_reg_i = 0.0
                                 if (startseries .or. flag_kinIsClimato) then
                                    if (debug_flag) print*, '--> Using regularization to start process'
                                    CkI_reg_i(0, 0) = CkI_reg(0, 0)
                                    CkI_reg_i(1, 1) = CkI_reg(1, 1)
                                    CkI_reg_i(2, 2) = CkI_reg(2, 2)
                                 Else
                                    if (debug_flag) print*, '--> Regularization is not used'
                                    CkI_reg_i(0, 0) = CkI_reg(0, 0)/1000.
                                    CkI_reg_i(1, 1) = CkI_reg(1, 1)/1000.
                                    CkI_reg_i(2, 2) = CkI_reg(2, 2)/1000.
                                 EndIf

                                 If (k_in(3) .eq. tau_badValue) then
                                    if (debug_flag) print*, 'tau_in is not used as there is no a priori AOD'
                                    CkI_in(3, 3) = CkI_in(3, 3)/1000.
                                 EndIf

501                              k(3) = 0.
                                 tau_0 = 0.
                                 tau_0_tilde = -1.
                                 tau_1 = -1.

                                 k(0) = -1.0
                                 k(1) = -1.0
                                 k(2) = -1.0

                                 Ck = 0.
                                 Ck_ite = 0.

                                 cpt_tau = 1
                                 cpt_tau_conv = 8
                                 cpt_tau_max = 30

                                 if (debug_flag) print*, 'Previous k:', k_in

                                 ! begin loop for solution convergence
                                 if (debug_flag) print *, 'Convergence loop starts'
                                 Do While (cpt_tau .le. cpt_tau_max)

                                    convergence_OK = (ABS(k(3)-tau_0) .lt. delta_tau) .and. (ABS(k(3)-tau_1) .lt. delta_tau)
                                    limit_AOD = (k(3) .eq. AOD_min) .or. (k(3) .eq. AOD_max)
                                    persistent_low_AOD = (k(3) .eq. AOD_min) .and. (tau_0 .eq. AOD_min) .and. (tau_1 .eq. AOD_min)
                                    persistent_high_AOD = (k(3) .eq. AOD_max) .and. (tau_0 .eq. AOD_max) .and. (tau_1 .eq. AOD_max)

                                    ! covariance matrix to propate is Ck
                                    if (cpt_tau .eq. cpt_tau_conv) Ck = Ck_ite

                                    ! Exit conditions
                                    if (cpt_tau .gt. cpt_tau_conv .and. persistent_low_AOD) then
                                       if (debug_flag) print*, 'Retrieval did not work: persistent AOD=AOD_min'
                                       k(3) = tau_badValue
                                       exit
                                    endif
                                    if (cpt_tau .gt. cpt_tau_conv .and. persistent_high_AOD) then
                                       if (debug_flag) print*, 'Retrieval did not work: persistent AOD=AOD_max'
                                       k(3) = tau_badValue
                                       exit
                                    endif
                                    if (cpt_tau .gt. cpt_tau_conv .and. convergence_OK .and. .not. limit_aod) then
                                       if (debug_flag) print*, 'Potential convergence!'
                                       exit
                                    endif
                                    if (cpt_tau .eq. cpt_tau_max) then
                                       if (.not. convergence_OK .or. limit_aod) then
                                           if (debug_flag) print*, 'Retrieval did not work: no convergence'
                                           k(3) = tau_badValue
                                           exit
                                       else
                                          if (debug_flag) print*, 'Potential convergence! AT THE LAST ITERATION???'
                                          exit
                                       endif
                                    endif

                                    ! Saving previous inversion
                                    tau_1 = tau_0
                                    tau_0 = k(3)

                                    ! calculation of AOD tilde
                                    aod_pos = MINLOC(ABS(AOD_max_steps-tau_0), DIM=1)
                                    tau_0_tilde = (1.-ssa(aod_pos,I)*eta(aod_pos,I))*tau_0

                                    ! fill matrix A and vector b
                                    NN = 0
                                    Do N = 1, N_scenes
                                       If (valid(N)) Then
                                          NN = NN+1

                                          phFunc = phaseFunc_value(scat_ang(N),I)

                                          ! RTM from Katsev et al. 2010
                                          us = Cos(theta_sol(N))
                                          uv = Cos(theta_sat(N))
                                          x1_tilde = 3.*g_tilde(aod_pos,I)
                                          ! Aerosol transmittances 
                                          TS = exp(-tau_0_tilde*(1.-ssa_tilde(aod_pos,I)*(1.-(1.-g_tilde(aod_pos,I))/2.))/us)
                                          TV = exp(-tau_0_tilde*(1.-ssa_tilde(aod_pos,I)*(1.-(1.-g_tilde(aod_pos,I))/2.))/uv)
                                          ! Spherical albedo
                                          Salb = tau_0_tilde/(tau_0_tilde+4./(3.-x1_tilde))
                                          ! Calculation of aerosol/surface coefficient
                                          trans_coeff = TS*TV/(1.-Salb*albed_b_in)
                                          ! Aerosol reflectance
                                          ! - Single scattering (without SSA and phase function)
                                          rho_1 = 1./(4.*(us+uv))*(1.-exp(-tau_0_tilde*(1./us+1./uv)))
                                          ! - Multiple scattering
                                          R_MS_us = 1.+1.5*us+(1.-1.5*us)*exp(-tau_0_tilde/us)
                                          R_MS_uv = 1.+1.5*uv+(1.-1.5*uv)*exp(-tau_0_tilde/uv)
                                          R_MS = 1.-R_MS_us*R_MS_uv/(4.+(3.-x1_tilde)*tau_0_tilde) &
                                             & +((3.+x1_tilde)*us*uv-2.*(us+uv))*rho_1

                                          A(NN, 0:MM) = brdfmodel_aerosol(theta_sat(N), theta_sol(N), phi_sat(N), phi_sol(N), phi_del(N), &
                                             & wspeed(N), wdir(N), I, model, ocean_flag, tau_0_tilde, g(aod_pos,I), ssa_tilde(aod_pos,I), &
                                             & trans_coeff, eta(aod_pos,I), phFunc(aod_pos))/sigrefl(N)
                                          b(NN) = (refl(N)-R_MS)/trans_coeff/sigrefl(N)
                                          if (A(NN, 3) .eq. 0) then
                                             A(NN, 0:MM) = 0.0000001
                                             b(NN) = 0.0000001
                                          endif

                                       End If
                                    End Do

                                    ! generate matrix ATA
                                    AT = Transpose(A)
                                    ATA = Matmul(AT(0:MM, 1:N_valid_obs), A(1:N_valid_obs, 0:MM))
                                    ! force the symmetry of the matrix ATA
                                    Do N = 0, MM
                                       Do NN = N+1, MM
                                          ATA(N, NN) = 0.5*(ATA(N, NN)+ATA(NN, N))
                                          ATA(NN, N) = ATA(N, NN)
                                       End Do
                                    End Do

                                    CkI_tau = CkI_in
                                    ! Weight of a priori increases with value of tau_0 and albed_b_in
                                    CkI_tau(0, 0) = CkI_in(0, 0)*min(exp(tau_0*(1.+albed_b_in)), 50.)
                                    CkI_tau(1, 1) = CkI_in(1, 1)*min(exp(tau_0*(1.+albed_b_in)), 50.)
                                    CkI_tau(2, 2) = CkI_in(2, 2)*min(exp(tau_0*(1.+albed_b_in)), 50.)
                                    ! a priori when we go beyond a certain number of iterations
                                    if (cpt_tau .gt. cpt_tau_conv) then
                                       CkI_tau(0, 0) = CkI_tau(0, 0)*min(exp(tau_0*(1.+albed_b_in)), 50.)
                                       CkI_tau(1, 1) = CkI_tau(1, 1)*min(exp(tau_0*(1.+albed_b_in)), 50.)
                                       CkI_tau(2, 2) = CkI_tau(2, 2)*min(exp(tau_0*(1.+albed_b_in)), 50.)
                                    endif
                                    if (cpt_tau .gt. cpt_tau_conv*2) then
                                       CkI_tau(0, 0) = CkI_tau(0, 0)*min(exp(tau_0*(1.+albed_b_in)), 50.)
                                       CkI_tau(1, 1) = CkI_tau(1, 1)*min(exp(tau_0*(1.+albed_b_in)), 50.)
                                       CkI_tau(2, 2) = CkI_tau(2, 2)*min(exp(tau_0*(1.+albed_b_in)), 50.)
                                    endif

                                    Do N = 0, 2
                                       CkI_tau(N, N) = Maxval((/CkI_tau(N, N), epsilon_var/))
                                       CkI_reg_i(N, N) = Maxval((/CkI_reg_i(N, N), epsilon_var/))
                                       ATA(N, N) = Maxval((/ATA(N, N), epsilon_var/))
                                    End Do

                                    ! Covariance matrix
                                    Ck_ite = invers(ATA+CkI_tau+CkI_reg_i)
!                                    print*, 'AT(0:MM, 1:N_valid_obs)', AT(0:MM, 1:N_valid_obs)
!                                    print*, 'A(1:N_valid_obs, 0:MM)', A(1:N_valid_obs, 0:MM)
!                                    print*, 'ATA', ATA
!                                    print*, 'CkI_tau', CkI_tau
!                                    print*, 'CkI_reg_i', CkI_reg_i
!                                    print*, 'b(1:N_valid_obs)', b(1:N_valid_obs)

                                    if (Ck_ite(1, 1) .eq. -111.) then
                                       tau_0 = -1.
                                       Ck_ite = 0.
                                       Do N = 0, 2
                                          Ck_ite(N, N) = 50000.
                                       End Do
                                       goto 301
                                    endif

                                    ! use of regulation factor to constrain lambertian reflectance over ocean and to impose positive matrix at the initialisation
                                    k = Matmul(Ck_ite, Matmul(AT(0:MM, 1:N_valid_obs), b(1:N_valid_obs)) + &
                                            & Matmul(CkI_tau, k_in) + Matmul(CkI_reg_i, k_reg))

!                                    print*, 'Ck_ite', Ck_ite
!                                    print*, 'Matmul(AT(0:MM, 1:N_valid_obs), b(1:N_valid_obs))', Matmul(AT(0:MM, 1:N_valid_obs), b(1:N_valid_obs))
!                                    print*, 'Matmul(CkI_tau, k_in)', Matmul(CkI_tau, k_in)
!                                    print*, 'Matmul(CkI_reg_i, k_reg)', Matmul(CkI_reg_i, k_reg)
!                                    print*, '************************'

                                    ! from tau_0_tilde to tau_0...
                                    ! aod_pos = MINLOC(ABS(AOD_max_steps-k(3)), DIM=1) -> Wrong, as 'aod_pos' should correspond to 'tau' and not 'tau_tilde'
                                    ! By doing nothing, we use 'tau_0', which is the previous value of AOD.
                                    ! This is not totally correct either but it is better because AOD should not change much from one iteration to the other
                                    k(3) = k(3)/(1.-ssa(aod_pos,I)*eta(aod_pos,I))
                                  
                                    ! calculating fitting error
                                    aod_pos = MINLOC(ABS(AOD_max_steps-k(3)), DIM=1)
                                    tau_0_tilde = (1.-ssa(aod_pos,I)*eta(aod_pos,I))*k(3)

                                    aod_pos_j = MINLOC(ABS(AOD_max_steps-(k(3)+delta_tau)), DIM=1)
                                    if (aod_pos_j .gt. N_AOD_max) aod_pos_j = N_AOD_max
                                    tau_0_tilde_j = (1.-ssa(aod_pos_j,I)*eta(aod_pos_j,I))*(k(3)+delta_tau)

                                    fit_error = 0.
                                    jacobian = 0.
                                    counter = 0
                                    Do N = 1, N_scenes
                                       If (valid(N)) Then
                                          ref_s = dot_product(brdfmodel(theta_sat(N), theta_sol(N), phi_sat(N), phi_sol(N), phi_del(N), &
                                             & wspeed(N), wdir(N), I, model, ocean_flag), k(0:2))
                                          P_tilde = phaseFunc_value(scat_ang(N),I)/(1.-eta(aod_pos,I))
                                          call calculate_ref_tol(ref_s, P_tilde(aod_pos), ssa_tilde(aod_pos,I), &
                                             & g_tilde(aod_pos,I), tau_0_tilde, albed_b_in, &
                                             & Cos(theta_sat(N)), Cos(theta_sol(N)), ref_tol)
                                          fit_error = fit_error + abs(refl(N)-ref_tol)

                                          P_tilde_j = phaseFunc_value(scat_ang(N),I)/(1.-eta(aod_pos_j,I))
                                          call calculate_ref_tol(ref_s, P_tilde_j(aod_pos_j), ssa_tilde(aod_pos_j,I), &
                                              & g_tilde(aod_pos_j,I), tau_0_tilde_j, albed_b_in, &
                                              & Cos(theta_sat(N)), Cos(theta_sol(N)), ref_tol_j)
 
                                          jacobian = jacobian + ABS(ref_tol_j-ref_tol)/delta_tau

                                          counter = counter+1
                                       End If
                                    End Do

                                    if (debug_flag) print *, '     i: k - Ck_ite(3,3) ; error ; jacobian: ', &
                                            & cpt_tau, ':', k, '-', Ck_ite(3, 3), ';', fit_error/counter, ';', jacobian/counter

                                    ! AOD must remain between AOD_min and AOD_max
                                    if (k(3) .gt. AOD_max) k(3) = AOD_max
                                    if (k(3) .lt. AOD_min) k(3) = AOD_min

                                    ! constraint for processing over ocean: k(0) is between 0 and k0_ocean_max and k(1) is 1
                                    If (ocean_flag) then
                                       if (k(0) .lt. 0.0) k(0) = 0.0
                                       if (k(0) .gt. k0_ocean_max) k(0) = k0_ocean_max
                                       k(1) = 1.0
                                       k(2) = 0.0
                                    EndIF

                                    cpt_tau = cpt_tau+1
                                 EndDo

                                 ! constraint for processing over ocean: k(0) is between 0 and k0_ocean_max and k(1) is 1
                                 If (ocean_flag) then
                                    if (k(0) .lt. 0.0) k(0) = 0.0
                                    if (k(0) .gt. k0_ocean_max) k(0) = k0_ocean_max
                                    k(1) = 1.0
                                    k(2) = 0.0
                                 EndIF

                                 goto 302
301                              print *, '***singular*301**'
302                              continue

                              Endif ! end of no procesingover snow pixels

                              flag_invFailed = .false.

                              ! Inversion fails if AOD is not within the limits
                              flag_invFailed = ((k(3) .le. AOD_min) .or. (k(3) .ge. AOD_max))
                              ! Inversion fails if BRDF is not realistic
                              If (.not. flag_invFailed .and. land_flag) flag_invFailed = ((k(0) .lt. -0.02) .or. (k(1) .lt. 0.) .or. (k(2) .lt. 0.))

                              ! BRDF is estimated only when AOD is low or mild
                              if (.not. flag_invFailed) then
                                 ! BRDF is not estimated for the first time if AOD > 0.5 over land (or AOD > 0.25 over ocean)
                                 if (land_flag .and. (k(3) .ge. 0.5) .and. flag_kinIsClimato) flag_invFailed = .true.
                                 if (ocean_flag .and. (k(3) .ge. 0.25) .and. flag_kinIsClimato) flag_invFailed = .true.

                                 ! BRDF is not updated if AOD > 1.0 over land (or AOD > 0.5 over ocean)
                                 if (land_flag .and. (k(3) .ge. 1.0)) flag_invFailed = .true.
                                 if (ocean_flag .and. (k(3) .ge. 0.5)) flag_invFailed = .true.
                              endif

                              ! Inversion is done with CAMS AOD when:
                              !  (1) retrieval was not possible using climatologic BRDF after 5 days
                              !      This corresponds to very bright regions for which AOD and BRDF cannot be retrieved simultaneously
                              !  (2) retrieval was not possible after 15 days
                              !      This corresponds to medium bright, high AOD, or cloudy regions for which AOD is hard to estimate
                              flag_invClimato = .false.
                              if (flag_invFailed .and. previous) Then
                                 if ((age_obs_in(X, Y, 1, I) .ge. 5-1 .and. flag_kinIsClimato) .or. (age_obs_in(X, Y, 1, I) .ge. 15-1)) then
                                    flag_invClimato = .true.
                                    k(3) = tau_force
                                 end if
                              end if

                              flag_failedInvExit = .false.

                              ! Filter on Ck(3,3) was removed since some good AOD retrievals were being discarded
                              if ( (flag_invFailed .and. .not. flag_invClimato) .or. (N_valid_snow .ge. N_snow_limit) ) then

                                 if (debug_flag .and. ((k(3) .le. AOD_min) .or. (k(3) .ge. AOD_max))) print *, '      Not good as AOD out of limits'
                                 if (debug_flag .and. land_flag .and. (k(0) .lt. -0.02)) print *, '      Not good as k(0) is negative'
                                 if (debug_flag .and. land_flag .and. (k(1) .lt. 0)) print *, '      Not good as k(1) is negative'
                                 if (debug_flag .and. land_flag .and. (k(2) .lt. 0)) print *, '      Not good as k(2) is negative'
                                 if (debug_flag .and. (N_valid_snow .ge. N_snow_limit)) print *, '      Not good as too many slots with snow'

                                 flag_failedInvExit = .true.
                              else

                                 if (ocean_flag) Then
                                    if (debug_flag) print *, '   daily AOD was retrieved. Skipping 2nd inversion for ocean surfaces'
                                    ! updating BRDF with good solution for next iteration
                                    k_(0)=k(0)
                                    k_(1)=k(1)
                                    k_(2)=k(2)
                                    age = 0
                                    goto 311
                                 EndIf

                                 if (debug_flag) Then
                                    if (.not. flag_invClimato) then
                                       print *, '   daily AOD was retrieved. Proceeding with 2nd inversion for land surfaces'
                                    else
                                       print *, '   daily AOD could not be retrieved over land. Using climatology for 2nd inversion'
                                    end if
                                 end if

                                 k_ = -1.0

                                 aod_pos = MINLOC(ABS(AOD_max_steps-k(3)), DIM=1)
                                 tau_0_tilde = (1.-ssa(aod_pos,I)*eta(aod_pos,I))*k(3)

                                 ! fill matrix A_ and vector b_
                                 NN = 0
                                 Do N = 1, N_scenes
                                    If (valid(N)) Then
                                       NN = NN+1

                                       phFunc = phaseFunc_value(scat_ang(N),I)

                                       ! RTM from Katsev et al. 2010
                                       us = Cos(theta_sol(N))
                                       uv = Cos(theta_sat(N))
                                       x1_tilde = 3.*g_tilde(aod_pos,I)
                                       ! Aerosol transmittances
                                       TS = exp(-tau_0_tilde*(1.-ssa_tilde(aod_pos,I)*(1.-(1.-g_tilde(aod_pos,I))/2.))/us)
                                       TV = exp(-tau_0_tilde*(1.-ssa_tilde(aod_pos,I)*(1.-(1.-g_tilde(aod_pos,I))/2.))/uv)
                                       ! Spherical albedo
                                       Salb = tau_0_tilde/(tau_0_tilde+4./(3.-x1_tilde))
                                       ! Calculation of aerosol/surface coefficient
                                       trans_coeff = TS*TV/(1.-Salb*albed_b_in_)
                                       ! Aerosol reflectance
                                       ! - Single scattering (without SSA and phase function)
                                       rho_1 = 1./(4.*(us+uv))*(1.-exp(-tau_0_tilde*(1./us+1./uv)))
                                       ! - Multiple scattering
                                       R_MS_us = 1.+1.5*us+(1.-1.5*us)*exp(-tau_0_tilde/us)
                                       R_MS_uv = 1.+1.5*uv+(1.-1.5*uv)*exp(-tau_0_tilde/uv)
                                       R_MS = 1.-R_MS_us*R_MS_uv/(4.+(3.-x1_tilde)*tau_0_tilde) &
                                          & +((3.+x1_tilde)*us*uv-2.*(us+uv))*rho_1

                                       A_(NN, 0:MM-1) = brdfmodel(theta_sat(N), theta_sol(N), phi_sat(N), phi_sol(N), phi_del(N), &
                                          & wspeed(N), wdir(N), I, model, ocean_flag) / sigrefl_(N)
                                       b_(NN) = (refl(N) - ssa_tilde(aod_pos,I)*phFunc(aod_pos)/(1.-eta(aod_pos,I))*rho_1 &
                                             & - R_MS) / trans_coeff / sigrefl_(N)
                                    End If
                                 End Do

                                 ! generate matrix ATA_
                                 AT_ = Transpose( A_ )
                                 ATA_ = Matmul( AT_(0:MM-1,1:N_valid_obs), A_(1:N_valid_obs,0:MM-1) )
                                 ! force the symmetry of the matrix ATA_
                                 Do N = 0, MM-1
                                    Do NN = N+1, MM-1
                                       ATA_(N,NN) = 0.5*( ATA_(N,NN)+ATA_(NN,N) )
                                       ATA_(NN,N) = ATA_(N,NN)
                                    End Do
                                 End Do

                                 CkI_reg_i_ = 0.0
                                 CkI_reg_i_(0, 0) = CkI_reg_i(0, 0)
                                 CkI_reg_i_(1, 1) = CkI_reg_i(1, 1)
                                 CkI_reg_i_(2, 2) = CkI_reg_i(2, 2)
                                 
                                 If (.not. flag_invClimato) then
                                    !CkI_tau_ = CkI_in_
                                    ! We use the same covariance matrix than for the first retrieval (TO BE IMPROVED?)
                                    CkI_tau_ = CkI_tau(0:2,0:2)
                                 Else
                                    CkI_tau_ = CkI_in(0:2,0:2)
                                 End if

                                 Do N = 0, 2
                                    CkI_tau_(N, N) = Maxval((/CkI_tau_(N, N), epsilon_var/))
                                    CkI_reg_i_(N, N) = Maxval((/CkI_reg_i_(N, N), epsilon_var/))
                                    ATA_(N, N) = Maxval((/ATA_(N, N), epsilon_var/))
                                 End Do

                                 ! Covariance matrix
                                 Ck_ = invers(ATA_+CkI_tau_+CkI_reg_i_)

                                 if (Ck_(1, 1) .eq. -111.) then
                                    print*, 'Singular cov matrix for INST BRDF -> STOP!'
                                    STOP
                                 endif

                                 ! use of regulation factor to constrain lambertian reflectance over ocean and to impose positive matrix at the initialisation
                                 k_ = Matmul(Ck_, Matmul(AT_(:, 1:N_valid_obs), b_(1:N_valid_obs)) + &
                                         & Matmul(CkI_tau_, k_in_(0:2)) + Matmul(CkI_reg_i_, k_reg_))

                                 ! 2nd inversion fails if retrieved BRDF is not realistic
                                 flag_invFailed = .false.
                                 flag_invFailed = ((k_(0) .lt. -0.02) .or. (k_(1) .lt. 0.) .or. (k_(2) .lt. 0.))

                                 ! 1st and 2nd inversions worked
                                 if (.not. flag_invClimato .and. .not. flag_invFailed) age = 0
                                 ! 1st inversion worked but the 2nd inversion failed
                                 if (.not. flag_invClimato .and. flag_invFailed) then
                                    ! updating BRDF with good solution for next iteration
                                    k_(0)=k(0)
                                    k_(1)=k(1)
                                    k_(2)=k(2)
                                    age = 0
                                 end if
                                 ! 1st and 2nd inversions failed
                                 if (flag_invClimato .and. flag_invFailed) flag_failedInvExit = .true.
                                 ! 1st inversion failed but 2nd inversion worked
                                 if (flag_invClimato .and. .not. flag_invFailed) then
                                    ! updating BRDF with good solution for next iteration
                                    k(0)=k_(0)
                                    k(1)=k_(1)
                                    k(2)=k_(2)
                                    !k(3)=tau_badValue -> for the time being we write the CAMS AOD values in the daily aerosol product
                                    age = 0
                                 end if

                                 if (debug_flag) Then

                                    aod_pos = MINLOC(ABS(AOD_max_steps-k(3)), DIM=1)
                                    tau_0_tilde = (1.-ssa(aod_pos,I)*eta(aod_pos,I))*k(3)

                                    print *, 'Final k : ', k
                                    print *, '   Ck(3,3): ', Ck(3, 3)
                                    print *, 'Final k_: ', k_
                                    print *, '   Ck_(3,3): none!'
                                    print *, ''

                                    print *, '(k case): fitting error for each observation is:'
                                    fit_error = 0.
                                    counter = 0
                                    print *, '     N,     refl(N),     ref_tol,     diff.,     ref_s,     xi,     poids'
                                    Do N = 1, N_scenes
                                       If (valid(N)) Then
                                          ref_s = dot_product(brdfmodel(theta_sat(N), theta_sol(N), phi_sat(N), phi_sol(N), phi_del(N), &
                                             & wspeed(N), wdir(N), I, model, ocean_flag), k(0:2))
                                          P_tilde = phaseFunc_value(scat_ang(N),I)
                                          P_tilde = P_tilde/(1.-eta(aod_pos,I)) ! truncation
                                          call calculate_ref_tol(ref_s, P_tilde(aod_pos), ssa_tilde(aod_pos,I), g_tilde(aod_pos,I), &
                                             & tau_0_tilde, albed_b_in, Cos(theta_sat(N)), Cos(theta_sol(N)), ref_tol)
                                          fit_error = fit_error+abs(refl(N)-ref_tol)
                                          counter = counter+1
                                          print *, N, refl(N), ref_tol, refl(N)-ref_tol, ref_s, scat_ang(N), 1.0/sigrefl(N)
                                       End If
                                    End Do
                                    print *, 'Average absolute fit error is: ', fit_error/counter
                                    print *, '******************************'
                                    print *, '(k_ case): fitting error for each observation is:'
                                    fit_error = 0.
                                    counter = 0
                                    print *, '     N,     refl(N),     ref_tol,     diff.,     ref_s_,     xi,     poids'
                                    Do N = 1, N_scenes
                                       If (valid(N)) Then
                                          ref_s_ = dot_product(brdfmodel(theta_sat(N), theta_sol(N), phi_sat(N), phi_sol(N), phi_del(N), &
                                             & wspeed(N), wdir(N), I, model, ocean_flag), k_(0:2))
                                          P_tilde = phaseFunc_value(scat_ang(N),I)
                                          P_tilde = P_tilde/(1.-eta(aod_pos,I)) ! truncation
                                          call calculate_ref_tol(ref_s_, P_tilde(aod_pos), ssa_tilde(aod_pos,I), g_tilde(aod_pos,I), &
                                             & tau_0_tilde, albed_b_in_, Cos(theta_sat(N)), Cos(theta_sol(N)), ref_tol)
                                          fit_error = fit_error+abs(refl(N)-ref_tol)
                                          counter = counter+1
                                          print *, N, refl(N), ref_tol, refl(N)-ref_tol, ref_s_, scat_ang(N), 1.0/sigrefl_(N)
                                       End If
                                    End Do
                                    print *, 'Average absolute fit error is: ', fit_error/counter
                                    print *, '******************************'
                                    print *, '***END DAILY DEBUG***'
                                    print *, ''

                                 End If
                              endif
                           Else ! no new observations available
                              if (debug_flag) print *, 'Not enough observations! Min. number is: ', N_obs_limit
                              flag_failedInvExit = .true.
                           End If

                           if (flag_failedInvExit) Then
                              if (previous) then
                                 if (debug_flag) print *, '   Propagating previous surface solution'
                                 k = k_in
                                 k(3) = tau_badValue
                                 Ck = Ck_in
                                 k_ = k_in_(0:2)
                                 If (age_obs_in(X, Y, 1, I) .le. age_max-days_last_in) Then
                                    age = age_obs_in(X, Y, 1, I)+days_last_in
                                 Else
                                    age = age_max
                                 End If
                                 If (BTest(quality_in(X, Y, 1, I), BIT_SNOW)) Then
                                    snow = .true.
                                    quality(X, Y, I) = quality(X, Y, I)+QUA_SNOW
                                 End If
                              Else
                                 if (debug_flag) print *, '   No previous surface solution. Taking climatological one'
                                 k = k_reg
                                 k(3) = tau_badValue
                                 Ck = 0.
                                 k_ = k_reg_
                                 Do N = 0, MM
                                    Ck(N, N) = sig_k_reg(N)
                                 End Do
                                 age = -1
                              Endif
                           Endif

311                        if (debug_flag) print*, 'End of inversion' ! we come here after succesful ocean AOD inversion

                           if (freezeDaily) Then
                              ! We propagate the last surface retrieval
                              Do N = 0, 2
                                 k(N) = k_in(N)
                                 Ck(N, N) = Ck_in(N, N)
                                 k_(N) = k_in_(N)
                              End Do
                           End If

                           If (.false.) Then
                              quality(X, Y, I) = quality(X, Y, I)+QUA_FAILS
                           Else
                              valids(I) = .true.
                              ! check range and store brdf parameters
                              Do N = 0, MM
                                 k(N) = Minval((/k(N), par_max/))
                                 k(N) = Maxval((/k(N), par_min/))
                              End Do
                              brdf(X, Y, I, :) = NInt(k*scale_par, Kind=parkind)
                              ! check range and store brdf parameters
                              Do N = 0, MM-1
                                 k_(N) = Minval((/k_(N), par_max/))
                                 k_(N) = Maxval((/k_(N), par_min/))
                              End Do
                              brdf1(X, Y, I, :) = NInt(k_*scale_par, Kind=parkind)
                              ! check range and store covariance matrix elements
                              Do N = 0, MM
                                 Ck(N, N) = Maxval((/Ck(N, N), 0./))
                                 sqcvm(N, N) = Sqrt(Ck(N, N))
                                 sqcvm(N, N) = Minval((/sqcvm(N, N), cxx_max/))
                                 sqcvm(N, N) = Maxval((/sqcvm(N, N), cxx_min/))
                              End Do
                              Do N = 0, MM-1
                                 Do NN = N+1, MM
                                    sqcvm(N, NN) = Sign(Sqrt(Abs(Ck(N, NN))), Ck(N, NN))
                                    sqcvm(N, NN) = Minval((/sqcvm(N, NN), cxy_max/))
                                    sqcvm(N, NN) = Maxval((/sqcvm(N, NN), cxy_min/))
                                 End Do
                              End Do
                              covariance(X, Y, I, :, :) = NInt(sqcvm*scale_cov, Kind=covkind)

                              ! store age of information
                              age_obs(X, Y, I) = age

                              ! albedo calculation
                              If (theta_sol_midi .le. theta_sol_midi_limit) Then
                                 ! set albedo quality flag
                                 quality_al(X, Y, I) = quality(X, Y, I)
                                 age_obs_al(X, Y, I) = age_obs(X, Y, I)
                                 ! calculate directional-hemispherical albedo
                                 albed_d(I) = Dot_Product(k, dihi_calc)
                                 albed_d(I) = Minval((/albed_d(I), alb_max/))
                                 albed_d(I) = Maxval((/albed_d(I), alb_min/))
                                 albedo_sdh(X, Y, I) = Nint(albed_d(I)*scale_alb, &
                                    & Kind=albkind)
                                 ! calculate the error of the di.-hemispherical albedo
                                 ! from the covariance matrix of the parameter errors
                                 ! sigma = Sqrt( I^T C I)
                                 variance = Dot_Product(dihi_calc, Matmul(Ck, dihi_calc))
                                 variance = Maxval((/variance, 0./))
                                 sigma_d(I) = Sqrt(variance)
                                 sigma_d(I) = Minval((/sigma_d(I), sig_max/))
                                 sigma_d(I) = Maxval((/sigma_d(I), sig_min/))
                                 sigma_albedo_sdh(X, Y, I) = Nint(sigma_d(I)*scale_sig, &
                                    & Kind=albkind)
                                 ! calculate bi-hemispherical albedo
                                 albed_b(I) = Dot_Product(k, bihi_calc)
                                 albed_b(I) = Minval((/albed_b(I), alb_max/))
                                 albed_b(I) = Maxval((/albed_b(I), alb_min/))
                                 albedo_sbh(X, Y, I) = Nint(albed_b(I)*scale_alb, &
                                    & Kind=albkind)
                                 ! calculate the error of the bi-hemispherical albedo
                                 ! from the covariance matrix of the parameter errors
                                 variance = Dot_Product(bihi_calc, Matmul(Ck, bihi_calc))
                                 variance = Maxval((/variance, 0./))
                                 sigma_b(I) = Sqrt(variance)
                                 sigma_b(I) = Minval((/sigma_b(I), sig_max/))
                                 sigma_b(I) = Maxval((/sigma_b(I), sig_min/))
                                 sigma_albedo_sbh(X, Y, I) = Nint(sigma_b(I)*scale_sig, &
                                    & Kind=albkind)
                              End If
                           End If
!!$ END daily_estimates

                        End If !!end if observations .or. previous

                     End If !!end if instantaneous .or. daily

                  End If !!end if recursion

               End Do

               ! if applicable calculate narrowband to broadband conversion
               If (recursion .and. Count(valids) .eq. N_channels .and. theta_sol_midi .le. theta_sol_midi_limit) Then
                  quality_al(X, Y, N_channels+1) = quality_al(X, Y, N_channels+1)+QUA_MSG
                  If (snow) quality_al(X, Y, N_channels+1) = quality_al(X, Y, N_channels+1)+QUA_SNOW
                  age_obs_al(X, Y, N_channels+1) = Maxval(age_obs_al(X, Y, 1:N_channels))
                  ! broadband directional-hemispherical albedo
                  albed = co_bb(0)+Dot_Product(albed_d(1:N_channels), co_bb(1:N_channels))
                  albed = Minval((/albed, alb_max/))
                  albed = Maxval((/albed, alb_min/))
                  albedo_bdh(X, Y) = Nint(albed*scale_alb, Kind=albkind)
                  ! error of broadband directional-hemispherical albedo
                  sigma = Sqrt((albed*sigma_co)**2+ &
                     & Dot_Product(co_bb(1:N_channels)**2, sigma_d(1:N_channels)**2))
                  sigma = Minval((/sigma, sig_max/))
                  sigma = Maxval((/sigma, sig_min/))
                  sigma_albedo_bdh(X, Y) = Nint(sigma*scale_sig, Kind=albkind)
                  ! broadband bi-hemispherical albedo
                  albed = co_bb(0)+Dot_Product(albed_b(1:N_channels), co_bb(1:N_channels))
                  albed = Minval((/albed, alb_max/))
                  albed = Maxval((/albed, alb_min/))
                  albedo_bbh(X, Y) = Nint(albed*scale_alb, Kind=albkind)
                  ! error of broadband bi-hemispherical albedo
                  sigma = Sqrt((albed*sigma_co)**2+ &
                     & Dot_Product(co_bb(1:N_channels)**2, sigma_b(1:N_channels)**2))
                  sigma = Minval((/sigma, sig_max/))
                  sigma = Maxval((/sigma, sig_min/))
                  sigma_albedo_bbh(X, Y) = Nint(sigma*scale_sig, Kind=albkind)
                  ! visual directional-hemispherical albedo
                  albed = co_vi(0)+Dot_Product(albed_d(1:N_channels), co_vi(1:N_channels))
                  albed = Minval((/albed, alb_max/))
                  albed = Maxval((/albed, alb_min/))
                  albedo_vdh(X, Y) = Nint(albed*scale_alb, Kind=albkind)
                  ! error of visual directional-hemispherical albedo
                  sigma = Sqrt((albed*sigma_co)**2+ &
                     & Dot_Product(co_vi(1:N_channels)**2, sigma_d(1:N_channels)**2))
                  sigma = Minval((/sigma, sig_max/))
                  sigma = Maxval((/sigma, sig_min/))
                  sigma_albedo_vdh(X, Y) = Nint(sigma*scale_sig, Kind=albkind)
                  ! infrared directional-hemispherical albedo
                  albed = co_ni(0)+Dot_Product(albed_d(1:N_channels), co_ni(1:N_channels))
                  albed = Minval((/albed, alb_max/))
                  albed = Maxval((/albed, alb_min/))
                  albedo_ndh(X, Y) = Nint(albed*scale_alb, Kind=albkind)
                  ! error of infrared directional-hemispherical albedo
                  sigma = Sqrt((albed*sigma_co)**2+ &
                     & Dot_Product(co_ni(1:N_channels)**2, sigma_d(1:N_channels)**2))
                  sigma = Minval((/sigma, sig_max/))
                  sigma = Maxval((/sigma, sig_min/))
                  sigma_albedo_ndh(X, Y) = Nint(sigma*scale_sig, Kind=albkind)
               End If

            End If
         End Do ! loop over X
      End Do ! loop over Y

      If (writeInversion_flag) Then
         print * , '*_*_*_*_*_*_*_*_*_*'
         print * ,'END - Debugging info for improvement of inversion method'
         print * , ''
      EndIf

      ! write the result to disk
      If (NB .eq. 1) Then
         status = 'Replace'
      Else
         status = 'Old'
      End If

      If (recursion) Then

        !! begin write instantaneous estimates !!
         if (instoutput .and. .not. startseries) Then

            ! YFILEINS(1) is the current slot that is processed

            ! AOD
            Do I = 1, N_channels
               If (.not. write_tmp_file(Trim(YFILEINS(1)), parkind, MSGpixX, status, NB, LinesBlock, N_Blocks, Lines, &
                  & ProcessingStatus_p, aodins(:, 1:Lines, I), "Ch"//Char(48+I))) Return
            EndDo

            ! scattering angle
            If (.not. write_tmp_file(Trim(YFILEINS(1)), parkind, MSGpixX, status, NB, LinesBlock, N_Blocks, Lines, &
               & ProcessingStatus_p, xiins(:, 1:Lines), "xi")) Return

            ! surface reflectance VIS06 (ref)
            If (.not. write_tmp_file(Trim(YFILEINS(1)), parkind, MSGpixX, status, NB, LinesBlock, N_Blocks, Lines, &
               & ProcessingStatus_p, refins(:, 1:Lines), "ref")) Return

            ! AOD interval with error lower than 10% for VIS06 (jacoAOD)
            If (.not. write_tmp_file(Trim(YFILEINS(1)), parkind, MSGpixX, status, NB, LinesBlock, N_Blocks, Lines, &
               & ProcessingStatus_p, jacoAODins(:, 1:Lines), "jacoAOD")) Return

            ! Confidence measures
            If (.not. write_tmp_file(Trim(YFILEINS(1)), agekind, MSGpixX, status, NB, LinesBlock, N_Blocks, Lines, &
               & ProcessingStatus_p, CMins(:, 1:Lines, 1), "age")) Return

               ! quality flag spectral
            If (.not. write_tmp_file(Trim(YFILEINS(1)), quakind, MSGpixX, status, NB, LinesBlock, N_Blocks, Lines, &
               & ProcessingStatus_p, quality(:, 1:Lines, 1), "qua")) Return
         EndIf
        !! end write instantaneous estimates !!

         Do I = 1, N_channels !correctif du 2006116

            J = 0
            Do N = 0, MM
               ! BRDF model parameters
               If (.not. write_tmp_file(Trim(YFILEOUT((I-1)*N_files_out_channel+1)), parkind, MSGpixX, status, NB, LinesBlock, N_Blocks, Lines, &
                  & ProcessingStatus_p, brdf(:, 1:Lines, I, N), "K"//Char(48+N))) Return

               ! covariance matrix elements
               Do NN = N, MM
                  J = J+1
                  If (.not. write_tmp_file(Trim(YFILEOUT((I-1)*N_files_out_channel+2)), covkind, MSGpixX, status, NB, LinesBlock, N_Blocks, Lines, &
                     & ProcessingStatus_p, covariance(:, 1:Lines, I, N, NN), "E"//Char(48+J))) Return
               End Do
            End Do

            Do N = 0, MM-1
            ! BRDF_ model parameters
               If (.not. write_tmp_file(Trim(YFILEOUT((I-1)*N_files_out_channel+1)), parkind, MSGpixX, status, NB, LinesBlock, N_Blocks, Lines, &
                  & ProcessingStatus_p, brdf1(:, 1:Lines, I, N), "R"//Char(48+N))) Return
            End Do

            ! age of the last observation used
            If (.not. write_tmp_file(Trim(YFILEOUT((I-1)*N_files_out_channel+1)), agekind, MSGpixX, status, NB, LinesBlock, N_Blocks, Lines, &
               & ProcessingStatus_p, age_obs(:, 1:Lines, I), "age")) Return

            ! age of the last observation used for spectral albedo files
            If (.not. write_tmp_file(Trim(YFILEOUT((I-1)*N_files_out_channel+3)), agekind, MSGpixX, status, NB, LinesBlock, N_Blocks, Lines, &
               & ProcessingStatus_p, age_obs_al(:, 1:Lines, I), "age")) Return

            ! spectral directional-hemispherical albedo
            If (.not. write_tmp_file(Trim(YFILEOUT((I-1)*N_files_out_channel+3)), albkind, MSGpixX, status, NB, LinesBlock, N_Blocks, Lines, &
               & ProcessingStatus_p, albedo_sdh(:, 1:Lines, I), "asdh")) Return

            ! error estimate for the spectral di.-hemispherical albedo
            If (.not. write_tmp_file(Trim(YFILEOUT((I-1)*N_files_out_channel+3)), albkind, MSGpixX, status, NB, LinesBlock, N_Blocks, Lines, &
               & ProcessingStatus_p, sigma_albedo_sdh(:, 1:Lines, I), "asdh_err")) Return

            ! spectral bi-hemispherical albedo
            If (.not. write_tmp_file(Trim(YFILEOUT((I-1)*N_files_out_channel+3)), albkind, MSGpixX, status, NB, LinesBlock, N_Blocks, Lines, &
               & ProcessingStatus_p, albedo_sbh(:, 1:Lines, I), "asbh")) Return

            ! error estimate for the spectral bi-hemispherical albedo
            If (.not. write_tmp_file(Trim(YFILEOUT((I-1)*N_files_out_channel+3)), albkind, MSGpixX, status, NB, LinesBlock, N_Blocks, Lines, &
               & ProcessingStatus_p, sigma_albedo_sbh(:, 1:Lines, I), "asbh_err")) Return

            ! quality flag spectral
            If (.not. write_tmp_file(Trim(YFILEOUT((I-1)*N_files_out_channel+1)), quakind, MSGpixX, status, NB, LinesBlock, N_Blocks, Lines, &
               & ProcessingStatus_p, quality(:, 1:Lines, I), "qua")) Return

            ! quality flag albedo spectral
            If (.not. write_tmp_file(Trim(YFILEOUT((I-1)*N_files_out_channel+3)), quakind, MSGpixX, status, NB, LinesBlock, N_Blocks, Lines, &
               & ProcessingStatus_p, quality_al(:, 1:Lines, I), "qua")) Return

         End Do

         ! broadband directional-hemispherical albedo
         If (.not. write_tmp_file(Trim(YFILEOUT(N_channels*N_files_out_channel+1)), albkind, MSGpixX, status, NB, LinesBlock, N_Blocks, Lines, &
            & ProcessingStatus_p, albedo_bdh(:, 1:Lines), "abdh")) Return

         ! error of the broadband directional-hemispherical albedo
         If (.not. write_tmp_file(Trim(YFILEOUT(N_channels*N_files_out_channel+1)), albkind, MSGpixX, status, NB, LinesBlock, N_Blocks, Lines, &
            & ProcessingStatus_p, sigma_albedo_bdh(:, 1:Lines), "abdh_err")) Return

         ! broadband bi-hemispherical albedo
         If (.not. write_tmp_file(Trim(YFILEOUT(N_channels*N_files_out_channel+1)), albkind, MSGpixX, status, NB, LinesBlock, N_Blocks, Lines, &
            & ProcessingStatus_p, albedo_bbh(:, 1:Lines), "abbh")) Return

         ! error of the broadband bi-hemispherical albedo
         If (.not. write_tmp_file(Trim(YFILEOUT(N_channels*N_files_out_channel+1)), albkind, MSGpixX, status, NB, LinesBlock, N_Blocks, Lines, &
              & ProcessingStatus_p, sigma_albedo_bbh(:, 1:Lines), "abbh_err")) Return

         ! visual directional-hemispherical albedo
         If (.not. write_tmp_file(Trim(YFILEOUT(N_channels*N_files_out_channel+1)), albkind, MSGpixX, status, NB, LinesBlock, N_Blocks, Lines, &
            & ProcessingStatus_p, albedo_vdh(:, 1:Lines), "avdh")) Return

         ! error of the visual directional-hemispherical albedo
         If (.not. write_tmp_file(Trim(YFILEOUT(N_channels*N_files_out_channel+1)), albkind, MSGpixX, status, NB, LinesBlock, N_Blocks, Lines, &
            & ProcessingStatus_p, sigma_albedo_vdh(:, 1:Lines), "avdh_err")) Return

         ! infrared directional-hemispherical albedo
         If (.not. write_tmp_file(Trim(YFILEOUT(N_channels*N_files_out_channel+1)), albkind, MSGpixX, status, NB, LinesBlock, N_Blocks, Lines, &
            & ProcessingStatus_p, albedo_ndh(:, 1:Lines), "andh")) Return

         ! error of the infrared directional-hemispherical albedo
         If (.not. write_tmp_file(Trim(YFILEOUT(N_channels*N_files_out_channel+1)), albkind, MSGpixX, status, NB, LinesBlock, N_Blocks, Lines, &
            & ProcessingStatus_p, sigma_albedo_ndh(:, 1:Lines), "andh_err")) Return

         ! quality flag broadband products
         If (.not. write_tmp_file(Trim(YFILEOUT(N_channels*N_files_out_channel+1)), quakind, MSGpixX, status, NB, LinesBlock, N_Blocks, Lines, &
            & ProcessingStatus_p, quality_al(:, 1:Lines, N_channels+1), "qua")) Return

         ! age of information for broadband albedo file
         If (.not. write_tmp_file(Trim(YFILEOUT(N_channels*N_files_out_channel+1)), agekind, MSGpixX, status, NB, LinesBlock, N_Blocks, Lines, &
            & ProcessingStatus_p, age_obs_al(:, 1:Lines, N_channels+1), "age")) Return

      End If
   End Do

   ! memory deallocation
   Deallocate (                                                      &
      & reflectance, lwcs_mask, zenith_sat, zenith_sol,              &
      & azimuth_sat, azimuth_sol, latitude, longitude,               &
      & mask_cma, ang_ok, processed, cloudy, snowy, bad_cma, valid,  &
      & refl, sigrefl, sigrefl_, wi_angular, wi_angular_,            &
      & scat_ang, sunglint, theta_sat,                               &
      & theta_sol, t_sat_rel, t_sol_rel, phi_sat, phi_sol, phi_del,  &
      & wspeed, wdir, A, AT, b, b_, aer_mod_1, aer_mod_2, aer_mod_3, &
      & bdw_mod_2, bdw_mod_3, tot_aod, wind_speed,                   &
      & wind_dir, coast_pix, w_ok, STAT=astat)
   If (IsError(astat .eq. 0, 'Error when deallocating variables!', ProcessingStatus_p)) Return

   If (recursion .and. .not. startseries) Then
      Deallocate (brdf_in, brdf_in_, covariance_in, quality_in, age_obs_in, STAT=astat)
      If (IsError(astat .eq. 0, 'Error when deallocating variables!', ProcessingStatus_p)) Return
   End If

   If (recursion) Then
      Deallocate (                                                                                     &
         & brdf, brdf1, aodins, xiins, refins, jacoAODins, covariance, age_obs, age_obs_al, quality,  &
         & quality_al, albedo_sdh, sigma_albedo_sdh, albedo_sbh, sigma_albedo_sbh,                     &
         & albedo_bdh, sigma_albedo_bdh, albedo_bbh, sigma_albedo_bbh,                                 &
         & albedo_vdh, sigma_albedo_vdh, albedo_ndh, sigma_albedo_ndh,                                 &
         & STAT=astat)
      If (IsError(astat .eq. 0, 'Error when deallocating variables!', ProcessingStatus_p)) Return
   End If

   Call stopping
   ProcessingStatus_p = PROCESS_OK

End Subroutine start_exe
