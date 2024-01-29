Module algoconf
!  Use nmlst
  Use py_ifc
!  Use hdf5_types

!!$  ! number of scenes
!!$  Integer, Parameter :: N_scenes_max = 96
!!$
!!$  ! number of channels
!!$  Integer, Parameter :: N_channels = 3
!!$
!!$  ! number of model parameters MM=2 -> 3 model parameters
!!$  Integer, Parameter :: MM = 3   

  ! lower limit for allowed values of variances and determinants
  Real, Parameter :: epsilon_var  = 1.E-6
  Real, Parameter :: epsilon_det2 = 1.E-12
  Real, Parameter :: epsilon_det3 = 1.E-18

!!$  ! number of files valid for all scenes
!!$  Integer, Parameter :: N_files_in_all = 1 + 2*N_channels
!!$
!!$  ! number of input files per scene
!!$  Integer, Parameter :: N_files_in_scene = 7
!!$
!!$  ! number of output files per channel ( recursion K012, CK, AL )
!!$  Integer, Parameter :: N_files_out_channel = 3
!!$
!!$  ! number of output files per channel ( composition K012, CK )
!!$  Integer, Parameter :: N_files_out_channel_compo = 2
!!$  
!!$  ! number of output files broadband ( ALBEDO )
!!$  Integer, Parameter :: N_files_out_broad = 1
!!$  
!!$  ! number of input file names to be read from the product configuration file
!!$  Integer, Parameter :: N_FILEINP = N_scenes_max * N_files_in_scene + &
!!$       &                            N_files_in_all
!!$  ! number of output file names ( recursion K012, CK, AL + ALBEDO )
!!$  Integer, Parameter :: N_fileout_rec = N_channels * N_files_out_channel + &
!!$       &                               N_files_out_broad
!!$  ! number of output file names to be read from the product configuration file
!!$  Integer, Parameter :: N_FILEOUT = N_fileout_rec + N_channels * N_files_out_channel_compo
!!$
!!$  ! default value for non-assigned filenames
!!$  Character(LEN=255), Parameter :: filename_undefined = '-'

!!$  ! variables for reading the ProductConfigurationFile
!!$  Integer :: MODE              ! algorithm mode, not implemented
!!$  Integer :: YFREQINSECS       ! frequency for calling GetStopStatus
!!$  Character(LEN=255), Dimension(1:N_FILEINP) :: YFILEINP = filename_undefined
!!$  Character(LEN=255), Dimension(1:N_FILEOUT) :: YFILEOUT

  Logical, Parameter :: NO_CLOSE          = .false.
  Logical, Parameter :: DO_CLOSE          = .true.

  ! variables containing the file attributes
  Character(LEN=12), Dimension(1:N_scenes_max*Max_Nfiles_in_scene) :: attr_product
  Character(LEN=12), Dimension(1:N_scenes_max*Max_Nfiles_in_scene) :: attr_im_aq_time
  Integer, Dimension(1:N_scenes_max*Max_Nfiles_in_scene)           :: attr_sp_ch_id = 0

  ! flag for I/O-status
  Integer :: ios

  ! flag for allocate-status
  Integer :: astat

  ! product acronyms in the HDF attribute list
  Character(LEN=3) :: attr_lat = 'LAT'
  Character(LEN=3) :: attr_lon = 'LON'
  Character(LEN=4) :: attr_brf = 'BRF'
  Character(LEN=3) :: attr_vaa = 'VAA'
  Character(LEN=3) :: attr_vza = 'VZA'
  Character(LEN=3) :: attr_saa = 'SAA'
  Character(LEN=3) :: attr_sza = 'SZA'
 
!  Character(LEN=10), Dimension(N_channels) :: attr_par = &
!       &       (/ 'AL-C1-K012', 'AL-C2-K012', 'AL-C3-K012' /)
!  Character(LEN=8), Dimension(N_channels) :: attr_cov = &
!       &      (/ 'AL-C1-CK', 'AL-C2-CK', 'AL-C3-CK' /)
!  Character(LEN=5), Dimension(N_channels) :: attr_asp = &
!       &      (/ 'AL-C1', 'AL-C2', 'AL-C3' /)
  Character(LEN=10), Dimension(MaxNChannels) :: attr_par = &
       &       (/ 'AL-C1-K012', 'AL-C2-K012', 'AL-C3-K012' /)
  Character(LEN=8), Dimension(MaxNChannels) :: attr_cov = &
       &      (/ 'AL-C1-CK', 'AL-C2-CK', 'AL-C3-CK' /)
  Character(LEN=5), Dimension(MaxNChannels) :: attr_asp = &
       &      (/ 'AL-C1', 'AL-C2', 'AL-C3' /)
  Character(LEN=6) :: attr_alb = 'ALBEDO'

  ! dataset name in the HDF file
  Character(LEN=3)  :: data_lat   = 'LAT'
  Character(LEN=3)  :: data_lon   = 'LON'
  Character(LEN=7)  :: data_brf_d = 'BRF-TOC'
  Character(LEN=10) :: data_brf_q = 'BRF_Q_Flag'
  Character(LEN=3)  :: data_vaa   = 'VAA'
  Character(LEN=3)  :: data_vza   = 'VZA'
  Character(LEN=3)  :: data_saa   = 'SAA'
  Character(LEN=3)  :: data_sza   = 'SZA'
  
  Character(LEN=2), Dimension(0:MM) :: data_par = (/ 'K0', 'K1', 'K2', 'K3'/)
  Character(LEN=6)  :: data_qua   = "Q-Flag"
  Character(LEN=5)  :: data_age   = "Z_Age"
  Character(LEN=3), Dimension(1:MM*MM+1) :: data_cov =  &
       &       (/ 'C00', 'C01', 'C02', 'C03', 'C11', 'C12', 'C13', 'C22', 'C23', 'C33' /)

!!$  Real              :: Scale_LAT                       ! scale factor latitude -> degrees
!!$  Real              :: Offset_LAT
!!$  Integer           :: Missing_LAT

!!$  Real, Dimension(1:N_channels,0:MM)    :: Scale_PAR_in   ! scale factor model parameters
!!$  Real, Dimension(1:N_channels,0:MM)    :: Offset_PAR_in
!!$  Integer, Dimension(1:N_channels,0:MM) :: Missing_PAR_in
!!$
!!$  Real, Dimension(1:N_channels,1:MM*(MM+1))    :: Scale_COV_in   ! scale factor covariance matrix
!!$  Real, Dimension(1:N_channels,1:MM*(MM+1))    :: Offset_COV_in
!!$  Integer, Dimension(1:N_channels,1:MM*(MM+1)) :: Missing_COV_in
!!$
!!$  Real, Dimension(:), Allocatable       :: Scale_REF    ! scale factor reflectance -> 1.
!!$  Real, Parameter    :: Scale_conv_REF = 100.          ! conversion for scale factor [% -> 1.]
!!$  Real, Dimension(:), Allocatable       :: Offset_REF
!!$  Integer, Dimension(:), Allocatable    :: Missing_REF
!!$
!!$  Real, Dimension(:), Allocatable       :: Scale_SAA    ! scale factor solar azimuth angle -> degrees
!!$  Real, Dimension(:), Allocatable       :: Offset_SAA
!!$  Integer, Dimension(:), Allocatable    :: Missing_SAA
!!$
!!$  Real, Dimension(:), Allocatable       :: Scale_SZA    ! scale factor solar zenith angle -> degrees
!!$  Real, Dimension(:), Allocatable       :: Offset_SZA
!!$  Integer, Dimension(:), Allocatable    :: Missing_SZA
!!$
!!$  Real, Dimension(:), Allocatable       :: Scale_VAA    ! scale factor view azimuth angle -> degrees
!!$  Real, Dimension(:), Allocatable       :: Offset_VAA
!!$  Integer, Dimension(:), Allocatable    :: Missing_VAA
!!$
!!$  Real, Dimension(:), Allocatable       :: Scale_VZA    ! scale factor view zenith angle -> degrees
!!$  Real, Dimension(:), Allocatable       :: Offset_VZA
!!$  Integer, Dimension(:), Allocatable    :: Missing_VZA

  ! input: land/water mask values
  Integer, Parameter :: MLW_SPACE    = B'00000010'
  Integer, Parameter :: MLW_OCEAN    = B'00000000'
  Integer, Parameter :: MLW_WATER    = B'00000011'
  Integer, Parameter :: MLW_LAND     = B'00000001'

  ! input: "AL1-cloud mask" values
  Integer, Parameter :: MCL_NOMASK   = B'00000000'
  Integer, Parameter :: MCL_CLEAR    = B'00000100'
  Integer, Parameter :: MCL_CONTAM   = B'00001000'
  Integer, Parameter :: MCL_CLOUD    = B'00001100'
  Integer, Parameter :: MCL_SNOW     = B'00010000'
  Integer, Parameter :: MCL_NOTCLASS = B'00010100'
  Integer, Parameter :: MCL_CLEAR_X  = B'00011000'
  Integer, Parameter :: MCL_SNOW_X   = B'00011100'

  ! input: atmospheric correction processing flag
  Integer, Parameter :: BIT_PROC     = 7

  ! output: quality flag values
  Integer, Parameter :: QUA_MSG      = B'00000100'
  Integer, Parameter :: QUA_EPS      = B'00001000'
  Integer, Parameter :: QUA_APRIORI  = B'00010000'
  Integer, Parameter :: QUA_SNOW     = B'00100000'
  Integer, Parameter :: QUA_FAILS    = B'10000000'

  Integer, Parameter :: BIT_MSG      = 2
  Integer, Parameter :: BIT_SNOW     = 5
  Integer, Parameter :: BIT_FAILS    = 7

  ! input data variables
!!$  Integer, Parameter :: refkind = 2
!!$  Integer, Parameter :: maskind = 1
!!$  Integer, Parameter :: angkind = 2
!!$  Integer, Parameter :: latkind = 2
!!$  Integer, Parameter :: inkind  = 2
  Integer(Kind=refkind), Dimension(:,:,:,:), Allocatable :: reflectance
  Integer(Kind=maskind), Dimension(:,:,:,:), Allocatable :: lwcs_mask  
  Integer(Kind=angkind), Dimension(:,:,:),   Allocatable :: zenith_sol
  Integer(Kind=angkind), Dimension(:,:,:),   Allocatable :: azimuth_sol
  Integer(Kind=angkind), Dimension(:,:,:),   Allocatable :: zenith_sat
  Integer(Kind=angkind), Dimension(:,:,:),   Allocatable :: azimuth_sat
  Integer(Kind=latkind), Dimension(:,:),     Allocatable :: latitude
  Integer(Kind=latkind), Dimension(:,:),     Allocatable :: longitude
!  Integer(Kind=aemkind), Dimension(:,:),     Allocatable :: aer_mod
  Integer(Kind=aemkind), Dimension(:,:),     Allocatable :: aer_mod_1
  Integer(Kind=aemkind), Dimension(:,:),     Allocatable :: aer_mod_2
  Integer(Kind=aemkind), Dimension(:,:),     Allocatable :: aer_mod_3
  Integer(Kind=taokind), Dimension(:,:),     Allocatable :: tot_aod
  Integer(Kind=bdwkind), Dimension(:,:),     Allocatable :: bdw_mod_2
  Integer(Kind=bdwkind), Dimension(:,:),     Allocatable :: bdw_mod_3
  Integer(Kind=cstkind), Dimension(:,:),     Allocatable :: coast_pix
  Integer(Kind=spdkind), Dimension(:,:,:),   Allocatable :: wind_speed
  Integer(Kind=angkind), Dimension(:,:,:),   Allocatable :: wind_dir

  ! in/output data variables
!!$  Integer, Parameter :: parkind = 2
!!$  Integer, Parameter :: covkind = 2
!!$  Integer, Parameter :: quakind = 1
!!$  Integer, Parameter :: agekind = 1
  Integer(Kind=parkind), Dimension(:,:,:,:,:),       Allocatable :: brdf_in
  Integer(Kind=parkind), Dimension(:,:,:,:,:),       Allocatable :: brdf_in_
  Integer(Kind=covkind), Dimension(:,:,:,:,:,:),     Allocatable :: covariance_in
  Integer(Kind=quakind), Dimension(:,:,:,:),         Allocatable :: quality_in
  Integer(Kind=agekind), Dimension(:,:,:,:),         Allocatable :: age_obs_in

  Integer :: age_max = (256**agekind)/2 - 1

  ! output data variables
!!$  Integer, Parameter :: albkind = 2
!!$  Integer, Parameter :: outkind = 2
  Integer(Kind=parkind), Dimension(:,:,:,:),   Allocatable :: brdf
  Integer(Kind=parkind), Dimension(:,:,:),     Allocatable :: aodins
  Integer(Kind=parkind), Dimension(:,:),       Allocatable :: xiins
  Integer(Kind=parkind), Dimension(:,:),       Allocatable :: refins
  Integer(Kind=parkind), Dimension(:,:),       Allocatable :: jacoAODins
  Integer(Kind=covkind), Dimension(:,:,:,:,:), Allocatable :: covariance
  Integer(Kind=agekind), Dimension(:,:,:),     Allocatable :: age_obs
  Integer(Kind=agekind), Dimension(:,:,:),     Allocatable :: age_obs_al
  Integer(Kind=agekind), Dimension(:,:,:),     Allocatable :: CMins
  Integer(Kind=quakind), Dimension(:,:,:),     Allocatable :: quality
  Integer(Kind=quakind), Dimension(:,:,:),     Allocatable :: quality_al
  Integer(Kind=albkind), Dimension(:,:,:),     Allocatable :: albedo_sdh
  Integer(Kind=albkind), Dimension(:,:,:),     Allocatable :: sigma_albedo_sdh
  Integer(Kind=albkind), Dimension(:,:,:),     Allocatable :: albedo_sbh
  Integer(Kind=albkind), Dimension(:,:,:),     Allocatable :: sigma_albedo_sbh
  Integer(Kind=albkind), Dimension(:,:),       Allocatable :: albedo_bdh
  Integer(Kind=albkind), Dimension(:,:),       Allocatable :: sigma_albedo_bdh
  Integer(Kind=albkind), Dimension(:,:),       Allocatable :: albedo_bbh
  Integer(Kind=albkind), Dimension(:,:),       Allocatable :: sigma_albedo_bbh
  Integer(Kind=albkind), Dimension(:,:),       Allocatable :: albedo_vdh
  Integer(Kind=albkind), Dimension(:,:),       Allocatable :: sigma_albedo_vdh
  Integer(Kind=albkind), Dimension(:,:),       Allocatable :: albedo_ndh
  Integer(Kind=albkind), Dimension(:,:),       Allocatable :: sigma_albedo_ndh
  
  Integer(Kind=parkind), Dimension(:,:,:,:),   Allocatable :: brdf1
  Integer(Kind=covkind), Dimension(:,:,:,:,:), Allocatable :: covariance1
  Integer(Kind=quakind), Dimension(:,:,:),     Allocatable :: quality1
  
  ! pixel time series variables
  Integer, Dimension(:), Allocatable :: mask_cma
  Logical, Dimension(:), Allocatable :: ang_ok
  Logical, Dimension(:), Allocatable :: w_ok
  Logical, Dimension(:), Allocatable :: processed
  Logical, Dimension(:), Allocatable :: cloudy
  Logical, Dimension(:), Allocatable :: snowy
  Logical, Dimension(:), Allocatable :: bad_cma
  Logical, Dimension(:), Allocatable :: valid
  Logical, Dimension(:), Allocatable :: valid_near
  Real,    Dimension(:), Allocatable :: refl
  Real,    Dimension(:), Allocatable :: sigrefl
  Real,    Dimension(:), Allocatable :: sigrefl_
  Real,    Dimension(:), Allocatable :: wi_angular
  Real,    Dimension(:), Allocatable :: wi_angular_
  Real,    Dimension(:), Allocatable :: scat_ang
  Real,    Dimension(:), Allocatable :: sunglint
  Real,    Dimension(:), Allocatable :: theta_sat
  Real,    Dimension(:), Allocatable :: theta_sol
  Real,    Dimension(:), Allocatable :: t_sat_rel
  Real,    Dimension(:), Allocatable :: t_sol_rel
  Real,    Dimension(:), Allocatable :: phi_sat
  Real,    Dimension(:), Allocatable :: phi_sol
  Real,    Dimension(:), Allocatable :: phi_del
  Real,    Dimension(:), Allocatable :: wspeed
  Real,    Dimension(:), Allocatable :: wdir

  ! variables for linear model inversion
  Real, Dimension(:,:), Allocatable :: A
  Real, Dimension(:,:), Allocatable :: AT
  Real, Dimension(:),   Allocatable :: b, b_

  ! variables for communication with HDF-routines
!!$  Character(LEN=3)                                      :: O_QUA_FLAG = "OK"
  Integer, Parameter                                    :: IO_OK   = 0
  Integer                                               :: file_id
!!$  ! input
!!$  Integer(Kind=inkind), Dimension(:,:), Allocatable     :: data_r
!!$  Type(general_attributes)                              :: g_attr_r
!!$  Type(dataset_attributes)                              :: d_attr_r
!!$  ! i/o binary
!!$  Integer(Kind=1), Dimension(:,:), Allocatable          :: data_b1
!!$  Integer(Kind=2), Dimension(:,:), Allocatable          :: data_b2
!!$  ! output
!!$  Integer                                               :: n_data_w
!!$  Integer(Kind=outkind), Dimension(:,:,:), Allocatable  :: data_w
!!$  Type(general_attributes)                              :: g_attr_w
!!$  Type(dataset_attributes), Dimension(:), Allocatable   :: d_attr_w

!  ! Below, this module includes the definition of the variables 
!  ! to be read from the algorithm configuration file via namelists.
!  ! Default values can be specified which are used if the
!  ! respective variables do not appear in the configuration file.
!
!  ! NAM_PROC
!  Integer :: LinesBlock                 ! number of lines processed in one "Block"
!  Integer :: Unit                       ! Unit number for data files
!  Integer :: LogLinFrequ                ! frequency of log messages in line numbers
!  Integer :: Compression                ! 0=not compressed 1=compressed
!  Character(LEN=255) :: Path_tmp        ! path for temporary files
!  Character(LEN=1)  :: Processing_Mode  ! "N"=Nominal, "B"=Backlog, "R"=Reprocessing, "V"=Validation
!  Character(LEN=1)  :: Disposition_Flag ! "T"=Testing, "O"=Operational, "C"=Commissioning
!
!  ! NAM_SCALE scale factors and limits for storage in integer variables
!  ! note that scale factors and limits need to be consistent with Kind-values
!  Real :: scale_par  ! model parameters
!  Real :: scale_cov  ! covariance matrix
!  Real :: scale_alb  ! albedo
!  Real :: scale_sig  ! albedo error
!  Real :: par_max    ! maximum value for kernel parameters
!  Real :: par_min    ! minimum value for kernel parameters
!  Real :: cxx_max    ! maximum value for diagonal covariance matrix elements
!  Real :: cxx_min    ! minimum value for diagonal covariance matrix elements
!  Real :: cxy_max    ! maximum value for other covariance matrix elements
!  Real :: cxy_min    ! minimum value for other covariance matrix elements
!  Real :: alb_max    ! maximum value for albedo
!  Real :: alb_min    ! minimum value for albedo
!  Real :: sig_max    ! maximum value for albedo error
!  Real :: sig_min    ! minimum value for albedo error
! 
!  Integer :: MissingValue      ! integer value assigned for data points not calculated
!  Integer :: MissingValue_kcov ! idem for brdf and covariance matrix
!
!  ! NAM_INV
!  Integer :: model             ! 0=Roujean et al. 1="LiRossHotspot"
!  Logical :: recursion         ! .true. = generate recursive result based on previous estimate
!  Logical :: startseries       ! .true. = start the recursive sequence
!  Logical :: composition       ! .true. = generate one-day composition result
!  Logical :: bad_CMa_elim      ! .true. = eliminate reflectances with bad cloud mask quality
!  Real    :: bad_CMa_factor    ! penalisation factor for observations with bad cloud mask quality
!  Integer :: N_slot_elim       ! number of slots eliminated next to cloudy slots
!  Integer :: N_obs_limit       ! minimal number of observations
!  Logical :: snow_flag_one     ! .true. = one "snowy" slot suffices to set snow flag
!  Real    :: timescale         ! characteristic time scale of temporal composition (in days)
!  Real :: theta_sat_limit      ! limit for satellite zenith angles (in degrees)
!  Real :: theta_sol_limit      ! limit for solar zenith angles (in degrees)
!  Real :: theta_sat_wlimit     ! limit for weighting equation (in degrees)
!  Real :: theta_sol_wlimit     ! limit for weighting equation (in degrees)
!  Real :: theta_ref_dh_limit   ! limit for dir.-hem. reference angle (in degrees)
!  Real :: theta_sol_midi_limit ! solar zenith angle limit for re-initialisation
!  ! regularisation for model parameters (same values for all channels)
!  Real, Dimension(0:MM)         :: k_reg
!  Real, Dimension(0:MM)         :: sig_k_reg
!  ! reflectance error estimates
!  Real, Dimension(1:N_channels) :: sig_nadir_a
!  Real, Dimension(1:N_channels) :: sig_nadir_b
!  Real :: sigrefl_min          ! lower limit for reflectance error estimates
!  Real :: sigrefl_max          ! upper limit for reflectance error estimates
!  ! spectral to band conversion relative error
!  Real :: sigma_co
!  ! spectral to band conversion coefficients
!  Real, Dimension(0:N_channels) :: co_bb  ! 300-4000 broadband
!  Real, Dimension(0:N_channels) :: co_vi  ! 400- 700 visible
!  Real, Dimension(0:N_channels) :: co_ni  ! 700-4000 NIR/SWIR
  
End Module algoconf
