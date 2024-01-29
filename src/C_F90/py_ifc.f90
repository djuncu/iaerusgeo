Module py_ifc

  ! Some parameters variables are repeated without '_' for compatibility with F2PY

  ! number of AOD value for i-retrieval
  Integer, Parameter :: nb_AODpoints = 26

  ! number of scenes
  Integer, Parameter :: NScenesMax = 96
  Integer, Parameter :: N_scenes_max = NScenesMax

  ! number of channels
  Integer, Parameter :: MaxNChannels = 3
  Integer, Parameter :: Max_Nchannels = MaxNChannels
!  Integer, Parameter :: NChannels = 1
!  Integer, Parameter :: NChannels = 3
!  Integer, Parameter :: N_channels = NChannels
  Integer :: N_channels

  ! number of history files
  Integer, Parameter :: MaxNHistFiles = 3
  Integer, Parameter :: Max_NHistFiles = MaxNHistFiles

  ! number of AOD steps
  Integer, Parameter :: N_AOD_max = 301 ! from 0.00 to 3.00 with 0.01 intervals
  Integer, Parameter :: N_AOD_min = 201 ! from 0.00 to 2.00 with 0.01 intervals

  ! number of model parameters MM=2 -> 3 model parameters
  Integer, Parameter :: MM = 3
  Integer, Parameter :: MMP1 = MM + 1
  Integer, Parameter :: MM2p1 = MM*(MM + 1)

  ! Size of spatial box for instantaneous processing
  Integer, Parameter :: boxSize = 3

  ! number of fixed input files valid for all scenes (lat, lon, aod model, total aod, coast_pixels)
  Integer, Parameter :: NFixedFiles = 5

  ! number of files valid for all scenes
  Integer, Parameter :: NFilesInAll = 1 + 2*MaxNChannels
  Integer, Parameter :: N_files_in_all = NFilesInAll

  ! number of input files per scene (angles:4, ecmwf-wind:2, reflectances:NChannels)
  Integer, Parameter :: MaxNFilesInScene = 4 + 2 + MaxNChannels
  Integer, Parameter :: Max_Nfiles_in_scene = MaxNFilesInScene
  Integer :: N_files_in_scene

  ! number of output files per channel ( recursion K012, CK, AL )
  Integer, Parameter :: NFilesOutChannel = 3
  Integer, Parameter :: N_files_out_channel = NFilesOutChannel

  ! number of output files per channel ( composition K012, CK )
  Integer, Parameter :: NFilesOutChannelCompo = 2
  Integer, Parameter :: N_files_out_channel_compo = NFilesOutChannelCompo

  ! number of output files broadband ( ALBEDO )
  Integer, Parameter :: NFilesOutBroad = 1
  Integer, Parameter :: N_files_out_broad = NFilesOutBroad

  ! number of input file names to be read from the product configuration file
!  Integer, Parameter :: N_FILEINP = NScenesMax * NFilesInScene + NFilesInAll
  Integer, Parameter :: N_FILEINP = MaxNHistFiles * NFilesInAll + NFixedFiles
  Integer, Parameter :: N_FILEIN  = NScenesMax*MaxNFilesInScene

  ! number of output file names ( recursion K012, CK, AL + ALBEDO )
  Integer, Parameter :: NFileoutRec = MaxNChannels * NFilesOutChannel + NFilesOutBroad
  Integer, Parameter :: N_fileout_rec = NFileoutRec

  ! number of instantaneous output file names
  Integer, Parameter :: NFilesoutInst = 1
  Integer, Parameter :: N_files_out_inst = NFilesoutInst

  ! number of output file names to be read from the product configuration file
  Integer, Parameter :: N_FILEOUT = NFileoutRec + MaxNChannels * NFilesOutChannelCompo
  Integer, Parameter :: N_FILEINS = NScenesMax * NFilesOutInst

  ! default value for non-assigned filenames
  Character(LEN=255), Parameter :: filename_undefined = '---'

  ! variables for reading the ProductConfigurationFile
  Integer :: MODE              ! algorithm mode, not implemented
  Logical :: OCEAN_OK          ! to process ocean pixels or not
  Logical :: LAND_OK           ! to process land pixels or not
  Logical :: COAST_OK          ! to process coast pixels or not
  Integer :: YFREQINSECS       ! frequency for calling GetStopStatus

  ! variables for reading the input files
  Integer :: NHistFiles
  Integer   (Kind=2), Dimension(1:Max_NHistFiles) :: age_hist_file
  Character(LEN=255), Dimension(1:N_FILEINP)      :: YFILEINP = filename_undefined
  Character(LEN=255), Dimension(1:N_FILEOUT)      :: YFILEOUT = filename_undefined
  Character(LEN=255), Dimension(1:N_FILEINS)      :: YFILEINS = filename_undefined

  ! variable containing the filenames for the complete scenes
  Character(LEN=255), Dimension(1:N_FILEIN ) :: YFILEIN

  ! offset of history files in input files list
  Integer :: prv_off

  Logical :: fullDisk

  ! Below, this module includes the definition of the variables
  ! to be read from the algorithm configuration file via namelists.
  ! Default values can be specified which are used if the
  ! respective variables do not appear in the configuration file.

  ! NAM_PROC
  Integer :: LinesBlock                  ! number of lines processed in one "Block"
  Integer :: Unit                        ! Unit number for data files
  Integer :: LogLinFrequ                 ! frequency of log messages in line numbers
  Integer :: Compression                 ! 0=not compressed 1=compressed
  Character(LEN=255) :: Path_tmp         ! path for temporary files
  Character(LEN=1)   :: Processing_Mode  ! "N"=Nominal, "B"=Backlog, "R"=Reprocessing, "V"=Validation
  Character(LEN=1)   :: Disposition_Flag ! "T"=Testing, "O"=Operational, "C"=Commissioning

  ! NAM_SCALE scale factors and limits for storage in integer variables
  ! note that scale factors and limits need to be consistent with Kind-values
  Real    :: scale_par            ! model parameters
  Real    :: scale_cov            ! covariance matrix
  Real    :: scale_alb            ! albedo
  Real    :: scale_sig            ! albedo error
  Real    :: par_max              ! maximum value for kernel parameters
  Real    :: par_min              ! minimum value for kernel parameters
  Real    :: cxx_max              ! maximum value for diagonal covariance matrix elements
  Real    :: cxx_min              ! minimum value for diagonal covariance matrix elements
  Real    :: cxy_max              ! maximum value for other covariance matrix elements
  Real    :: cxy_min              ! minimum value for other covariance matrix elements
  Real    :: alb_max              ! maximum value for albedo
  Real    :: alb_min              ! minimum value for albedo
  Real    :: AOD_max              ! maximum value for AOD
  Real    :: AOD_min              ! minimum value for AOD
  Real    :: sig_max              ! maximum value for albedo error
  Real    :: sig_min              ! minimum value for albedo error
  Integer :: MissingValue         ! integer value assigned for data points not calculated
  Integer :: MissingValue_kcov    ! idem for brdf and covariance matrix

  ! NAM_INV
  Integer :: model_land           ! 1=Roujean, 2=RTLS, 3=Maignan
  Integer :: model_ocean          ! 1=Cox&Munk
  Logical :: recursion            ! .true. = generate recursive result based on previous estimate
  Logical :: startseries          ! .true. = start the recursive sequence
  Logical :: composition          ! .true. = generate one-day composition result
  Logical :: instantaneous        ! .true. = instantaneous estimates (vs daily estimates)
  Logical :: instoutput           ! .true. = writes instantaneous estimations
  Logical :: bad_CMa_elim_inst    ! .true. = eliminate reflectances with bad cloud mask quality for INST
  Logical :: bad_CMa_elim_daily   ! .true. = eliminate reflectances with bad cloud mask quality for DAILY
  Real    :: bad_CMa_factor       ! penalisation factor for observations with bad cloud mask quality
  Integer :: N_slot_elim          ! number of slots eliminated next to cloudy slots
  Integer :: N_obs_limit          ! minimal number of observations
  Logical :: snow_flag_one        ! .true. = one "snowy" slot suffices to set snow flag
  Real    :: timescale            ! characteristic time scale of temporal composition (in days)
  Real    :: theta_sat_limit      ! limit for satellite zenith angles (in degrees)
  Real    :: theta_sol_limit      ! limit for solar zenith angles (in degrees)
  Real    :: theta_sat_wlimit     ! limit for weighting equation (in degrees)
  Real    :: theta_sol_wlimit     ! limit for weighting equation (in degrees)
  Real    :: theta_ref_dh_limit   ! limit for dir.-hem. reference angle (in degrees)
  Real    :: theta_sol_midi_limit ! solar zenith angle limit for re-initialisation
  Real    :: xi_limit             ! scattering angle limit
  ! regularisation for model parameters (same values for all channels)
  Real, Dimension(0:MM)         :: k_reg_land
  Real, Dimension(0:MM)         :: sig_k_reg_land
  Real, Dimension(0:MM)         :: k_reg_ocean
  Real, Dimension(0:MM)         :: sig_k_reg_ocean
  ! reflectance error estimates
!  Real, Dimension(1:NChannels) :: sig_nadir_a
!  Real, Dimension(1:NChannels) :: sig_nadir_b
  Real, Dimension(1:MaxNChannels) :: sig_nadir_a
  Real, Dimension(1:MaxNChannels) :: sig_nadir_b
  Real    :: sigrefl_min          ! lower limit for reflectance error estimates
  Real    :: sigrefl_max          ! upper limit for reflectance error estimates
  ! spectral to band conversion relative error
  Real    :: sigma_co
  ! spectral to band conversion coefficients
!  Real, Dimension(0:NChannels) :: co_bb  ! 300-4000 broadband
!  Real, Dimension(0:NChannels) :: co_vi  ! 400- 700 visible
!  Real, Dimension(0:NChannels) :: co_ni  ! 700-4000 NIR/SWIR
  Real, Dimension(0:MaxNChannels) :: co_bb  ! 300-4000 broadband
  Real, Dimension(0:MaxNChannels) :: co_vi  ! 400- 700 visible
  Real, Dimension(0:MaxNChannels) :: co_ni  ! 700-4000 NIR/SWIR

  ! scientific algorithm version
  Character(LEN=4), Parameter    :: Version = '1.0'

  ! number of columns and lines
  Integer :: MSGpixX
  Integer :: MSGpixY

  ! number of useful scenes
  Integer :: N_scenes

  ! input data variables type size
  Integer, Parameter :: refkind = 2
  Integer, Parameter :: maskind = 1
  Integer, Parameter :: angkind = 2
  Integer, Parameter :: latkind = 2
  Integer, Parameter :: inkind  = 2
  Integer, Parameter :: aemkind = 1
  Integer, Parameter :: taokind = 2
  Integer, Parameter :: bdwkind = 2
  Integer, Parameter :: spdkind = 2
  Integer, Parameter :: cstkind = 1

  ! in/output data variables type size
  Integer, Parameter :: parkind = 2
  Integer, Parameter :: covkind = 2
  Integer, Parameter :: quakind = 1
  Integer, Parameter :: agekind = 1

  ! output data variables type size
  Integer, Parameter :: albkind = 2
  Integer, Parameter :: outkind = 2

  Real               :: Scale_LAT                                   ! scale factor latitude  -> degrees
  Real               :: Offset_LAT
  Integer            :: Missing_LAT
  Real               :: Scale_LON                                   ! scale factor longitude -> degrees
  Real               :: Offset_LON
  Integer            :: Missing_LON

  Real               :: Scale_TAO                                   ! scale factor total AOD climato
  Real               :: Offset_TAO
  Integer            :: Missing_TAO
  Integer            :: Missing_AEM
  Real               :: Scale_BDW
  Real               :: Offset_BDW
  Integer            :: Missing_BDW

  Integer            :: Missing_ANG
  Integer            :: Missing_COV
  Integer            :: Missing_AGE
  Integer            :: Missing_QUA_REF
  Integer            :: Missing_QUA
  Integer            :: Missing_PAR
  Integer            :: Missing_CST
  
  Real, Dimension(1:Max_Nchannels,0:MM)       :: Scale_PAR_in          ! scale factor model parameters
  Real, Dimension(1:Max_Nchannels,0:MM)       :: Offset_PAR_in
  Integer, Dimension(1:Max_Nchannels,0:MM)    :: Missing_PAR_in

  Real, Dimension(1:Max_Nchannels,1:MM2p1)    :: Scale_COV_in          ! scale factor covariance matrix
  Real, Dimension(1:Max_Nchannels,1:MM2p1)    :: Offset_COV_in
  Integer, Dimension(1:Max_Nchannels,1:MM2p1) :: Missing_COV_in

  Real, Dimension(:), Allocatable          :: Scale_REF             ! scale factor reflectance -> 1.
  Real, Dimension(:), Allocatable          :: Offset_REF
  Integer, Dimension(:), Allocatable       :: Missing_REF
  Real, Parameter                          :: Scale_conv_REF = 100. ! conversion for scale factor [% -> 1.]

  Real, Dimension(:), Allocatable          :: Scale_SAA    ! scale factor solar azimuth angle -> degrees
  Real, Dimension(:), Allocatable          :: Offset_SAA
  Integer, Dimension(:), Allocatable       :: Missing_SAA

  Real, Dimension(:), Allocatable          :: Scale_SZA    ! scale factor solar zenith angle -> degrees
  Real, Dimension(:), Allocatable          :: Offset_SZA
  Integer, Dimension(:), Allocatable       :: Missing_SZA

  Real, Dimension(:), Allocatable          :: Scale_VAA    ! scale factor view azimuth angle -> degrees
  Real, Dimension(:), Allocatable          :: Offset_VAA
  Integer, Dimension(:), Allocatable       :: Missing_VAA

  Real, Dimension(:), Allocatable          :: Scale_VZA    ! scale factor view zenith angle -> degrees
  Real, Dimension(:), Allocatable          :: Offset_VZA
  Integer, Dimension(:), Allocatable       :: Missing_VZA

  Real, Dimension(:), Allocatable          :: Scale_SWND    ! scale factor wind speed       -> m/s
  Real, Dimension(:), Allocatable          :: Offset_SWND
  Integer, Dimension(:), Allocatable       :: Missing_SWND

  Real, Dimension(:), Allocatable          :: Scale_DWND    ! scale factor wind direction   -> rad (north)
  Real, Dimension(:), Allocatable          :: Offset_DWND
  Integer, Dimension(:), Allocatable       :: Missing_DWND

  ! variables for communication with HDF-routines
  Character(LEN=3) :: O_QUA_FLAG = "OK"

Contains
  subroutine SetFileArray(ifil, input, infile)
    integer, intent(in) :: ifil, input
    character (len=*), intent(in) :: infile
    if (input.eq.0) then
       YFILEOUT(ifil+1) = infile
    elseif (input.eq.1) then
       YFILEINP(ifil+1) = infile
    elseif (input.eq.2) then
       YFILEIN (ifil+1) = infile
    else
       YFILEINS(ifil+1) = infile
    endif
  end subroutine SetFileArray

End Module py_ifc

Subroutine reportLog(LogInfoMsg_p,LogStatus_p)
  Implicit None

  Character(LEN=*), Intent(In) :: LogInfoMsg_p
  Integer, Intent(Out) :: LogStatus_p

  Integer, Parameter :: LOG_PROCESSED_OK = 0
  Integer, Parameter :: UNABLE_TO_LOG = 1

!  Print '(2I3,TR1,A)', M, D, Trim(LogInfoMsg_p)
  Print *, Trim(LogInfoMsg_p)
  LogStatus_p = LOG_PROCESSED_OK

End Subroutine reportLog

Subroutine getStopStatus(StopStatus_p)
  Implicit None

  Integer, Intent(Out) :: StopStatus_p

  Integer, Parameter :: KEEP_EXECUTION = 0
  Integer, Parameter :: STOP_EXECUTION = 1

!  Print '(2I3,TR1,A)', M, D, 'getStopStatus'
!  Print *, 'getStopStatus'
  StopStatus_p = KEEP_EXECUTION

End Subroutine getStopStatus


Subroutine stopping
  Implicit None
End Subroutine stopping
