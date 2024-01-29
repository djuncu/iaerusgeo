Module start_mod
  Use py_ifc

  Integer            :: StopStatus, LogStatus
  Character(LEN=255) :: LogInfoMsg

  Integer, Parameter :: PROCESS_OK        = 0
  Integer, Parameter :: UNABLE_TO_PROCESS = 1
  Integer, Parameter :: STOP_EXECUTION    = 1

  ! variables for reading the system clock
  Integer :: Count_New, Count_Rate, Count_Max, Count_Old, Count_Wait

  ! total number of pixels
  Integer :: N_pix

  ! number of blocks
  Integer :: N_Blocks, N_Blocks_F

  ! number of lines in the last block
  Integer :: LinesRest

  ! character variable definitions
  Character(LEN=6)   :: outstr
  Character(LEN=3)   :: timescale_str

  ! variables for date_and_time subroutine
  Character(LEN=8)   :: aaaammdd
  Character(LEN=6)   :: hhmmss
  Character(LEN=10)  :: hhmmssxxxx
  Character(LEN=5)   :: shhmm
  Integer, Dimension(1:8) :: dtvalues

  Integer :: year_in         = 1900 ! year
  Integer :: month_in        =    1 ! month
  Integer :: day_of_month_in =    1 ! day of the month

  Integer :: year            ! year
  Integer :: month           ! month
  Integer :: day_of_month    ! day of the month

  ! mathematical constants
  Real   , Parameter :: pi      = 3.141592653589793238462643383279502884197
  Real   , Parameter :: two_pi  = 2.*pi
  Real   , Parameter :: pi_half = pi/2.
  Real   , Parameter :: rad     = pi/180.
  Real   , Parameter :: deg     = 180./pi

  ! variables for linear model inversion
  Real   , Dimension(0:MM)      :: k_force   ! model parameters
  Real   , Dimension(0:MM,0:MM) :: Ck_force  ! covariance matrix of the k estimate

  Real   , Dimension(0:MM)      :: k         ! model parameters
  Real   , Dimension(0:MM,0:MM) :: Ck        ! covariance matrix of the k estimate
  Real   , Dimension(0:MM,0:MM) :: Ck_ite    ! covariance matrix of the k estimate
  Real   , Dimension(0:MM,0:MM) :: CkI       ! inverse of the covariance matrix of the k estimate
  Real   , Dimension(0:MM,0:MM) :: ATA       ! matrix ATA
  Real   , Dimension(0:MM)      :: k_in      ! a priori information for parameters
  Real   , Dimension(0:MM,0:MM) :: Ck_in     ! a priori info covariance matrix
  Real   , Dimension(0:MM,0:MM) :: CkI_tau   ! inverse of the covariance matrix of the k estimate

  Real   , Dimension(0:MM-1)        :: k_         ! model parameters
  Real   , Dimension(0:MM-1,0:MM-1) :: Ck_        ! covariance matrix of the k estimate
  Real   , Dimension(0:MM-1,0:MM-1) :: CkI_       ! inverse of the covariance matrix of the k estimate
  Real   , Dimension(0:MM-1)        :: k_reg_
  Real   , Dimension(0:MM)          :: k_in_      ! a priori information for parameters
  Real   , Dimension(0:MM-1,0:MM-1) :: CkI_in_    ! inverse of a priori info covariance matrix
  Real   , Dimension(0:MM-1,0:MM-1) :: CkI_reg_i_ ! inverse of initialisation covariance matrix
  Real   , Dimension(0:MM-1,0:MM-1) :: CkI_tau_   ! inverse of the covariance matrix of the k estimate
  Real   , Dimension(0:MM-1,0:MM-1) :: ATA_       ! matrix ATA
  Real   , Dimension(0:MM-1,0:MM-1) :: Ck_in_     ! a priori info covariance matrix

!***DAILY
  Real   , Dimension(0:MM,0:MM)     :: CkI_surf   ! covariance matrix of the k estimate
  Real   , Dimension(0:MM,0:MM)     :: CkI_in     ! inverse of a priori info covariance matrix
  Real   , Dimension(0:MM,0:MM)     :: CkI_reg    ! inverse of initialisation covariance matrix
  Real   , Dimension(0:MM,0:MM)     :: CkI_reg_i  ! inverse of initialisation covariance matrix
  Real   , Dimension(0:MM,0:MM)     :: Ck_1       ! covariance matrix of the k estimate
  Real   , Dimension(0:MM,0:MM)     :: sqcvm      ! square root of cov. matrix elements
  Real   , Dimension(0:MM)          :: k_sp
  Real   , Dimension(0:MM)          :: k_surf     ! model parameters

  Real   , Dimension(1:Max_Nchannels)  :: albed_d    ! auxiliary variable for dir. albedo
  Real   , Dimension(1:Max_Nchannels)  :: albed_b    ! auxiliary variable for bi-hem. albedo
  Real   , Dimension(1:Max_Nchannels)  :: sigma_d    ! auxiliary variable for albedo error
  Real   , Dimension(1:Max_Nchannels)  :: sigma_b    ! auxiliary variable for albedo error
  Logical, Dimension(1:Max_Nchannels)  :: valids     ! validity of spectral albedo estimates

  ! directional-hemispherical integrals of the model kernel functions
  Real   , Dimension(0:MM)          :: dihi_calc

  Real   , Dimension(:,:), Allocatable :: AT_
  Real   , Dimension(:,:), Allocatable :: A_
  Real   , Dimension(:)  , Allocatable :: trans

  ! interface for matrix inversion
  Interface
     Function invers(M)
       Real, Dimension(:,:), Intent(In)     :: M
       Real, Dimension(Size(M,1),Size(M,1)) :: invers
     End Function invers
  End Interface

End Module start_mod


Character(LEN=255) Function tmp_file_path(filename, extent, varkind)
  Use py_ifc
  Use algoconf

  Character(LEN=*), intent(in)   :: filename, extent
  Integer, intent(in)            :: varkind
  Character(LEN=4)               :: kindext

  kindext = '.gb2'
  if (varkind.eq.1) kindext = '.gb1'
  tmp_file_path = filename(Index(filename,"/",.true.)+1:Len(filename))
  if (Len_Trim(extent).gt.0) tmp_file_path = Trim( Trim(tmp_file_path) // '.' // Trim(extent))
  tmp_file_path = Trim(Path_tmp) // Trim(Trim(tmp_file_path) // kindext)
end Function tmp_file_path


Logical Function IsError(test, msg, ProcessingStatus_p)

  Character(LEN=*), intent(in) :: msg
  Logical, Intent(In)          :: test
  Integer, Intent(Out)         :: ProcessingStatus_p
  Integer                      :: LogStatus

  IsError = ( .not. test )
  If (IsError) Then
     Call reportLog(msg, LogStatus)
     Call stopping
     ProcessingStatus_p = UNABLE_TO_PROCESS
     Return
  End If

end Function IsError


Logical Function read_tmp_file(filename, varkind, npix, rec, linblk, maxrec, ProcessingStatus_p, data, extent, fillv)
  Use py_ifc

  Character(LEN=*), intent(in)   :: filename, extent
  Integer, intent(in)            :: varkind, linblk, npix, rec, maxrec, fillv
  Integer, Intent(Out)           :: ProcessingStatus_p
  Integer(Kind=1), Dimension(npix*varkind, 0:linblk+1), intent(out) :: data

  Character(LEN=255)             :: tfilename, tmp_file_path
  Logical read_tmp_file_bis

  tfilename = tmp_file_path(filename, extent, varkind)
  data = fillv
  read_tmp_file = read_tmp_file_bis(tfilename, varkind, linblk*npix, rec, ProcessingStatus_p, data(:,1:linblk))
  If ( fullDisk ) Then
     If (read_tmp_file .and. rec.gt.1     ) read_tmp_file = read_tmp_file_bis(tfilename, &
          & varkind, npix, (rec-1)*linblk, ProcessingStatus_p, data(:, 0))
     If (read_tmp_file .and. rec.lt.maxrec) read_tmp_file = read_tmp_file_bis(tfilename, &
          & varkind, npix,  rec*linblk+1 , ProcessingStatus_p, data(:,linblk+1))
  EndIf

End Function read_tmp_file


Logical Function read_tmp_file_bis(filepath, varkind, npix, rec, ProcessingStatus_p, data)
  Use py_ifc

  Character(LEN=*), intent(in)   :: filepath
  Integer, intent(in)            :: varkind, npix, rec
  Integer, Intent(Out)           :: ProcessingStatus_p
  Integer(Kind=1), Dimension(varkind*npix), intent(out) :: data

  Character(LEN=255)             :: LogInfoMsg
  Integer                        :: LogStatus
  Logical                        :: IsError

  read_tmp_file_bis = .False.
  Open (UNIT=Unit, IOSTAT=ios, STATUS='Old', ACTION='Read', FILE=filepath, &
       & FORM='Unformatted', ACCESS='Direct', RECL=varkind*npix)
  If (IsError(ios.eq.0, 'Error when R-opening '//Trim(filepath), ProcessingStatus_p)) Return
  Read (UNIT=Unit, IOSTAT=ios, REC=rec) data
  Close (UNIT=Unit)
  If (IsError(ios.eq.0, 'Error when reading '//Trim(filepath), ProcessingStatus_p)) Return
  read_tmp_file_bis = .True.

End Function read_tmp_file_bis


Logical Function write_tmp_file(filename, varkind, npix, status, rec, linblk, maxrec, nlines, ProcessingStatus_p, data, extent)
  Use py_ifc

  Character(LEN=*), intent(in) :: filename, status, extent
  Integer, intent(in)          :: varkind, npix, rec, linblk, nlines
  Integer, Intent(Out)         :: ProcessingStatus_p
  Integer(Kind=1), Dimension(varkind*npix, 1:linblk), intent(in) :: data

  Character(LEN=255)           :: tfilename, tmp_file_path
  Integer                      :: ilin, irec
  Logical write_tmp_file_bis

  tfilename = tmp_file_path(filename, extent, varkind)
  If (rec.lt.maxrec .or. .not. FullDisk) Then
     write_tmp_file = write_tmp_file_bis(tfilename, varkind, linblk*npix, status, rec, ProcessingStatus_p, data(:,1:linblk))
  Else
     Do ilin = 1, nlines
        irec = (rec-1)*linblk + ilin
        write_tmp_file = write_tmp_file_bis(tfilename, varkind, npix, status, irec, ProcessingStatus_p, data(:,ilin))
        if (.not. write_tmp_file) return
     End Do
  EndIf
  write_tmp_file = .True.

End Function write_tmp_file


Logical Function write_tmp_file_bis(filepath, varkind, npix, status, rec, ProcessingStatus_p, data)
  Use py_ifc

  Character(LEN=*), intent(in) :: filepath, status
  Integer, intent(in)          :: varkind, npix, rec
  Integer, Intent(Out)         :: ProcessingStatus_p
  Integer(Kind=1), Dimension(varkind*npix), intent(in) :: data

  Character(LEN=255)           :: LogInfoMsg
  Integer                      :: LogStatus
  Logical                      :: IsError

  write_tmp_file_bis = .False.
  Open (UNIT=Unit, IOSTAT=ios, STATUS=status, ACTION='Write',&
       & FILE=filepath, FORM='Unformatted',&
       & ACCESS='Direct', RECL=varkind*npix)
  If (IsError(ios.eq.0, 'Error when W-opening '//Trim(filepath), ProcessingStatus_p)) Return
  Write (UNIT=Unit, IOSTAT=ios, REC=rec) data
  Close (UNIT=Unit)
  If (IsError(ios.eq.0, 'Error when writing '//Trim(filepath), ProcessingStatus_p)) Return
  write_tmp_file_bis = .True.

End Function write_tmp_file_bis
