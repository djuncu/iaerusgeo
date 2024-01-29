Subroutine smoothing(status, AOD, n, CM, m, QFLG, p, coast_in, q, nx, ny, win_size1, win_size2, scl, AOD_min, AOD_max)
  implicit none

  Integer        , Intent(out)                 :: status
  Integer(kind=2), intent(inout), Dimension(n) :: AOD
  Integer(kind=1), intent(in)   , Dimension(m) :: CM
  Integer(kind=1), intent(in)   , Dimension(p) :: QFLG
  Integer(kind=1), intent(in)   , Dimension(q) :: coast_in
  Integer        , intent(in)                  :: n, m, p, q, nx, ny, win_size1, win_size2
  Real           , intent(in)                  :: scl, AOD_min, AOD_max

  Real           , Dimension(:,:), Allocatable :: AOD_0, AOD_s1, AOD_s2, AOD_s3, AOD_s4
  Logical        , Dimension(:,:), Allocatable :: mask
  Integer(kind=1), Dimension(:,:), Allocatable :: coast_pix, CM_0, QFLG_0
  Real           , Dimension(:)  , Allocatable :: aod_w
  Logical        , Dimension(:)  , Allocatable :: msk_w, bad_w

  Integer  :: x, y, rsize1, rsize2, rsize_cst, thrshld2 = 3, thrshld_cst = 9, astat, win_size_cst
  Real     :: max_w, mean_w, thrshld1 = 1500. ! 0.15 but note that it was 0.2 for 550nm in Lyasputin et al. (2018)

  ! *****************************************************
  Integer, dimension(3):: ligs, cols
  Integer i, l, c
  data ligs /1300, 1300, 3050/, cols /1300, 3000, 1300/
  ! *****************************************************

  win_size_cst = win_size2 * 3.0 ! e.g., win_size2 = 3 -> win_size_cst = 9

  Allocate ( AOD_0(1:nx,1:ny), AOD_s1(1:nx,1:ny), AOD_s2(1:nx,1:ny), mask(1:nx,1:ny), &
       &     coast_pix(1:nx,1:ny), aod_w(1:win_size1*win_size1),                      &
       &     msk_w(1:win_size1*win_size1), bad_w(1:win_size1*win_size1), STAT=status)
  if (status .ne. 0) then
     print *, 'Error when allocating variables'
     return
  endif

  AOD_0     = reshape(AOD     , (/nx,ny/))/scl
  CM_0      = reshape(CM      , (/nx,ny/))/scl
  QFLG_0    = reshape(QFLG    , (/nx,ny/))/scl
  coast_pix = reshape(coast_in, (/nx,ny/))/scl
  mask      = ((AOD_0 .gt. AOD_min) .and. (AOD_0 .lt. AOD_max))

! 1st filter: eliminating isolated pixels with high AOD (Emili et al., 2011)
  AOD_s1 = AOD_0
  rsize1 = (win_size1-1)/2
  Do Y = 1+rsize1, ny-rsize1
     Do X = 1+rsize1, nx-rsize1
        if (.not. mask(x, y)) cycle
        aod_w = RESHAPE(AOD_0(x-rsize1:x+rsize1, y-rsize1:y+rsize1),(/win_size1*win_size1/))
        msk_w = ((aod_w .gt. AOD_min) .and. (aod_w .lt. AOD_max))
        max_w = MAXVAL(aod_w, msk_w)
        bad_w = (msk_w .and. (aod_w .eq. max_w))
        where(bad_w) msk_w = .False.
        mean_w = SUM(aod_w, msk_w)/COUNT(msk_w)
        If (max_w .gt. (mean_w + thrshld1)) Then
           where(bad_w) aod_w = mean_w
           AOD_s1(X-rsize1:X+rsize1,Y-rsize1:Y+rsize1) = RESHAPE(aod_w,(/win_size1,win_size1/))
        EndIf
     End Do
  End Do

! 2nd filter: running averaging window over coast pixels to compensate the lack of pixels for smoothing over these regions
  AOD_s2 = AOD_s1
  mask = ((AOD_s1 .gt. AOD_min) .and. (AOD_s1 .lt. AOD_max))
  rsize_cst = (win_size_cst-1)/2
  Do Y = 1+rsize_cst, ny-rsize_cst
     Do X = 1+rsize_cst, nx-rsize_cst
        if (coast_pix(X, Y) .ne. 1) cycle
        aod_w = RESHAPE(AOD_s1(x-rsize_cst:x+rsize_cst, y-rsize_cst:y+rsize_cst),(/win_size_cst*win_size_cst/))
        msk_w = ((aod_w .gt. AOD_min) .and. (aod_w .lt. AOD_max))
        If (COUNT(msk_w) .ge. thrshld_cst) Then ! if at least thrshld_cst pixels are good
           mean_w = SUM(aod_w, msk_w)/COUNT(msk_w)
           AOD_s2(x,y) = mean_w
        Endif
     End Do
  End Do

! 3rd filter: running averaging window (Lyapustin et al. 2018)
  AOD_s3 = AOD_s2
  mask = ((AOD_s2 .gt. AOD_min) .and. (AOD_s2 .lt. AOD_max))
  rsize2 = (win_size2-1)/2
  Do Y = 1+rsize2, ny-rsize2
     Do X = 1+rsize2, nx-rsize2
        if (.not. mask(x, y)) cycle
        aod_w = RESHAPE(AOD_s2(x-rsize2:x+rsize2, y-rsize2:y+rsize2),(/win_size1*win_size1/))
        msk_w = ((aod_w .gt. AOD_min) .and. (aod_w .lt. AOD_max))
        If (COUNT(msk_w) .ge. thrshld2) Then ! if at least thrshld2 pixels are good
           mean_w = SUM(aod_w, msk_w)/COUNT(msk_w)
           AOD_s3(x,y) = mean_w
        Endif
     End Do
  End Do

! 4th filter: filling NaN values
  AOD_s4 = AOD_s3
  mask = ((AOD_s3 .gt. AOD_min) .and. (AOD_s3 .lt. AOD_max))
  rsize2 = (win_size2-1)/2
  Do Y = 1+rsize2, ny-rsize2
     Do X = 1+rsize2, nx-rsize2
        if (mask(x, y)) cycle
        aod_w = RESHAPE(AOD_s3(x-rsize2:x+rsize2, y-rsize2:y+rsize2),(/win_size1*win_size1/))
        msk_w = ((aod_w .gt. AOD_min) .and. (aod_w .lt. AOD_max))
        If (COUNT(msk_w) .ge. thrshld2) Then ! if at least thrshld2 pixels are good
           mean_w = SUM(aod_w, msk_w)/COUNT(msk_w)
           AOD_s4(x,y) = mean_w
        Endif
     End Do
  End Do

! Update AOD
  AOD = RESHAPE(nint(AOD_s4*scl), (/nx*ny/))

end Subroutine smoothing

