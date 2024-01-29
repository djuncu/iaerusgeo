!
! routine adapted from Numerical Recipes in Fortran90
!
Function julday(mm,id,iyyy)
  Implicit None

  Integer, Intent(In) :: mm,id,iyyy
  Integer :: julday
  Integer, Parameter :: IGREG=15+31*(10+12*1582)
  Integer :: ja,jm,jy

  jy=iyyy
  If (jy == 0) jy=1 ! there is no year zero
  If (jy < 0) jy=jy+1
  If (mm > 2) Then
     jm=mm+1
  Else
     jy=jy-1
     jm=mm+13
  End If

  julday=int(365.25*jy)+int(30.6001*jm)+id+1720995

  If (id+31*(mm+12*iyyy) >= IGREG) Then
     ja=int(0.01*jy)
     julday=julday+2-ja+int(0.25*ja)
  End If

End Function julday
