Function invers(M)
  ! calculate the inverse of a symmetric matrix M
  Real, Dimension(:,:), Intent(In)     :: M
  Real, Dimension(Size(M,1),Size(M,1)) :: invers

  Real, Dimension(Size(M,1),Size(M,1)) :: A
  Real, Dimension(Size(M,1)) :: b
  Integer, Dimension(Size(M,1)) :: indx
  Integer :: I, J, N
  Real :: d

  Interface
     Subroutine ludcmp(a,indx,d)
       Use nrtype
       Real(SP), Dimension(:,:), Intent(INOUT) :: A
       Integer(I4B), Dimension(:), Intent(OUT) :: indx
       Real(SP), Intent(OUT) :: d
     End Subroutine ludcmp
     Subroutine lubksb(a,indx,b)
       Use nrtype
       Real(SP), Dimension(:,:), Intent(IN)   :: A
       Integer(I4B), Dimension(:), Intent(IN) :: indx
       Real(SP), Dimension(:), Intent(INOUT)  :: b
     End Subroutine lubksb
  End Interface

  A = M
  N = Size(M,1)

  ! LU decomposition
  Call ludcmp(A,indx,d)
  if (A(1,1) .eq. -111.123) then
     invers(1:N,1:N)=-111.
     goto 501 
  endif

  ! calculate the inverse
  Do I = 1, N
     b = 0.
     b(I) = 1.
     Call lubksb(A,indx,b)
     invers(:,I) = b
  End Do

  ! restore symmetry of the matrix
  Do I = 1, N
     Do J = I+1, N
        invers(I,J) = 0.5*( invers(I,J)+invers(J,I) )
        invers(J,I) = invers(I,J)
     End Do
  End Do

!!$  If (1 .eq. 2) Then 
!!$501  print *, '---singular matrix no inversion----'
!!$  EndIf
  goto 502
501 print *, '---singular matrix no inversion----'
502 continue

End Function invers
