MODULE areamod

  USE prms
  USE data
  USE membrane

  IMPLICIT none

CONTAINS

  SUBROUTINE adjustareas(X)
    real, dimension(2,Nb,Nm)    :: X
    real, dimension(Nm)      :: A

    real                     :: Atot

    real, parameter          :: fac=1.0001

!!$    call areas(X,A)
!!$
!!$    Atot = SUM(A)
!!$
!!$    if (Atot/tau0.lt.0.55) then
!!$       print *,"grow",Atot/tau0
!!$       dSo = dSo*fac
!!$       Soe = Soe*fac
!!$       len0 = len0*fac
!!$       rad = rad*fac
!!$    else
!!$       print *,Atot/tau0
!!$    end if

  END SUBROUTINE adjustareas

  SUBROUTINE areas(X,A)
    real, dimension(2,Nb,Nm)    :: X
    real, dimension(Nm)      :: A

    real, dimension(Nb) :: x1,x2,dy
    integer :: m

    do m = 1,Nm
       x1 = X(1,:,m)
       x2 = X(2,:,m)

       call dersp(x2,len0(m),dy)

       A(m) = SUM(x1*dy)*dSo((m-1)*Nb + 1)
    end do

  END SUBROUTINE areas

  SUBROUTINE correctareas(X)
    real, dimension(2,Nb,Nm) :: X

    real, dimension(Nb)      :: dx1,dy1,D,x1,x2
    real, dimension(2,Nb)    :: nvec
    real, dimension(Nm)      :: A
    real , dimension(Nm)     :: Ao


    real, parameter :: eps = 0.0001
    integer         :: m

    call areas(X,A)
!    print *,MAXVAL(A),MINVAL(A)

    Ao = Pi*rad**2

    do m = 1,Nm
       !  if (ABS(A(m)-Ao)/Ao.gt.4.*eps) then

       x1 = X(1,:,m)
       x2 = X(2,:,m)
       call dersp(x1,len0(m),dx1)
       call dersp(x2,len0(m),dy1)

       D = SQRT(dx1**2 + dy1**2)

       nvec(1,:) = dy1/D
       nvec(2,:) = -dx1/D

       X(:,:,m) = X(:,:,m) + (Ao(m)-A(m))/SUM(D)*dSo((m-1)*Nb+1)*nvec

    end do

    call areas(X,A)
!    print *,"AFTER",MAXVAL(A),MINVAL(A)

  END SUBROUTINE correctareas
  
END MODULE areamod
