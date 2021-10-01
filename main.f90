PROGRAM stokes2

  USE prms
  USE data
  USE rkmod
  USE membrane
  USE outmod
  USE icmod
  USE areamod
  USE pme
  USE expint

  IMPLICIT none


  call init
!  call testQflow; stop
!  call testNb2Nf; stop
!  call testderNbb; stop
!  call testMp; stop
!  call testb; stop
!  call testFFT; stop
  call solve
  call final

CONTAINS

  SUBROUTINE init

    call initEi!; call testEi; stop
    call initprms
    call initdata
    call initRK
    call initmembrane
    call ic
    call initout(lt0)
!    call initQflow; print *,"INITIALIZED QFLOW"
    call initsurfor; print *,"SURFOR DONE"
    call nlister(X,lt0*dt); print *,"NLISTER DONE"
    call initrhs; print *,"INITIALIZED RHS"
    call mknlist(X,lt0*dt); print *,"MADE INITIAL NBRLIST"
    call initPME; print *,"INITIALIZED PME"
    call output(lt0,time0); print *,"IC OUTPUT"

  END SUBROUTINE init

  SUBROUTINE solve

    real     :: time
    integer  :: lt,k,l


    do lt = lt0+1,lt0+Nt

       print *,lt
       time = time0 + REAL(lt-1-lt0)*dt

       if (Nti.gt.0.and.MOD(lt,Nti).eq.0) then
          print *,"SHOCK SHOCK SHOCK SHOCK"
          Umax = 0.
          X(1,:) = X(1,:) + Ueps_in*COS(4.*Pi*X(2,:)/Lb(1))
!          X(1,:) = X(1,:) + Ueps_in*(X(2,:)-Lb(2)/2)/Lb(2)
       end if

       call rk2(time,lt)

       call adjustareas(X)
       call correctareas(X)

       if (MOD(lt,inlist).eq.0) then
          call nlister(X,time)
          call mknlist(X,time)
       end if

       call output(lt,time)

    end do

  END SUBROUTINE solve

  SUBROUTINE final

    call finaldata
    call finalmembrane
    call finalout

  END SUBROUTINE final

  SUBROUTINE testFFT

    real, dimension(NcP(1),NcP(2),2) :: QQ
    real     :: xx,yy
    integer  :: i,j

    do j = 1,NcP(2)
       yy = Lb(2)*REAL(j-1)/REAL(NcP(2))
       do i = 1,NcP(1)
          xx = Lb(1)*REAL(i-1)/REAL(NcP(1))
!          QVm(i,j,:) = SIN(2.*Pi*xx/Lb) + 1. + SIN(4.*Pi*xx/Lb + 4.*Pi*yy/Lb)
          QVm(i,j,:) = EXP(SIN(4.*Pi*xx/Lb(1) + 4.*Pi*yy/Lb(2)))
!          QQ(i,j,:) = 8.*Pi**2/Lb**2*EXP(SIN(4.*Pi*xx/Lb(1) + 4.*Pi*yy/Lb(2)))*(1.+COS(8.*Pi*(xx+yy)/Lb)-2.*SIN(4.*Pi*(xx/Lb(1)+yy/Lb(2))))
!          QQ(i,j,:) = - 16.*Pi**2/Lb**2*SIN(4.*Pi*xx/Lb + 4.*Pi*yy/Lb)
       end do
    end do

!   QQ = QVm

    call poissonG

    !print *,MAXVAL((QVm-QQ)**2)

    do i = 1,NcP(1)
       write(42,*)i,QQ(i,1,1)
    end do

    do j = 1,NcP(2)
       write(50,"(16F8.4)")(QQ(i,j,2),i=1,NcP(1))
    end do
    do j = 1,NcP(2)
       write(51,"(16F8.4)")(QVm(i,j,2),i=1,NcP(1))
    end do

  END SUBROUTINE testFFT


  SUBROUTINE testb

    integer :: i
    integer, parameter :: N = 32
    integer, parameter :: p = 6

    complex, dimension(-N/2:N/2-1) :: b
    complex                        :: a

    real               :: u
    integer            :: k,m

!    call makeb(p,N,b)

    do k = -N/2,N/2-1
       write(61,"(I5,2F20.10)")k,b(k)
    end do

    k = 13

    do i = 1,1001
       u = -20. + 40.*REAL(i-1)/1000
       a = 0.
       do m = INT(u)-p,INT(u)
          a = a + Mp(u-REAL(m),p)*EXP(2.*Pi*(0.,1.)*REAL(k*m)/REAL(N))
       end do
       write(60,"(5F20.10)")u,EXP(2.*Pi*(0.,1.)*REAL(k)*u/REAL(N)),b(k)*a
    end do


  END SUBROUTINE testb


  SUBROUTINE testMp

    integer :: i
    real    :: u,b

    do i = 1,101
       u = -5. + 10.*REAL(i-1)/100
       b = Mp(u,3)
       write(60,*)u,b
    end do

  END SUBROUTINE testMp


  SUBROUTINE testder

    real, dimension(Nb)  :: f,df
    integer              :: i
    real                 :: th,l0

    l0 = 2.*Pi
    do i = 1,Nb
       th = 2.*Pi*REAL(i-1)/REAL(Nb)
       f(i) = EXP(COS(th))
    end do

    call dersp(f,l0,df)

    do i = 1,Nb
       th = 2.*Pi*REAL(i-1)/REAL(Nb)
       write(2,"(I4,3F20.10)")i,f(i),df(i),-SIN(th)*f(i)
    end do

  END SUBROUTINE testder


!!$  SUBROUTINE testderNbb
!!$
!!$    real, dimension(Nbb)  :: f,df
!!$    integer              :: i
!!$    real                 :: th,l0
!!$
!!$    l0 = 2.*Pi
!!$    do i = 1,Nbb
!!$       th = 2.*Pi*REAL(i-1)/REAL(Nbb)
!!$       f(i) = EXP(COS(th))
!!$    end do
!!$
!!$    call derspNbb(f,l0,df)
!!$
!!$    do i = 1,Nbb
!!$       th = 2.*Pi*REAL(i-1)/REAL(Nbb)
!!$       write(20,"(I4,3F20.10)")i,f(i),df(i),-SIN(th)*f(i)
!!$    end do
!!$
!!$  END SUBROUTINE testderNbb


END PROGRAM stokes2
