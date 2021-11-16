MODULE outmod

  USE prms
  USE data
  USE membrane
  USE stokes

  IMPLICIT none

CONTAINS

  SUBROUTINE initout(lt)
    integer               :: lt

    character(50)         :: fn

    write(fn,"('D/D.out.',I10.10)")lt
    if (D_out.gt.0) open(D_unit,file=fn)

    write(fn,"('D/th.out.',I10.10)")lt
    if (th_out.gt.0) open(th_unit,file=fn)

!!$    write(fn,"('D/def.out.',I10.10)")lt
!!$    if (def_out.gt.0) open(def_unit,file=fn)

  END SUBROUTINE initout

  SUBROUTINE finalout

    if (D_out.gt.0) close(D_unit)
    if (th_out.gt.0) close(th_unit)
!!$    if (def_out.gt.0) close(def_unit)

  END SUBROUTINE finalout

  SUBROUTINE output(lt,t)
    integer      :: lt
    real         :: t

    if (x_out.gt.0.and.MOD(lt,x_out).eq.0) then
       call xout(X,X(1,Ns+1),lt,t)
    end if
    if (s_out.gt.0.and.MOD(lt,s_out).eq.0) then
       call specout(X,lt)
    end if
    if (r_out.gt.0.and.MOD(lt,r_out).eq.0) then
       call writerestart(X,lt,t)
    end if
    if (D_out.gt.0.or.th_out.gt.0) then
       call outD(X,lt)
    end if
    if (def_out.gt.0.and.MOD(lt,def_out).eq.0) then
       call outDef(X,lt,t)
    end if
!!$    if (A_out.gt.0.and.MOD(lt,A_out).eq.0) then
!!$       call outAmat(lt)
!!$    end if

  END SUBROUTINE output

  SUBROUTINE writerestart(X,lt,t)
    real, dimension(2,Np)    :: X
    integer                  :: lt
    real                     :: t

    character(50)         :: fn
    integer               :: i,m

    real, dimension(2,Np)    :: Xpert

    Xpert(:,:) = X(:,:)

    do m=1,Nm
      do i=1,Nb
        Xpert(2,(m-1)*Nb+i) = Xpert(2,(m-1)*Nb+i) + 0.01*Sin(m*2*Pi/Nm)
        Xpert(1,(m-1)*Nb+i) = Xpert(1,(m-1)*Nb+i) - 0.01*Sin(m*2*Pi/Nm)
         !Xpert(1,(m-1)*Nb+i) = Xpert(1,(m-1)*Nb+i) + 0.001*m
      end do
    end do

    write(fn,"('./D/restart.',I10.10)")lt
    open(1,file=fn,form='UNFORMATTED')
    write(1)Np,Nb,Nm,lt
    write(1)t,rad,dSo,Soe,len0,Tc,Mc,wX,wDx,wD
    write(1)Xpert
    close(1)

    write(fn,"('./D/Newrestart.',I10.10)")lt
    open(1,file=fn,form='UNFORMATTED')
    write(1)Np,Nb,Nm,lt
    write(1)t,rad
    write(1)dSo
    write(1)Soe,len0,Tc,Mc
    write(1)wX,wDx,wD
    write(1)X
    close(1)

    print *, 'len0', len0(1),len0(10)
    print *, 'dSo', dSo(1),dSo(10)
    print *, 'Soe', Soe(1),Soe(10)

  END SUBROUTINE writerestart

  SUBROUTINE xout(X,Xw,lt,t)
    real, dimension(2,Nb,Nm) :: X
    real, dimension(2,Nw)    :: Xw
    integer                  :: lt, i, m 
    real                     :: t
    real, dimension(2,Nb,Nm) :: Xe
    character(50)            :: fn,fn2
    real, dimension(Nb)      :: sum, sum2

     Xe = X
     Xe(1,:,:) = Xe(1,:,:) - FLOOR( Xe(2,:,:)/Lb(2) )*MOD(Ueps*(t-T_sh_on),Lb(1))

     !! write cell points
     write(fn,"('./D/xeout.',I10.10)")lt
     open(1,file=fn)
     do m = 1,Nm
        do i = 1,Nb
           write(1,"(2F20.10)") Xe(:,i,m) - FLOOR(Xe(:,i,m)/Lb)*Lb
        end do
     end do
     do i = 1,Nw
        write(1,"(2F20.10)") Xw(:,i)
     end do
     close(1)

     !! write wall points
     write(fn,"('./D/wout.',I10.10)")lt
     open(1,file=fn)
     do i = 1,Nw
        write(1,"(2F20.10)") Xw(:,i)
     end do
     close(1)


    !write(fn,"('./D/qeout.',I10.10)")lt
    !open(1,file=fn)

    !do m = 1,Nm
    !      write(1,"(A)") "ZONE"
    !   do i = 1,Nb
    !      if (Abs( ( Xe(1,i,m) - Floor(Xe(1,i,m)/Lb(1))*Lb(1) ) - ( Xe(1,i-1,m) - Floor(Xe(1,i-1,m)/Lb(1))*Lb(1)) ) .gt. 10 ) Then
    !        write(1,"(A)") "ZONE"
    !      end if
    !   write(1,"(2F20.10)") Xe(:,i,m) - FLOOR(Xe(:,i,m)/Lb)*Lb
    !      if (i .eq. Nb) Then
    !        write(1,"(2F20.10)") Xe(:,1,m) - FLOOR(Xe(:,1,m)/Lb)*Lb
    !      end if
    !   end do
    !end do
    
    !write(1,"(A)") "ZONE"

    !do i = 1,Nw
    !   write(1,"(2F20.10)") Xw(:,i)
    !   if (i .eq. Int(Nw/2)) then
    !     write(1,"(A)") "ZONE"
    !   end if
    !end do
    !close(1)

  END SUBROUTINE xout

  SUBROUTINE specout(X,lt)
    real, dimension(2,Nb,Nm) :: X
    integer                  :: lt

    real, dimension(Nb/2+1,2) :: xspec
    real, dimension(Nb/2+1,2) :: sspec

    real, dimension(2,2,Nm)  :: W,Xe
    real, dimension(Nb)      :: x1

    character(50)         :: fn
    integer               :: i,m,n,k1,k2,j


    write(fn,"('D/specout.',I10.10)")lt
    open(1,file=fn)
    do m = 1,Nm
       x1 = X(1,:,m)
       call spectrumNb(x1,xspec(1,1),m)
       x1 = X(2,:,m)
       call spectrumNb(x1,xspec(1,2),m)
       do i = 2,Nb/2+1
          write(1,"(I4,2E20.10)")i-1,xspec(i,:)
       end do
    end do
    close(1)

  END SUBROUTINE specout

  SUBROUTINE outDef(X,lt,t)
    real, dimension(2,Nb,Nm)   :: X
    integer                    :: lt
    real                       :: t

    real, dimension(2,Nb,Nm)   :: Xe
    real, dimension(Nb,Nm)     :: Tout,Mout

    integer               :: m,i
    character(50)         :: fn

    real      :: xmax, xc, dx,r1,r2
!!$
!!$    m = 1
!!$    xmax = MAXVAL(X(:,2,m))
!!$    xc = SUM(X(:,1,m))/REAL(Nb)
!!$    dx = MAXVAL(X(:,1,m))-MINVAL(X(:,1,m))
!!$    write(def_unit,"(6F20.10)")REAL(lt)*dt,xmax,xc,dx

!    print *,"XD",X(1,1,1)
    call getTM(X,Tout,Mout)

    Xe = X

!    print *,"UEPS",Ueps*(t-T_sh_on)

    Xe(1,:,:) = Xe(1,:,:) - FLOOR( Xe(2,:,:)/Lb(2) )*MOD(Ueps*(t-T_sh_on),Lb(1))
!    print *,"XD",Xe(1,1,1)
    write(fn,"('./D/Defeout.',I10.10)")lt
    open(1,file=fn)
    do m = 1,Nm
       do i = 1,Nb
          write(1,"(4F20.10)") Xe(:,i,m) - FLOOR(Xe(:,i,m)/Lb)*Lb,Tout(i,m),Mout(i,m)
       end do
    end do
!!$    r1 = SQRT(1.+1.0*(t-100.)/Pi)
!!$    r2 = SQRT(1.-1.0*(t-100.)/Pi)
!!$    do i = 1,Nb
!!$       write(1,"(4F20.10)") Lb/4.+ r1*COS(2.*Pi*REAL(i-1)/REAL(Nb)),Lb/2. + r1*SIN(2.*Pi*REAL(i-1)/REAL(Nb)),0.,0.
!!$    end do
!!$    do i = 1,Nb
!!$       write(1,"(4F20.10)") 3.*Lb/4.+ r2*COS(2.*Pi*REAL(i-1)/REAL(Nb)),Lb/2. + r2*SIN(2.*Pi*REAL(i-1)/REAL(Nb)),0.,0.
!!$    end do
    close(1)

  END SUBROUTINE outDef


  SUBROUTINE outD(X,lt)
    real, dimension(2,Nb,Nm)   :: X
    integer                    :: lt

    real, dimension(2,Nb,Nm)   :: Xe
    real, dimension(2,Nb,Nm)   :: nvec
    real, dimension(Nb,Nm)     :: D

    real, dimension(2,Nm)      :: Xc,Xw
    real, dimension(2,2,Nm)    :: Iij
    real                       :: tl
    real                       :: I11,I12,I21,I22

    real                       :: det,r
    real, dimension(2)         :: ev
    real                       :: lam1,lam2
    real, dimension(Nm)        :: A,B,th
    integer                    :: i,j,m,lm

    tl = REAL(lt)*dt

    Xe = X
    call normal(Xe,D,nvec)

    Xc = 0.; Xw = 0.
    do m = 1,Nm
       do j = 1,2
          do lm = 1,Nb
             Xc(j,m) = Xc(j,m) &
                  + SUM(Xe(:,lm,m)**2)*nvec(j,lm,m)*D(lm,m)*dSo((m-1)*Nm+lm)
          end do
       end do
    end do
    do m = 1,Nm
       do j = 1,2
          do lm = 1,Nb
             Xw(j,m) = Xw(j,m) &
                  + SUM(Xe(:,lm,m)*nvec(:,lm,m))*D(lm,m)*dSo((m-1)*Nm+lm)
          end do
       end do
    end do
    Xc = Xc/Xw
    write(20,*)tl,Xc

    A = -1.
    B = Lb(1)
    do m = 1,Nm
       do i = 1,Nb
          r = SQRT( SUM( (Xe(:,i,m)-Xc(:,m))**2) )
          A(m) = MAX(A(m),r)
          B(m) = MIN(B(m),r)
       end do
    end do

    Iij = 0.
    do m = 1,Nm
       do j = 1,2
          do i = 1,2
             do lm = 1,Nb
                Iij(i,j,m) = Iij(i,j,m) &
                     + 0.25*(SUM(Xe(:,lm,m)**2)*del(i,j) - Xe(i,lm,m)*Xe(j,lm,m)) &
                           *SUM(Xe(:,lm,m)*nvec(:,i,m))*D(lm,m)*dSo((m-1)*Nm+lm)
             end do
          end do
       end do
    end do

    do m = 1,Nm
       I11 = Iij(1,1,m); I12 = Iij(1,2,m); I21 = Iij(2,1,m); I22 = Iij(2,2,m)
       det = SQRT(I11**2+4.*I12*I21-2.*I11*I22+I22**2)
       lam1 = 0.5*(I11+I22-det)
       lam2 = 0.5*(I11+I22+det)
       write(23,*)tl,lam1,lam2
       if (lam1.gt.lam2) then
          ev(1) = (I11-I22-det)/2./I21
          ev(2) = 1.
       else
          ev(1) = (I11-I22+det)/2./I21
          ev(2) = 1.
       end if
 !      th(m) = ATAN2(ev(2),ev(1))
    end do
    write(22,*)tl,0.,0.
    write(22,*)tl,ev(1),ev(2)
    write(22,*)tl,0.,0.,th(2)

    if (th_out.gt.0) write(th_unit,"(20F10.5)")tl,MAXVAL(th)
!    if (th_out.gt.0) write(th_unit,"(20F10.5)")tl,(th(m),m=1,Nm)


  END SUBROUTINE outD

END MODULE outmod
