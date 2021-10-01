MODULE membrane

  USE prms
  USE data

  IMPLICIT none

! Nb
  integer(8)                             :: fp,rp
  real, allocatable, dimension(:)        :: fw
  complex, allocatable, dimension(:)     :: fhat

  integer(8)                             :: ffp
  real, allocatable, dimension(:)        :: ffw
  complex, allocatable, dimension(:)     :: ffh

  real, allocatable, dimension(:)        :: kb

  real, allocatable, dimension(:,:,:)    :: bend,stretch,totE

CONTAINS

  SUBROUTINE initmembrane

    allocate(fw(Nb),fhat(Nb/2+1))
    allocate(ffw(Nf),ffh(Nf/2+1))
    allocate(kb(Nb/2+1))
    allocate(bend(2,Nb,Nm),stretch(2,Nb,Nm),totE(2,Nb,Nm))

    call initdersp

  END SUBROUTINE initmembrane

  SUBROUTINE finalmembrane

 !   deallocate(fw,fhat,kb,ffw,ffh)

  END SUBROUTINE finalmembrane

  SUBROUTINE initdersp

    include 'fftw3.f'

    integer  :: n

    do n = 0,Nb/2
       kb(n+1) = 2.*Pi*REAL(n)
    end do

! Nb
    call dfftw_plan_dft_r2c_1d(fp,Nb,fw,fhat,FFTW_ESTIMATE)     ! fp:   fw(Nb) -> fhat(Nb/2+1)
    call dfftw_plan_dft_c2r_1d(rp,Nb,fhat,fw,FFTW_ESTIMATE)     ! rp:   fhat(Nb/2+1) -> fw(Nb)


    call dfftw_plan_dft_r2c_1d(ffp,Nf,ffw,ffh,FFTW_ESTIMATE)    ! ffp:  ffw(Nf) -> ffh(Nf/2+1)


  END SUBROUTINE initdersp

  SUBROUTINE spectrumNb(f,fh,m)
    real, dimension(Nb)     :: f
    real, dimension(Nb/2+1) :: fh
    integer                 :: m

    fw = f
    call dfftw_execute(fp)

    fhat = fhat/REAL(Nb)

    fh = fhat*CONJG(fhat)

  END SUBROUTINE spectrumNb

  SUBROUTINE dersp(f,ul,df)
    real, dimension(Nb) :: f,df
    real                :: ul

    integer  :: n

    fw = f
    call dfftw_execute(fp)

    fhat(Nb/2+1) = 0.
    fhat = fhat*(0.,1.)*kb

    call dfftw_execute(rp)
    df = fw/ul/REAL(Nb)

  END SUBROUTINE dersp

  SUBROUTINE filtknots(Nv,xx)
    integer                   :: Nv
    real, dimension(Nv,Nb,Nm) :: xx

    integer  :: n,m

    character(30) :: fn

    do m = 1,Nm
       do n = 1,Nv
          fw = xx(n,:,m)
          call dfftw_execute(fp)
          fhat(Nf/2+1:Nb/2+1) = 0.
          call dfftw_execute(rp)
          xx(n,:,m) = fw/REAL(Nb)
       end do
    end do

  END SUBROUTINE filtknots


  SUBROUTINE filt1(f,m)
    real, dimension(Nb) :: f

    integer  :: m

    fw = f
    call dfftw_execute(fp)
    fhat(Nf/2+2:Nb/2+1) = 0.
    call dfftw_execute(rp)
    f = fw/REAL(Nb)

  END SUBROUTINE filt1



  SUBROUTINE Nf2Nb(Xf,Xb)
    real, dimension(2,Nf,Nm) :: Xf
    real, dimension(2,Nb,Nm) :: Xb

    integer  :: i,m

    character(30) :: fn

    do m = 1,Nm
       do i = 1,2
          ffw = Xf(i,:,m)

          call dfftw_execute(ffp)
          fhat(1:Nf/2+1) = ffh
          fhat(Nf/2+2:Nb/2+1) = 0.
          call dfftw_execute(rp)
          Xb(i,:,m) = fw/REAL(Nf)
       end do
    end do

  END SUBROUTINE Nf2Nb

  SUBROUTINE normal(Xe,D,nvec)
    real, dimension(2,Nb,Nm) :: Xe,nvec
    real, dimension(Nb,Nm)   :: D

    real, dimension(Nb)      :: x1
    real, dimension(Nb)      :: dx1,dy1
    integer                  :: m

    do m = 1,Nm
       x1 = Xe(1,:,m)
       call dersp(x1,len0(m),dx1)
       x1 = Xe(2,:,m)
       call dersp(x1,len0(m),dy1)
       D(:,m) = SQRT(dx1**2 + dy1**2)

       nvec(1,:,m) = dy1/D(:,m)
       nvec(2,:,m) = -dx1/D(:,m)
    end do

  END SUBROUTINE normal

  SUBROUTINE wall(X,dv,DD,Dx)
    real, dimension(2,Np) :: X,dv,Dx,DDx2,nvec
    real, dimension(Np)   :: DD
    integer :: i

    DD(Ns+1:Np) = wD
    Dx(:,Ns+1:Np) = wDx

    dv(:,Ns+1:Np) = Sw*(X(:,Ns+1:Np)-wX)


    !open(1,file='Wall.Forces',form='FORMATTED')
    !Do i = Ns+1, Np
        !write(1,*) dv(:,i)
    !End Do

    !close(1)


  END SUBROUTINE wall

  SUBROUTINE construct_jumps(X,dv,DD,Dx,DDx2,nvec)
    real, dimension(2,Nb,Nm) :: X,dv,Dx,DDx2,nvec
    real, dimension(Nb,Nm)   :: DD

    real, dimension(Nb)   :: D,D3i,Di
    real, dimension(Nb)   :: x1,y1,dx1,dy1,dx2,dy2,dvx,dvy
    real, dimension(Nb)   :: dx1SQ,dy1SQ
    real, dimension(Nb)   :: kap,w1
    real, dimension(Nb,2) :: w2,tvec,nvecl,w3,u1,ttau

    character(20)  :: fn
    integer :: i,m

    call filtknots(2,X)

    do m = 1,Nm
       x1 = X(1,:,m)
       y1 = X(2,:,m)

       call dersp(x1,len0(m),dx1)
       call dersp(y1,len0(m),dy1)
       Dx(1,:,m) = dx1
       Dx(2,:,m) = dy1

       call dersp(dx1,len0(m),dx2)
       call dersp(dy1,len0(m),dy2)

       DDx2(1,:,m) = dx2
       DDx2(2,:,m) = dy2

       dx1SQ = dx1*dx1
       dy1SQ = dy1*dy1
       call filt1(dx1SQ,m)
       call filt1(dy1SQ,m)

       D = SQRT(dx1SQ + dy1SQ)
       call filt1(D,m)
       DD(:,m) =  D

       Di = 1./D
       call filt1(Di,m)

       D3i = Di**3
       call filt1(D3i,m)

       kap = (dy1*dx2-dx1*dy2)*D3i

       call filt1(kap,m)

       nvecl(:,1) = dy1*Di
       nvecl(:,2) = -dx1*Di

       call filt1(nvecl(1,1),m)
       call filt1(nvecl(1,2),m)

       nvec(1,:,m) = nvecl(:,1)
       nvec(2,:,m) = nvecl(:,2)

       tvec(:,1) = dx1*Di
       tvec(:,2) = dy1*Di

       call filt1(tvec(1,1),m)
       call filt1(tvec(1,2),m)

       call dersp(kap,len0(m),w1)
       w2(:,1) = nvecl(:,1)*Di*w1
       w2(:,2) = nvecl(:,2)*Di*w1

       call filt1(w2(1,1),m)
       call filt1(w2(1,2),m)

       call dersp(w2(1,1),len0(m),w3(1,1))
       call dersp(w2(1,2),len0(m),w3(1,2))
       w2(:,1) = w3(:,1)*Di
       w2(:,2) = w3(:,2)*Di

       call filt1(w2(1,1),m)
       call filt1(w2(1,2),m)

       ttau(:,1) = tvec(:,1)*(D-1.)
       ttau(:,2) = tvec(:,2)*(D-1.)
       call dersp(ttau(1,1),len0(m),u1(1,1))
       call dersp(ttau(1,2),len0(m),u1(1,2))

       call filt1(u1(1,1),m)
       call filt1(u1(1,2),m)

       dvx = -Tc(m)*u1(:,1)*Di + Mc(m)*w2(:,1)
       dvy = -Tc(m)*u1(:,2)*Di + Mc(m)*w2(:,2)


!!$       dvx = -T0*nvec(1,:,m)*kap
!!$       dvy = -T0*nvec(2,:,m)*kap

       call filt1(dvx,m)
       call filt1(dvy,m)

       bend(1,:,m) = Mc(m)*w2(:,1)
       bend(2,:,m) = Mc(m)*w2(:,2)
       stretch(1,:,m) = -Tc(m)*u1(:,1)*Di
       stretch(2,:,m) = -Tc(m)*u1(:,2)*Di
       totE(1,:,m) = -Tc(m)*u1(:,1)*Di + Mc(m)*w2(:,1)
       totE(2,:,m) = -Tc(m)*u1(:,2)*Di + Mc(m)*w2(:,2)

       dv(1,:,m) = dvx
       dv(2,:,m) = dvy


    end do


    !print *,"DD",MAXVAL(DD),MINVAL(DD)


  END SUBROUTINE construct_jumps

  SUBROUTINE getTM(X,Tout,Mout)
    real, dimension(2,Nb,Nm) :: X
    real, dimension(Nb,Nm)   :: Tout,Mout

    real, dimension(Nb)   :: D,D3i,Di
    real, dimension(Nb)   :: x1,y1,dx1,dy1,dx2,dy2
    real, dimension(Nb)   :: dx1SQ,dy1SQ
    real, dimension(Nb)   :: kap

    character(20)  :: fn
    integer :: i,m

    call filtknots(2,X)

    do m = 1,Nm
       x1 = X(1,:,m)
       y1 = X(2,:,m)

       call dersp(x1,len0(m),dx1)
       call dersp(y1,len0(m),dy1)

       call dersp(dx1,len0(m),dx2)
       call dersp(dy1,len0(m),dy2)

       dx1SQ = dx1*dx1
       dy1SQ = dy1*dy1
       call filt1(dx1SQ,m)
       call filt1(dy1SQ,m)

       D = SQRT(dx1SQ + dy1SQ)
       call filt1(D,m)

       Di = 1./D
       call filt1(Di,m)

       D3i = Di**3
       call filt1(D3i,m)

       kap = (dy1*dx2-dx1*dy2)*D3i

       Tout(:,m) = Tc(m)*(D-1)
       Mout(:,m) = Mc(m)*kap

    end do


  END SUBROUTINE getTM

END MODULE membrane
