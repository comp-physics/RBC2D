MODULE pme

  USE prms
  USE nlistmod
  USE data
  USE membrane

  IMPLICIT none 

  integer(8)   :: Vplan_r2c                            ! plan for forward vector transforms
  integer(8)   :: Vplan_c2r                            ! plan for inverse vector transforms
  integer(8)   :: Tplan_r2c                            ! plan for forward r2 tensor transforms
  integer(8)   :: Tplan_c2r                            ! plan for inverse r2 tensor transforms
                                                       ! NOTE THESE SHOULD HAVE THE SIZE OF
                                                       ! A GCC POINTER.  8 ON ALPHA, 4 ON INTEL
  real, allocatable, dimension(:,:,:,:)   :: Ghat      ! influence function Stokeslet
  real, allocatable, dimension(:,:,:,:,:) :: Hhat      ! influence function Stresslet
  real, allocatable, dimension(:)         :: kap1,kap2 ! all wavenumbers

  real, allocatable, dimension(:,:,:)     :: kapp      ! modified wavenumber for slant cell
  real                                    :: Utlkap    ! time of creation of modwavenumber
  real, allocatable, dimension(:,:,:,:)   :: Ghatp     ! influence function Stokeslet (slant)
  real, allocatable, dimension(:,:,:,:,:) :: Hhatp     ! influence function Stresslet (slant)

  real, allocatable, dimension(:,:,:)     :: QVm       ! mesh forces/velocities
  complex, allocatable, dimension(:,:,:)  :: FFhat,Vhat ! transformed force
 
  real, allocatable, dimension(:,:,:,:)    :: Rm       ! mesh u dot n 
  real, allocatable, dimension(:,:,:)      :: Am       ! mesh  Au
  complex, allocatable, dimension(:,:,:,:) :: Uhat     ! transformed mesh velocities
  complex, allocatable, dimension(:,:,:)   :: Ahat     ! transformed A

CONTAINS
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolation/distribution of P to M and M to P

  SUBROUTINE weights(X,Icell,Wgts)
    real,dimension(2,Np)                  :: X       ! particle positions
    integer,dimension(2,Np)               :: Icell   ! atom cell list
    real,dimension(-1:1,-1:1,Np)          :: Wgts    ! interpolation/distribution weights

    real,dimension(2)                     :: Xp      ! distance of particle from cell center
    real,dimension(2,-1:1)                :: w       ! weights for distributing charge

    real                                  ::  T2
    integer                               ::  it1,it2
    integer                               ::  i,n,nset

    Wgts = 0.
    do n = 1,Np
       Xp = (X(:,n) - FLOOR(X(:,n)/Lb)*Lb - (HP(:)*REAL(Icell(:,n)-1) + 0.5*HP(:)))/HP(:)
       w(:, 1) = 0.5*(0.5+Xp(:))*(0.5+Xp(:))
       w(:, 0) = 0.75 - Xp(:)*Xp(:)
       w(:,-1) = 0.5*(0.5-Xp(:))*(0.5-Xp(:))

       do it2 = -1,1
          T2 = w(2,it2)
          do it1 = -1,1
             Wgts(it1,it2,n) = T2*w(1,it1)
          end do
       end do

    end do

  END SUBROUTINE weights

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolation/distribution of P to M and M to P

  SUBROUTINE distribute(Icell,F,Wgts,Q,Nv)
    integer,dimension(2,Np)             :: Icell    ! atom cell list
    real,dimension(Nv,Np)               :: F        ! point forces
    real,dimension(NcP(1),NcP(2),Nv)    :: Q        ! mesh forces
    real,dimension(-1:1,-1:1,Np)        :: Wgts     ! interpolation/distribution weights
    integer                             :: Nv       ! number of variables 

    integer                             ::  it1,it2
    integer                             ::  ip1,ip2
    integer                             ::  i,j,n,nset

    Q = 0.
    do j = 1,Nv
       do n = 1,Np
          do it2 = -1,1
             ip2 = MOD(Icell(2,n)+it2-1+NcP(2),NcP(2)) + 1
             do it1 = -1,1
                ip1 = MOD(Icell(1,n)+it1-1+NcP(1),NcP(1)) + 1
                Q(ip1,ip2,j) = Q(ip1,ip2,j) + F(j,n)*Wgts(it1,it2,n)
             end do
          end do
       end do
    end do

  END SUBROUTINE distribute

  SUBROUTINE interpolate(Icell,Vm,Wgts,V,Nv)
    integer,dimension(2,Np)             :: Icell ! atom cell list
    real,dimension(Nv,Np)               :: V     ! Velocity on points
    real,dimension(NcP(1),NcP(2),Nv)    :: Vm    ! Velocity on mesh
    real,dimension(-1:1,-1:1,Np)        :: Wgts  ! interpolation/distribution weights
    integer                             :: Nv

    integer                             :: it1,it2
    integer                             :: ip1,ip2
    integer                             :: i,j,n,i1,nset

    V = 0.
    do j = 1,Nv
       do n = 1,Np
          do it2 = -1,1
             ip2 = MOD(Icell(2,n)+it2-1+NcP(2),NcP(2)) + 1
             do it1 = -1,1
                ip1 = MOD(Icell(1,n)+it1-1+NcP(1),NcP(1)) + 1
                V(j,n) = V(j,n) + Vm(ip1,ip2,j)*Wgts(it1,it2,n)
             end do
          end do
       end do
    end do

    call filtknots(Nv,V)

  END SUBROUTINE interpolate

  SUBROUTINE xtoxi(X,Xi,tl)
    real,dimension(2,Np)    :: X        ! particle positions
    real,dimension(2,Np)    :: Xi       ! transformed particle positions
    real                    :: tl       ! time local

    Xi(1,:) = X(1,:) - MOD(Ueps*(tl-T_sh_on),Lb(1))/Lb(1)*X(2,:)
    Xi(2,:) = X(2,:)

  END SUBROUTINE xtoxi

  
  SUBROUTINE computePMET(X,D,U,nvec,Au,tl)
    real,dimension(2,Np)    :: X        ! particle positions
    real,dimension(2,Np)    :: Au       ! A.u
    real,dimension(2,Np)    :: U        ! velocity
    real,dimension(2,Np)    :: nvec     ! normal
    real,dimension(Np)      :: D        ! stretch
    real                       :: tl       ! local time

    real,dimension(2,Np)    :: Xi       ! transformed particle positions
    real,dimension(2,2,Np)  :: UU       ! U dot n 

    integer,dimension(2,Np)         :: Icell  ! atom cell list
    real,dimension(-1:1,-1:1,Np)    :: Wgts   ! interpolation/distribution weights

    integer                   :: i,j,m

    do m = 1,Np
       do i = 1,2
          do j = 1,2
             UU(i,j,m) = U(i,m)*nvec(j,m)*D(m)*dSo(m)
          end do
       end do
    end do

    call xtoxi(X,Xi,tl)

    call atomcell(Xi,Icell,HP,NcP)
    call weights(Xi,Icell,Wgts)
    call distribute(Icell,UU,Wgts,Rm,4)

    call poissonH

    call interpolate(Icell,Am,Wgts,Au,2)
    
  END SUBROUTINE computePMET

  SUBROUTINE computePMEV(Xe,F,V,tl)
    real,dimension(2,Np)       :: Xe       ! particle positions
    real,dimension(2,Np)       :: F        ! forces
    real,dimension(2,Np)       :: V        ! induced velocity
    real                       :: tl       ! time local

    real,dimension(2,Np)            :: Xi     ! transformed particle positions
    integer,dimension(2,Np)         :: Icell  ! atom cell list
    real,dimension(-1:1,-1:1,Np)    :: Wgts   ! interpolation/distribution weights

    integer                   :: i,j

    call xtoxi(Xe,Xi,tl)
    call makeGhatp(tl)

    call atomcell(Xi,Icell,HP,NcP)
    call weights(Xi,Icell,Wgts)
    call distribute(Icell,F,Wgts,QVm,2)

    call poissonG

    call interpolate(Icell,QVm,Wgts,V,2)


  END SUBROUTINE computePMEV

  SUBROUTINE poissonH

    integer :: i,j,k

    ! Solve "Poisson" -- multiply by optimal influence function
    call dfftw_execute(Tplan_r2c)  !  Rm -> Uhat

    Ahat = 0.
    do k = 1,2
       do j = 1,2
          do i = 1,2
             Ahat(:,:,i) = Ahat(:,:,i) + Hhatp(:,:,i,j,k)*Uhat(:,:,j,k) 
          end do
       end do
    end do

   Ahat = Ahat*(0.,1.)
! Inverse transform 
    call dfftw_execute(Tplan_c2r)  !  Ahat -> Am

  END SUBROUTINE poissonH

  SUBROUTINE poissonG
    integer :: k1,k2

    ! Solve "Poisson" -- multiply by optimal influence function
    call dfftw_execute(Vplan_r2c)  !  Q -> FFhat
!!$
    Vhat(:,:,1) = Ghatp(:,:,1,1)*FFhat(:,:,1) + Ghatp(:,:,1,2)*FFhat(:,:,2)
    Vhat(:,:,2) = Ghatp(:,:,2,1)*FFhat(:,:,1) + Ghatp(:,:,2,2)*FFhat(:,:,2)

! Inverse transform 
    call dfftw_execute(Vplan_c2r)  !  Vhat -> Vm 

  END SUBROUTINE poissonG

  SUBROUTINE initPME

    integer :: k1,k2

    allocate(kap1(NcP(1)))
    allocate(kap2(NcP(2)))

    do k1 = 1,NcP(1)/2+1
       kap1(k1) = REAL(k1-1)/Lb(1)
    end do
    do k1 = NcP(1)/2+2,NcP(1)
       kap1(k1) = -REAL(NcP(1)-k1+1)/Lb(1)
    end do

    do k2 = 1,NcP(2)/2+1
       kap2(k2) = REAL(k2-1)/Lb(2)
    end do
    do k2 = NcP(2)/2+2,NcP(2)
       kap2(k2) = -REAL(NcP(2)-k2+1)/Lb(2)
    end do

    allocate(Ghat(NcP(1)/2+1,NcP(2),2,2))
    allocate(Hhat(NcP(1)/2+1,NcP(2),2,2,2))
    call makeGHhat

    allocate(QVm(NcP(1),NcP(2),2))
    allocate(FFhat(NcP(1)/2+1,NcP(2),2))
    allocate(Vhat(NcP(1)/2+1,NcP(2),2))

    allocate(Rm(NcP(1),NcP(2),2,2))
    allocate(Am(NcP(1),NcP(2),2))
    allocate(Uhat(NcP(1)/2+1,NcP(2),2,2))
    allocate(Ahat(NcP(1)/2+1,NcP(2),2))

    call initffts

    allocate(Ghatp(NcP(1)/2+1,NcP(2),2,2))
    allocate(Hhatp(NcP(1)/2+1,NcP(2),2,2,2))
    allocate(kapp(NcP(1)/2+1,NcP(2),2))
    Utlkap = -1.

  END SUBROUTINE initPME

  SUBROUTINE initffts

    include 'fftw3.f'

    integer, dimension(2)  :: NcPcmp 

    NcPcmp = (/ NcP(1)/2+1, NcP(2) /)

                             !  plan    ,rank, n ,howmany, in,nemb,stride,dist
    call dfftw_plan_many_dft_r2c(Vplan_r2c,2   ,NcP,2      ,QVm,NcP ,1     ,NcP(1)*NcP(2), &
        !                       out ,nemb,stride,dist  ,flags
                                FFhat,NcPcmp ,1     ,(NcP(1)/2+1)*NcP(2),FFTW_ESTIMATE)

                             !  plan    ,rank, n ,howmany, in,nemb,stride,dist
    call dfftw_plan_many_dft_c2r(Vplan_c2r,2   ,NcP,2      ,Vhat,NcPcmp ,1     ,(NcP(1)/2+1)*NcP(2), &
        !                       out ,nemb,stride,dist  ,flags
                                QVm,NcP ,1     ,NcP(1)*NcP(2),FFTW_ESTIMATE)


                             !  plan    ,rank, n ,howmany, in ,nemb,stride,dist
    call dfftw_plan_many_dft_r2c(Tplan_r2c,2   ,NcP,4      ,Rm,NcP ,1     ,NcP(1)*NcP(2), &
        !                       out ,nemb,stride,dist  ,flags
                                Uhat,NcPcmp ,1     ,(NcP(1)/2+1)*NcP(2),FFTW_ESTIMATE)

                             !  plan    ,rank, n ,howmany, in   ,nemb  ,stride,dist
    call dfftw_plan_many_dft_c2r(Tplan_c2r,2   ,NcP,2      ,Ahat,NcPcmp,1     ,(NcP(1)/2+1)*NcP(2), &
        !                       out,nemb,stride,dist  ,flags
                                Am ,NcP ,1     ,NcP(1)*NcP(2),FFTW_ESTIMATE)
  END SUBROUTINE initffts


  SUBROUTINE makeGhatp(tl)
    real    :: tl

    integer :: k1,k2
    integer :: i,j,k

    real                 :: kapM2
    real, dimension(2)   :: kap

    real                 :: bb

    complex, dimension(NcP(1))      :: b1
    complex, dimension(NcP(2))      :: b2

    call makeb(2,NcP(1),b1,kap1,Lb(1))
    call makeb(2,NcP(2),b2,kap2,Lb(2))

    if ((tl-T_sh_on)*Ueps.ne.Utlkap) then
       Utlkap = (tl-T_sh_on)*Ueps
       call makekapprime(tl)

       do k2 = 1,NcP(2)
          bb = b2(k2)*CONJG(b2(k2))
          do k1 = 1,NcP(1)/2+1
             bb = bb*b1(k1)*CONJG(b1(k1))
             kapM2 = SUM(kapp(k1,k2,:)**2)
             do i = 1,2
                do j = 1,2
                   Ghatp(k1,k2,i,j) = bb*(kapM2*del(i,j) - kapp(k1,k2,i)*kapp(k1,k2,j))/kapM2**2 &
                        *(1.+Pi*alpha*kapM2)*EXP(-Pi*alpha*kapM2)/(Pi*tau0)
                end do
             end do
             do k = 1,2
                do j = 1,2
                   do i = 1,2
                      Hhatp(k1,k2,i,j,k) = &
                           bb*(2.*(del(i,k)*kapp(k1,k2,j)+del(j,k)*kapp(k1,k2,i)+del(i,j)*kapp(k1,k2,k))/kapM2 &
                           -4.*kapp(k1,k2,i)*kapp(k1,k2,j)*kapp(k1,k2,k)/kapM2**2*(1.+Pi*alpha*kapM2)) &
                           *EXP(-Pi*alpha*kapM2)/tau0
                   end do
                end do
             end do
          end do
       end do
    end if
    Ghatp(1,1,:,:) = 0.
    Hhatp(1,1,:,:,:) = 0.

  END SUBROUTINE makeGhatp


  SUBROUTINE makekapprime(tl)
    real    :: tl   ! local time

    integer :: k1,k2

     do k2 = 1,NcP(2)
       do k1 = 1,NcP(1)/2+1
          kapp(k1,k2,1) = kap1(k1)
          kapp(k1,k2,2) = kap2(k2) - kap1(k1)*MOD(Ueps*(tl-T_sh_on),Lb(1))/Lb(1)
       end do
    end do

  END SUBROUTINE makekapprime


  SUBROUTINE makeGHhat

    real, dimension(NcP(1)/2+1,NcP(2),2,2)                :: BB
    real, dimension(NcP(1)/2+1,NcP(2),2,2,2)              :: CC

    complex, dimension(NcP(1))                            :: b1
    complex, dimension(NcP(2))                            :: b2

    real    :: bb1,bb2
    integer :: k1,k2
    integer :: i,j

    call makeBigB(BB)
    call makeBigC(CC)

    call makeb(2,NcP(1),b1,kap1,Lb(1))
    call makeb(2,NcP(2),b2,kap2,Lb(2))

    do k2 = 1,NcP(2)
       bb2 = b2(k2)*CONJG(b2(k2))
       do k1 = 1,NcP(1)/2+1
          bb1 = b1(k1)*CONJG(b1(k1))
          Ghat(k1,k2,:,:)   = bb1*bb2*BB(k1,k2,:,:)/(Pi*tau0)
          Hhat(k1,k2,:,:,:) = bb1*bb2*CC(k1,k2,:,:,:)/tau0
       !   write(70,"(2I4,4E12.2)")k1,k2,Ghat(k1,k2,:,:)
       end do
    end do
    Ghat(1,1,:,:)   = 0.
    Hhat(1,1,:,:,:) = 0.

  END SUBROUTINE makeGHhat


! kernel for stresselet integral
  SUBROUTINE makeBigC(CC)
    real, dimension(NcP(1)/2+1,NcP(2),2,2,2)  :: CC

    integer :: k1,k2
    integer :: i,j,k

    real                 :: kapM2
    real, dimension(2)   :: kap

    do k2 = 1,NcP(2)
       kap(2) = kap2(k2)
       do k1 = 1,NcP(1)/2+1
          kap(1) = kap1(k1)
          kapM2 = SUM(kap**2)
          do k = 1,2
             do j = 1,2
                do i = 1,2
                   CC(k1,k2,i,j,k) = &
                        (2.*(del(i,k)*kap(j)+del(j,k)*kap(i)+del(i,j)*kap(k))/kapM2 &
                        -4.*kap(i)*kap(j)*kap(k)/kapM2**2*(1.+Pi*alpha*kapM2)) &
                        *EXP(-Pi*alpha*kapM2)
                end do
             end do
          end do

       end do
    end do

  END SUBROUTINE makeBigC

! kernel for force convolution
  SUBROUTINE makeBigB(BB)
    real, dimension(NcP(1)/2+1,NcP(2),2,2)  :: BB

    integer :: k1,k2
    integer :: i,j

    real                 :: kapM2
    real, dimension(2)   :: kap

    do k2 = 1,NcP(2)
       kap(2) = kap2(k2)
       do k1 = 1,NcP(1)/2+1
          kap(1) = kap1(k1)
          kapM2 = SUM(kap**2)
          do i = 1,2
             do j = 1,2
                BB(k1,k2,i,j) = (kapM2*del(i,j) - kap(i)*kap(j))/kapM2**2 &
                     *(1.+Pi*alpha*kapM2)*EXP(-Pi*alpha*kapM2)
             end do
          end do

       end do
    end do

  END SUBROUTINE makeBigB

  SUBROUTINE makeb(p,N,b,kap,Lbox)
    integer                        :: p,N
    real, dimension(N)             :: kap
    complex, dimension(N)          :: b
    real                           :: Lbox
    
    integer :: k,kk,m

    b = 0.
    do k = 1,N
       kk = NINT(kap(k)*Lbox)
       do m = 0,p-2
          b(k) = b(k) + Mp(REAL(m+1),p)*EXP(2.*Pi*(0.,1.)*REAL(kk*m)/REAL(N))
       end do
       b(k) = EXP(2.*Pi*(0.,1.)*REAL(k*(p-1))/REAL(N))/b(k)
    end do

  END SUBROUTINE makeb

  FUNCTION Mp(u,p) RESULT (r)
    real       :: u,r
    integer    :: p

    integer    :: k

    if (p.le.1) stop 'bad p'

!    print *,fac(6); stop
    r = 0
    do k = 0,p
       r = r + (-1.)**k/fac(k)/fac(p-k)*MAX(u-k,0.)**(p-1)
    end do
    r = r*p

  END FUNCTION Mp
 
  FUNCTION fac(u) RESULT (fu)
    integer       :: u
    real          :: fu

    integer :: n

    fu = 1
    do n = 1,u
       fu = fu*n
    end do

  END FUNCTION fac

END MODULE pme




