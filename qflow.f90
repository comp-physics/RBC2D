MODULE  qflow

  USE prms

  IMPLICIT none

  integer                              :: NsQ
  real, allocatable, dimension(:)      :: Qi
  real, allocatable, dimension(:,:)    :: Xi

  real                                 :: alphaQ 

  integer                              :: Nk1,Nk2
  complex, allocatable, dimension(:,:) :: Qhat
  real , allocatable, dimension(:,:)   :: km2,km2i
  real, allocatable, dimension(:)      :: kap1,kap2

CONTAINS

  SUBROUTINE initQflow

    integer      :: l,k,k1,k2

    open(1,file='qflow.in')
    read(1,*) alphaQ
    read(1,*) Nk1,Nk2
    read(1,*) NsQ

    allocate(Qi(NsQ))
    allocate(Xi(2,NsQ))

    do l = 1,NsQ
       read(1,*)Qi(l),Xi(:,l)
    end do
    close(1)
    Xi(1,:) = Xi(1,:)*Lb(1)
    Xi(2,:) = Xi(2,:)*Lb(2)

    allocate(Qhat(-Nk1:Nk1,-Nk2:Nk2))
    allocate(kap1(-Nk1:Nk1),kap2(-Nk2:Nk2),km2(-Nk1:Nk1,-Nk2:Nk2),km2i(-Nk1:Nk1,-Nk2:Nk2))

 
    do k = -Nk1,Nk1
       kap1(k) = 2.*Pi*REAL(k)/Lb(1)
    end do
    do k = -Nk2,Nk2
       kap2(k) = 2.*Pi*REAL(k)/Lb(2)
    end do
    do k2 = -Nk2,Nk2
       do k1 = -Nk1,Nk1
          km2(k1,k2) = kap1(k1)**2 + kap2(k2)**2
       end do
    end do
    km2i = 1./km2
    km2i(0,0) = 0.

    print *,km2(Nk1,Nk2)
    write(6,"('Q FOURIER: ',E20.10)")EXP(-km2(Nk1,Nk2)/4./alphaQ)
    write(6,"('Q REAL:    ',E20.10)")MAXVAL(EXP(-Lb**2*alphaQ))

    call getQhat

  END SUBROUTINE initQflow

  SUBROUTINE testQflow

    integer, parameter     :: Nx = 200
    integer, parameter     :: Ny = 200
    integer, parameter     :: Np = Nx*Ny

    real, dimension(2,Nx,Ny)  :: Xt
    real, dimension(2,Nx,Ny)  :: u

    integer                 :: i,j,k

    do j = 1,Ny
       do i = 1,Nx
          Xt(1,i,j) = Lb(1)*REAL(i-1)/REAL(Nx-1)
          Xt(2,i,j) = Lb(2)*REAL(j-1)/REAL(Ny-1)
       end do
    end do

    call getQflow(Np,Xt,u)

    open(1,file='uout.f')
    write(1,*)Nx,Ny,2
    write(1,*)u(1,:,:),u(2,:,:)
    close(1)

    open(1,file='uout.g')
    write(1,*)Nx,Ny
    write(1,*)Xt(1,:,:),Xt(2,:,:)
    close(1)

  END SUBROUTINE testQflow


  SUBROUTINE getQhat

    integer    :: k1,k2
    integer    :: i

    Qhat = 0.

    do k2 = -Nk2,Nk2
       do k1 = -Nk1,Nk1
          do i = 1,NsQ
             Qhat(k1,k2) = Qhat(k1,k2)  &
                  + Qi(i)*EXP(-(0.,1.)*(kap1(k1)*Xi(1,i) + kap2(k2)*Xi(2,i)))
          end do
       end do
    end do

  END SUBROUTINE getQhat

  SUBROUTINE getQflow(Np,X,u)
    integer                :: Np
    real, dimension(2,Np)  :: X,u

    integer                :: i,j,k1,k2

    real                   :: r2
    real, dimension(2)     :: xij

    u = 0.
    do j = 1,Np 
       do i = 1,NsQ
          xij = X(:,j)-Xi(:,i)       
          xij = xij - NINT(xij/Lb)*Lb
          r2 = SUM(xij*xij)
          u(:,j) = u(:,j) + Qi(i)*EXP(-alphaQ*r2)*xij/r2/(2.*Pi)
       end do

       do k2 = -Nk2,Nk2
          do k1 = -Nk1,Nk1
             u(:,j) = u(:,j) - REAL(Qhat(k1,k2)* &
                  (0.,1.)*(/ kap1(k1),kap2(k2) /)*km2i(k1,k2) &
                  * EXP(-km2(k1,k2)/4./alphaQ+(0.,1.)*(kap1(k1)*X(1,j) + kap2(k2)*X(2,j))))/tau0
          end do
       end do
    end do

  END SUBROUTINE getQflow

END MODULE qflow
