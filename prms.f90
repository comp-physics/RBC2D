MODULE prms

  IMPLICIT none

  integer              :: Nt       ! number of timesteps
  real                 :: dt       ! timestep

  integer              :: which_ic ! sets inital condition

  integer              :: Nb       ! number of control points
  integer              :: Nbfac    ! = Nb/Nf
  integer              :: Nm       ! number of cells
  integer              :: Nf       ! number to filter down to
  integer              :: Ns       ! number of cell points Ns = Nm*Nb
  integer              :: Nw       ! number of wall points
  integer              :: Np       ! total number of points Np = Nw + Nc


  real, dimension(2)   :: Lb       ! box size
  real                 :: alpha    ! alpha:  PME
  integer,dimension(2) :: NcP      ! PME mesh
  real,dimension(2)    :: HP       ! PME mesh spaceing
  real                 :: tau0     ! box volume
  real                 :: rc       ! short-range cutoff radius

  real                 :: T0       ! tension coeficient
  real                 :: M0       ! moment coeficient
  real                 :: T1       ! tension coeficient
  real                 :: M1       ! moment coeficient
  real                 :: mu       ! viscosity
  real                 :: lambda   ! vis inside/vis outside

  real                 :: Ueps     ! shear velocity parameter
  real                 :: Umax     ! other velocity parameter
  integer              :: Nti      ! number time steps without shear
  real                 :: T_sh_on  ! time Ueps turns on
  real                 :: Ueps_in  ! input value for shear velocity

  real                 :: ro       ! initial radius
  real                 :: r1       ! initial radius leukocyte
  real                 :: l0       ! 1/initial stretch

  integer              :: x_out    ! points output interval
  integer              :: s_out    ! x-spectrum output interval
  integer              :: r_out    ! restart output interval
  integer              :: f_out    ! surface force output interval
  integer              :: D_out    ! deformation output flag
  integer              :: th_out   ! inclination output flag
  integer              :: A_out    ! output A matrix
  integer              :: def_out  ! deformation measures

  real                 :: Sf           ! surface force strength
  real                 :: sigma        ! attractive force length scale
  real                 :: Sw           ! surface force strength (wall)
  real                 :: Sa           ! attractive force strength
  real                 :: dnl          ! distance for neighborlist
  real                 :: lfac         ! *dxfac*dx for neighborlist
  real                 :: lfacdnl      !
  integer              :: inlist       !  how often to do neighborlist


  real                 :: Pi
  real, parameter      :: Egam = 0.57721566490153286061

  real, dimension(2,2) :: del

  integer, parameter   :: D_unit = 10
  integer, parameter   :: th_unit = 11
  integer, parameter   :: def_unit = 12
  
  real  :: pert
CONTAINS

  SUBROUTINE initprms

    del(1,1) = 1.; del(1,2) = 0.
    del(2,1) = 0.; del(2,2) = 1.


    open(1,file='stokes.in')
    read(1,*) Nt
    read(1,*) dt

    read(1,*) which_ic

    read(1,*) Nm
    read(1,*) Nbfac
    read(1,*) Nf;  Nb = Nf * Nbfac; Ns = Nb*Nm
    read(1,*) Nw;  Np = Ns + Nw
    read(1,*) Lb  ;  tau0 = PRODUCT(Lb)
    read(1,*) alpha
    read(1,*) NcP;    HP = Lb/REAL(NcP)
    read(1,*) rc
    read(1,*) mu
    read(1,*) lambda

    read(1,*) Ueps_in; Ueps = 0.
    read(1,*) Umax
    read(1,*) Nti; T_sh_on = 0.

    read(1,*) ro
    read(1,*) r1
    read(1,*) l0
    read(1,*) T0
    read(1,*) M0
    read(1,*) T1
    read(1,*) M1

    read(1,*) x_out
    read(1,*) s_out
    read(1,*) r_out
    read(1,*) f_out
    read(1,*) D_out
    read(1,*) th_out
    read(1,*) A_out
    read(1,*) def_out

    read(1,*) Sf
    read(1,*) sigma
    read(1,*) Sw
    read(1,*) Sa

    read(1,*) dnl
    read(1,*) lfacdnl
    read(1,*) lfac
    read(1,*) inlist
   
    read(1,*) pert

    close(1)
 !   Lb(2) = Lb(2) + 0.4
    Pi = 4.*ATAN(1.0)
    print *,"REAL:  PME",EXP(-rc*rc/alpha)
    print *,"FOURIER:  PME",EXP(-Pi*alpha*SUM( (REAL(NcP)/Lb)**2) )


  END SUBROUTINE initprms


END MODULE prms

