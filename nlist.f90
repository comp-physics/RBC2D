MODULE nlistmod

  USE prms
  USE data

  IMPLICIT none

  TYPE :: listelm
    integer                 :: ll    ! point index of neighbor
    integer                 :: mm    ! cell index of neighbor
  END TYPE listelm

CONTAINS

  SUBROUTINE atomcell(X,Icell,hspace,Nmax)
    real,dimension(2,Np)          :: X        ! particle positions
    integer,dimension(2,Np)       :: Icell    ! atom cell list
    real,dimension(2)             :: hspace   ! mesh cell spacing for distribution
    integer,dimension(2)          :: Nmax     ! max cell
    integer                       :: l

    do l = 1,2
       Icell(l,:) = INT( (X(l,:) - FLOOR(X(l,:)/Lb(l))*Lb(l))/hspace(l)) + 1
    end do

  END SUBROUTINE atomcell


  SUBROUTINE chainlist(LL,HOC,N,X,tl)
    integer,dimension(Np)                            :: LL    ! next in list pointer
    integer,dimension(2)                             :: N     ! number of cells
    integer,dimension(0:N(1)+1,0:N(2)+1)             :: HOC   ! Head-of-Chain pointer
    ! with perioidic continuation
    real   ,dimension(2,Np)                          :: X     ! particle positions
    real                                             :: tl    ! time

    integer,dimension(2,Np)                          :: Icell ! atom cell list
    integer                                          :: i
    real,dimension(2)                                :: HMf   ! cell spacing

    HMf = Lb/REAL(N)
    call atomcell(X,Icell,HMf,N)

    HOC = 0
    do i = 1,Np
       LL(i) = HOC(Icell(1,i),Icell(2,i))
       HOC(Icell(1,i),Icell(2,i)) = i
    end do

    call makeper(HOC,N,tl)

  END SUBROUTINE chainlist


  SUBROUTINE makeper(HOC,N,tl)
    integer,dimension(2)                             :: N     ! number of cells
    integer,dimension(0:N(1)+1,0:N(2)+1)             :: HOC   ! Head-of-Chain pointer
    real                                             :: tl    ! time

    integer               :: ic,icsh
    integer, dimension(2) :: i


    if (N(2).gt.2) then
       if (Ueps.eq.0) then
          HOC(:,0)        = HOC(:,N(2))
          HOC(:,N(2)+1)   = HOC(:,1)
       else
          do ic = 0,N(1)+1
             icsh = MOD(ic-1-INT(MOD(Ueps*(tl-T_sh_on),Lb(2))*REAL(N(1))/Lb(1))+2*N(1),N(1)) + 1
             HOC(ic,N(2)+1) = HOC(icsh,1)
             icsh = MOD(ic-1+INT(Ueps*(tl-T_sh_on)*REAL(N(1))/Lb(1)),N(1)) + 1
             HOC(ic,0) = HOC(icsh,N(2))
          end do
       end if
    end if

    if (N(1).gt.2) then
       HOC(0,:)        = HOC(N(1),:)
       HOC(N(1)+1,:)   = HOC(1,:)
    end if

  END SUBROUTINE makeper

END MODULE nlistmod




