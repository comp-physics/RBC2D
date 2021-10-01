MODULE stokes

  USE prms
  USE data
  USE membrane
  USE surfor
  USE expint
  USE nlistmod
  USE pme
  USE qflow
   
  IMPLICIT none

  real            :: acoef
  real            :: bcoef
  real            :: ccoef
  real            :: Co

  integer, dimension(2)                      :: Nc
  integer, allocatable, dimension(:)         :: lstS
  integer, allocatable,dimension(:,:)        :: nlstS

CONTAINS

  SUBROUTINE initrhs

    integer :: k1,k2

!    if (lambda.ne.1) stop
    acoef = -1./(2.*Pi*mu*(lambda+1.))
    bcoef = -4.*(1.-lambda)/( 2.*Pi*(lambda+1.) )
    ccoef = 2./(1. + lambda)
    Co = (1.-lambda)/(1.+lambda)/(2.*Pi)

    allocate(lstS(Nm*Nb*Nb*2),nlstS(Np,2))

    Nc = MAX(1,INT(Lb/rc/lfac))
    print *,"Nc",Nc

  END SUBROUTINE initrhs

  SUBROUTINE rhs(X,UU,t,lt)
    real, dimension(2,Np)  :: X
    real, dimension(2,Np)  :: UU
    real                      :: t
    integer                   :: lt   ! time step

    real, dimension(2,Np) :: Df    ! membrane force

    real, dimension(Np)   :: D     ! local stretch ds/dso
    real, dimension(2,Np) :: nvec,Dx2,Dx

    real, dimension(2,Np) :: Sr,Sf

    real, dimension(2,Np) :: UTMP,Umean
    real, dimension(2,Nm) :: Um

    integer :: Nmat
    integer :: i,m,n

    call construct_jumps(X,Df,D,Dx,Dx2,nvec)

    call wall(X,Df,D,Dx)

    call surforce(X,Df,D,lt,t)

    Sr = 0.; Sf = 0.  

    call computeSrp(X,alpha,Df,D,Dx,Sr,t)

    Df(1,:) = Df(1,:)*D*dSo
    Df(2,:) = Df(2,:)*D*dSo
    call computePMEV(X,Df,Sf,t)

    UU = acoef*(Sr + Sf)

    !call meanU(X,Umean)

    !UU = UU + Umean

    PRINT *,MAXVAL(UU)
    if (lambda.ne.1) then
       UTMP = UU 

       Nmat = Np*2/Nbfac
       call gmressolve(Nmat,10,X,D,nvec,Dx,Dx2,UU,t)

       call computeTrp(X,D,UU,nvec,Dx,Dx2,Sr,t)
       call computePMET(X,D,UU,nvec,Sf,t)
       call computeTo(X,D,UU,nvec,Dx)

       print *,"SSSS",SUM((UU + Co*(Sr+Sf+Dx) - UTMP)**2)/SUM(UU*UU)
    end if
    print *,"UMIN/MAX",MAXVAL(UU),MINVAL(UU)

!!$    do m = 1,Nm
!!$       n = 1 + (m-1)*Nb
!!$       Um(1,m) = SUM(UU(1,n:n+Nb-1))/REAL(Nb)
!!$       Um(2,m) = SUM(UU(2,n:n+Nb-1))/REAL(Nb)
!!$       do i = 0,Nb-1
!!$          UU(1,n+i) = UU(1,n+i) - Um(1,m)
!!$          UU(2,n+i) = UU(2,n+i) - Um(2,m)
!!$       end do
!!$    end do
!!$    Df(1,:) = nvec(1,:)*(UU(1,:)*nvec(1,:)+UU(2,:)*nvec(2,:))
!!$    Df(2,:) = nvec(2,:)*(UU(1,:)*nvec(1,:)+UU(2,:)*nvec(2,:))
!!$    UU = Df
!!$
!!$    do m = 1,Nm
!!$       n = 1 + (m-1)*Nb
!!$       do i = 0,Nb-1
!!$          UU(1,n+i) = UU(1,n+i) + Um(1,m)
!!$          UU(2,n+i) = UU(2,n+i) + Um(2,m)
!!$       end do
!!$    end do

  END SUBROUTINE rhs

  SUBROUTINE meanU(X,UU,time)
    real, dimension(2,Np) :: X
    real, dimension(2,Np) :: UU
real :: time
real :: Umaxt

!  If (time < Real(Nt)*dt/15. ) Then
!  If (time < Real(Nt)*dt ) Then
!  Umaxt = 0
!  Else
!  Umaxt = Umax
!  End If

       UU = 0.
  UU(1,:) = Umax
  UU(2,:) = 0

  if (which_ic .eq. 0) Then
    UU(1,:) = Umax
    UU(2,:) = 0
  end if
    ! UU(1,:) = UU(1,:) + ccoef*Umax !*COS(4.*Pi*X(2,:)/Lb)
  !  UU(1,:) = UU(1,:) + ccoef*Umax*(1.+SIN(2.*Pi*X(2,:)/Lb))
 !   UU(2,:) = UU(2,:) ! + ccoef*Umax*SIN(4*2.*Pi*X(1,:)/Lb(2))
!    UU(1,:) = UU(1,:) + ccoef*Ueps*(X(2,:)-Lb(2)/2)/Lb(2)

!    call getQflow(Np,X,UU)
!    UU = UU *ccoef

  END SUBROUTINE meanU

  SUBROUTINE computeSrp(X,alf,Df,D,Dx,Sr,tl)
    real, dimension(2,Np)     :: X
    real, dimension(2,Np)     :: Df,Dx  
    real, dimension(Np)       :: D
    real, dimension(2,Np)     :: Sr
    real                         :: alf
    real                         :: tl

    real, dimension(2)           :: xnm
    real                         :: r2
    integer                      :: ln,lm,i,j,n,m,Is,jj,ll,lm1,mset
    real, parameter              :: Eo = 1.7724538509
    real                         :: DS

    Sr = 0.

    do m = 1,Np
       do ll = nlstS(m,1),nlstS(m,2)
          n = lstS(ll)

          xnm = X(:,n) - X(:,m) 
          xnm(1) = xnm(1) - NINT(xnm(2)/Lb(2))*Ueps*(tl-T_sh_on)
          xnm = xnm - NINT(xnm/Lb)*Lb

          r2 = SUM(xnm*xnm) 
          if (n.eq.m) then
             Sr(:,m) = Sr(:,m) + (Dx(:,m)*SUM(Df(:,m)*Dx(:,m)) &
                  - SUM(Dx(:,m)**2)*Df(:,m))/D(m)*dSo(m) &
                  + Df(:,m)*SQRT(alf/Pi)*Eo

          else 
             Sr(:,m) = Sr(:,m) &
                  + Ei(Pi*r2/alf)*Df(:,n)/2.*D(n)*dSo(n) &
                  + EXP(-Pi*r2/alf)*(xnm(:)*SUM(xnm(:)*Df(:,n))/r2 &
                  - Df(:,n))*D(n)*dSo(n) 
          end if

       end do
    end do

    do m = 1,Np
       Is = INT(rc/dSo(m))
       do jj = 1,Is
          DS = REAL(jj)*dSo(m)
          Sr(:,m) = Sr(:,m) &
               - Ei(Pi*DS**2*D(m)**2/alf)*Df(:,m)*D(m)*dSo(m)
       end do
    end do
!    Sr = Sr*Soe(2)

    call filtknots(2,Sr)

  END SUBROUTINE computeSrp


  SUBROUTINE computeTo(X,D,UU,nvec,Bo)
    real, dimension(2,Np)     :: X
    real, dimension(2,Np)     :: nvec
    real, dimension(2,Np)     :: UU
    real, dimension(Np)       :: D
    real, dimension(2,Np)     :: Bo

    integer :: n,k,i

    real, dimension(2)   :: Si
    real, dimension(2,2) :: Sik

!    real, dimension(Nf,2,Nb)  :: Xb

    !    Xb = X - FLOOR(X/Lb)*Lb 

    Si = 0.; Sik = 0.
    do n = 1,Np
       Si = Si + SUM(X(:,n)*nvec(:,n))*UU(:,n)*D(n)
       do k = 1,2
          Sik(:,k) = Sik(:,k) + nvec(k,n)*UU(:,n)*D(n)
       end do
    end do

    do n = 1,Np
       do i = 1,2
          Bo(i,n) = 4.*Pi/tau0*Si(i)! - SUM(X(lm,:,m)*Sik(i,:)))
       end do
    end do
    Bo = 0.

  END SUBROUTINE computeTo

  SUBROUTINE computeTrp(X,D,UU,nvec,Dx,Dx2,Br,tl)
    real, dimension(2,Np)     :: X
    real, dimension(2,Np)     :: nvec
    real, dimension(2,Np)     :: UU
    real, dimension(Np)       :: D
    real, dimension(2,Np)     :: Dx    ! dx/dso
    real, dimension(2,Np)     :: Dx2
    real                         :: tl     ! local time

    real, dimension(2,Np)     :: Br

    real, dimension(2)           :: xnm
    real                         :: r2
    integer                      :: ln,lm,i,n,m,ll,lm1
    real                         :: DS

    real, dimension(Np,2)        :: dUU,dnvec
    real, dimension(Np)          :: DD

    real, dimension(Nb)          :: utmp


    real                         :: rfac

    Br = 0.
    do n = 1,Nm

       m = 1 + (n-1)*Nb

       utmp = UU(1,m:m+Nb-1)
       call dersp(utmp,len0(n),dUU(m,1))
       utmp = UU(2,m:m+Nb-1)
       call dersp(utmp,len0(n),dUU(m,2))

       utmp = nvec(1,m:m+Nb-1)
       call dersp(utmp,len0(n),dnvec(m,1))
       utmp = nvec(2,m:m+Nb-1)
       call dersp(utmp,len0(n),dnvec(m,2))

       call dersp(D(m),len0(n),DD(m))

    end do

    do m = 1,Np
       do ll = nlstS(m,1),nlstS(m,2)
          n = lstS(ll)


          xnm = X(:,n) - X(:,m) 
          xnm(1) = xnm(1) - NINT(xnm(2)/Lb(2))*Ueps*(tl-T_sh_on)
          xnm = xnm - NINT(xnm/Lb)*Lb

          r2 = SUM(xnm*xnm) 
          if (n.ne.m) then
             rfac = 4.*EXP(-Pi*r2/alpha)*(1.+Pi*r2/alpha) &
                  *SUM(xnm*UU(:,n))*SUM(xnm*nvec(:,n)) &
                  /r2**2*D(n)*dSo(n)
             do i = 1,2
                Br(i,m) = Br(i,m) + xnm(i)*rfac !&
             end do

          else if (n.eq.m) then
             Br(:,m) = Br(:,m) &
                  +2.*(2.*Dx(:,n)*SUM(Dx(:,n)*UU(:,n)) &
                  *SUM(Dx(:,n)*dnvec(n,:)) &
                  + 2.*Dx(:,n)*SUM(Dx(:,n)*dUU(n,:)) &
                  *SUM(Dx(:,n)*nvec(:,n)) &
                  + 2.*Dx(:,n)*SUM(Dx(:,n)*UU(:,n)) &
                  *SUM(Dx(:,n)*nvec(:,n))*DD(n)/D(n) &
                  + Dx(:,n)*SUM(Dx(:,n)*UU(:,n)) &
                  *SUM(Dx2(:,n)*nvec(:,n)) &
                  + Dx(:,n)*SUM(Dx2(:,n)*UU(:,n)) &
                  *SUM(Dx(:,n)*nvec(:,n)) &
                  + Dx2(:,n)*SUM(Dx(:,n)*UU(:,n)) &
                  *SUM(Dx(:,n)*nvec(:,n))) &
                  /D(n)**3*dSo(m)

          end if

       end do

    end do

    call filtknots(2,Br)


  END SUBROUTINE computeTrp


  SUBROUTINE mknlist(Xin,tl)
    real   ,dimension(2,Np)                              :: Xin   ! particle positions
    real                                                    :: tl    ! local time

    real   ,dimension(2,Np)                              :: X  
    real   ,dimension(2,Np)                              :: Xw 
    integer,dimension(Np)                                :: LL  ! next in list pointer
    integer,dimension(0:Nc(1)+1,0:Nc(2)+1)               :: HOC ! Head-of-Chain pointer
    ! with perioidic continuation

    real,dimension(2)     :: xx
    integer               :: ic1,ic2
    real                  :: rij2,rc2l
    integer               :: i,j
    integer               :: k,l,m,n,n1,n2,n3
    integer               :: nbr_cnt

    rc2l = rc*rc*lfac*lfac           ! list cutoff for PME

    Xw(1,:) = Xin(1,:) - FLOOR( Xin(2,:)/Lb(2) )*Ueps*(tl-T_sh_on)
    Xw(2,:) = Xin(2,:)
    call chainlist(LL,HOC,Nc,Xw,tl)

    nbr_cnt = 1
    nlstS(:,2) = -1
    nlstS(:,1) = 0

    do ic2 = 1,Nc(2)
       do ic1 = 1,Nc(1) 
          i = HOC(ic1,ic2)                ! start with atom at the head of the linked list
          do while(i.ne.0)                ! loop over all atoms in the linked list if any
             nlstS(i,1) = nbr_cnt

             ! loop over current and all nieghboring cells
             do n1 = -1,1
                do n2 = -1,1
                   j = HOC(ic1+n1,ic2+n2)
                   do while(j.ne.0)

                      xx = Xin(:,i) - Xin(:,j) 
                      xx(1) = xx(1) - NINT(xx(2)/Lb(2))*Ueps*(tl-T_sh_on)
                      xx = xx - NINT(xx/Lb)*Lb

                      rij2 = SUM(xx*xx)

                      lstS(nbr_cnt) = j
                      nbr_cnt = nbr_cnt + 1

                      j = LL(j)

                   end do

                end do
             end do

             nlstS(i,2) = nbr_cnt-1
             i = LL(i)

          end do

       end do
    end do

    print *,"stokes:",nbr_cnt
!!$    print *,Nm*Nb*Nb*20
  END SUBROUTINE mknlist


  SUBROUTINE gmressolve(n,m,X ,D ,nvec ,Dx ,Dx2, B,tl)
    integer                     :: n,m
    real, dimension(Nbfac*n)          :: X
    real, dimension(Nbfac*n)          :: nvec
    real, dimension(Nbfac*n)          :: Dx    ! dx/dso
    real, dimension(Nbfac*n)          :: Dx2
    real, dimension(Nbfac*n/2)        :: D

    real, dimension(Nbfac*n)          :: B
    real                        :: tl

    real, dimension(Nbfac*n)          :: Br,Bf,Bo
    real, dimension(Nbfac*n)          :: UTMP

    integer            :: lwork

    integer            :: i, j
    integer            :: revcom, colx, coly, colz, nbscal

    integer, dimension(5) :: irc
    integer, dimension(8) :: icntl
    integer, dimension(3) :: info

    integer, parameter :: matvec = 1
    integer, parameter :: precondLeft = 2
    integer, parameter :: precondRight = 3
    integer, parameter :: dotProd = 4

    integer            :: nout
    
    real, dimension(m**2 + m*(n+5) + 6*n + 1) :: work
    real, dimension(5)     :: cntl
    real, dimension(2)     :: rinfo

    real, parameter :: ZERO = 0.
    real, parameter :: ONE = 1.

    lwork = m**2 + m*(n+5) + 6*n + 1

!     print *,lwork; stop

!    call DGEMV('N',n,n,ONE,A,n,work(1),1,ZERO,work(n+1),1)
     
 !    print *,Xlast(1)
!!$    work(1:n) = Xlast
!!$    i = 4*n + n*m + (m+1)**2 + m 
!!$    work(i+1:i+1+n) = Xlast
    work(n+1:2*n:2) = B(1:n*Nbfac:2*Nbfac)
    work(n+1+1:2*n:2) = B(2:n*Nbfac:2*Nbfac)


    call init_dgmres(icntl,cntl)

    cntl = 0. ; icntl = 0
    ! Tolerance
    cntl(1) = 1.E-7
 !   icntl(2) = 0
 !   icntl(6) = 0
    ! Save the convergence history in file fort.20
    icntl(3) = 0
    ! No preconditioning
    icntl(4) = 1
    ! ICGS orthogonalization
    icntl(5) = 3
    ! Maximum number of iterations
    icntl(7) = 100 

    !*****************************************
    !** Reverse communication implementation
    !*****************************************
    
10  call drive_dgmres(n,n,m,lwork,work,irc,icntl,cntl,info,rinfo)
    revcom = irc(1)
    colx   = irc(2)
    coly   = irc(3)
    colz   = irc(4)
    nbscal = irc(5)
    
    if (revcom.eq.matvec) then
       ! perform the matrix vector product work(colz) <-- A * work(colx)
       call Nf2Nb(work(colx),UTMP)
       call computeTrp(X,D,UTMP,nvec,Dx,Dx2,Br,tl)
       call computePMET(X,D,UTMP,nvec,Bf,tl)

       work(colz:colz+n-1:2) = work(colx:colx+n-1:2) + Co*(Br(1:Nbfac*n:2*Nbfac)+Bf(1:Nbfac*n:2*Nbfac)+Bo(1:Nbfac*n:2*Nbfac))
       work(colz+1:colz+n:2) = work(colx+1:colx+n:2) + Co*(Br(2:Nbfac*n:2*Nbfac)+Bf(2:Nbfac*n:2*Nbfac)+Bo(2:Nbfac*n:2*Nbfac))

!       call sgemv('N',n,n,ONE,a,n,work(colx),1,ZERO,work(colz),1)
       goto 10
    else if (revcom.eq.precondLeft) then
       ! perform the left preconditioning work(colz) <-- M^{-1} * work(colx)
       call dcopy(n,work(colx),1,work(colz),1)
       goto 10
    else if (revcom.eq.precondRight) then
       ! perform the right preconditioning
       call dcopy(n,work(colx),1,work(colz),1)
       goto 10
    else if (revcom.eq.dotProd) then
       ! perform the scalar product work(colz) <-- work(colx) work(coly)
       call dgemv('C',n,nbscal,ONE,work(colx),n,work(coly),1,ZERO,work(colz),1)
       goto 10
    endif

    call Nf2Nb(work(1),B)

  END SUBROUTINE gmressolve

END MODULE stokes
