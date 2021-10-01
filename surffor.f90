MODULE surfor 

  USE prms
  USE data
  USE nlistmod
  USE membrane

  IMPLICIT none

  integer, dimension(2)                      :: NcS
  integer, allocatable, dimension(:)         :: lst
  integer, allocatable,dimension(:,:)        :: nlst

CONTAINS

  SUBROUTINE initsurfor

    allocate(lst(Nb*Np*2),nlst(Np,2))

    NcS = MAX(1,INT(Lb/dnl/lfacdnl))
    print *,"Nc-surfor",NcS

  END SUBROUTINE initsurfor
 

  SUBROUTINE nlister(Xin,tl)
    real   ,dimension(2,Np)                              :: Xin   ! particle positions
    real                                                    :: tl    ! local time

    real   ,dimension(2,Np)                              :: X
    real   ,dimension(2,Np)                              :: Xw
    integer,dimension(Np)                                :: LL  ! next in list pointer
    integer,dimension(0:NcS(1)+1,0:NcS(2)+1)             :: HOC ! Head-of-Chain pointer
    ! with perioidic continuation

    real,dimension(2)     :: xx
    integer               :: ic1,ic2,i0,ii,jj
    real                  :: rij2,rc2l
    integer               :: i,j
    integer               :: k,l,m,n,n1,n2

    integer               :: nbr_cnt


    rc2l = (dnl*lfacdnl)**2

    i0 = INT(lfacdnl*2.*dnl/MINVAL(dSo))+8


    Xw(1,:) = Xin(1,:) - FLOOR( Xin(2,:)/Lb(2) )*Ueps*(tl-T_sh_on)
    Xw(2,:) = Xin(2,:)

    call chainlist(LL,HOC,NcS,Xw,tl)

!!$    do m = 1,Nm
!!$       do l = 1,Nb
!!$          write(77,*)LL(l,:,m)
!!$       end do
!!$    end do
!!$    stop

    nbr_cnt = 1
    nlst(:,2) = -1
    nlst(:,1) = 0

    do ic2 = 1,NcS(2)
       do ic1 = 1,NcS(1)
          i = HOC(ic1,ic2)           ! start with atom at the head of the linked list
          whilei: do while(i.ne.0)          ! loop over all atoms in the linked list if any

             nlst(i,1) = nbr_cnt

             do n1 = -1,1
                do n2 = -1,1
                   j = HOC(ic1+n1,ic2+n2)
                   whilej: do while(j.ne.0)

                      xx = Xin(:,i) - Xin(:,j) 
                      xx(1) = xx(1) - NINT(xx(2)/Lb(2))*Ueps*(tl-T_sh_on)
                      xx = xx - NINT(xx/Lb)*Lb


                      rij2 = SUM(xx*xx)

                      if (rij2.lt.rc2l.and. (cell(i).ne.cell(j).or.&
                         (cell(i).eq.cell(j).and.cell(i).gt.0.and.ABS(j-i).gt.i0.and.ABS(i-j-1).lt.(Nb-i0))) ) then

                         lst(nbr_cnt) = j
                         nbr_cnt = nbr_cnt + 1

                      end if

                      j = LL(j)
!                      print *,j(1)

                   end do whilej

                end do
             end do

             nlst(i,2) = nbr_cnt-1
             i = LL(i)

          end do whilei

       end do
    end do
    print *,"surfor:",nbr_cnt
!!$    print *, Nm*Nb*Nb*20

  END SUBROUTINE nlister

  SUBROUTINE surforce(X,ff,D,lt,tl)
    real, dimension(2,Np)     :: X
    real, dimension(2,Np)     :: ff
    real, dimension(Np)       :: D
    real                      :: tl       ! local time
    integer                   :: lt       ! time step index

    integer :: n,m,i,j,l,ia,iia,ll
    real    :: fac,d2,s2,r2,r4,r6,nxi,nxj,fm,rr,r3,r5,f1
    real, dimension(2)  :: xij,fo,xij0
    integer, parameter :: na = 2

    real, dimension(2,Np)     :: ft,Xeo
    character(20)  :: fn

    integer :: nn = 1


    fm = 0.;ft = 0.
    d2 = dnl*dnl

    do n = 1,Np
       do ll= nlst(n,1),nlst(n,2)
          m = lst(ll)

          xij = X(:,n) - X(:,m) 
          xij(1) = xij(1) - NINT(xij(2)/Lb(2))*Ueps*(tl-T_sh_on)
          xij = xij - NINT(xij/Lb)*Lb

          r2 = SUM( xij**2 )
          if (r2.lt.d2) then
             rr = SQRT(r2)
             !      fo = Sf*(rr-dnl)*xij/rr
             if (rr.lt.dnl-sigma) then
                f1 = -(EXP(dnl-sigma-rr)-1.)/(EXP(dnl-sigma)-1.) 
             else if (cell(n).ne.0 .and. cell(m).ne.0.and.cell(n).ne.cell(m)) then
                f1 = Sa*( 1. - ((rr-dnl+sigma/2.)/(sigma/2.))**2 )
                !f1 = f1 - 2.*Sa*sigma*(rr-dnl)*EXP(-sigma*(rr-dnl)**2)
             else
                f1 = 0.
             end if
!             write(89,*)rr,f1
             fo = Sf*f1*xij/rr
             fm = MAX(fm,SQRT(SUM(fo*fo)))
             ft(:,n) = ft(:,n) + fo

!!$          else if (cell(n).eq.cell(m).and.r2.lt.4.*d2) then
!!$
!!$             rr = SQRT(r2)
!!$             !      fo = Sf*(rr-dnl)*xij/rr
!!$             fo = -Sf*(EXP(2.*dnl-rr)-1.)/(EXP(2.*dnl)-1.)*xij/rr
!!$             fm = MAX(fm,SQRT(SUM(fo*fo)))
!!$             ft(:,n) = ft(:,n) + fo-

          end if

       end do
    end do

    call filtknots(2,ft)

    ff(1,:) = ff(1,:) + ft(1,:)/D
    ff(2,:) = ff(2,:) + ft(2,:)/D 

    call filtknots(2,ff)

    !    print *,"fm",i_whls,fm,MAXVAL(ft)

    !    if (f_out.gt.0.and.MOD(lt,f_out).eq.0) then
!!$       Xeo(:,1,:) = Xe(:,1,:) - FLOOR( Xe(:,2,:)/Lb )*MOD(Ueps*tl,Lb)
!!$       Xeo(:,2,:) = Xe(:,2,:)
!!$       write(fn,"('D/nn.fout.',I10.10)")nn ! lt
!!$       nn = nn + 1
!!$       open(1,file=fn)
!!$       do n = 1,Nm
!!$          do i = 1,Nb
!!$             write(1,"(2F20.10)")Xeo(i,:,n)- FLOOR(Xeo(i,:,n)/Lb)*Lb
!!$             write(1,"(2F20.10)")Xeo(i,:,n)- FLOOR(Xeo(i,:,n)/Lb)*Lb-ft(i,:,n)/MAXVAL(ft)*0.1
!!$             write(1,"(2F20.10)")Xeo(i,:,n)- FLOOR(Xeo(i,:,n)/Lb)*Lb
!!$          end do
!!$       end do
!!$       close(1)
    !    end if

  END SUBROUTINE surforce

END MODULE surfor
