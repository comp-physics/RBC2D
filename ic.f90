MODULE icmod

  USE prms
  USE data

  IMPLICIT none

  integer  :: lt0  ! initial time step
  real     :: time0   ! initial time

CONTAINS

  SUBROUTINE ic

    real    :: th,r,dSow
    integer :: i,m,n,i1,i2,ncnt=0
    integer, dimension(1) :: luk
    integer :: Np_l,Nb_l,Nm_l

    len0 = l0*2.*Pi*ro
    rad = ro
    Tc = T0
    Mc = M0
    !rad(1) = r1
    !Tc(1) = T1
    !Mc(1) = M1
    do i = 1,Nb
       Soe(i) = l0*2.*Pi*ro*REAL(i-1)/REAL(Nb)
    end do
    do m = 1,Nm
       do i = 1,Nb
          dSo((m-1)*Nb+i) = len0(m)/REAL(Nb)
       end do
    end do

    lt0 = 0      ! initial time step :: zero for most cases
    time0 = 0.   ! initial time      :: zero for most cases

    cell = 0
    do m = 1,Nm
       cell(1+(m-1)*Nb:m*Nb) = m
    end do

    select case (which_ic)

     ! case(113)
! print*, 'pert ****', pert
    ! Tc = T0; Mc = M0
     !  ! if (Nm.ne.18) stop 'bad Nm'
       ! do m = 1,Nm
       !    do n = 1,1
       !       call RANDOM_NUMBER(r)
       !       do i = 1,Nb
       !          th = 2.*Pi*REAL(i-1)/REAL(Nb)
       !          X(1,i+Nb*(n - 1 + (m-1)*3)) = &
       !              Lb/5.4*REAL(m-.5) + rad(n + (m-1)*5)*COS(th) + 0.2*(r-0.5)
       !          X(2,i+Nb*(n - 1+ (m-1)*1)) = &
       !              2*rad(2)*SIN(th) - 0.8*(r-0.5) + Real(Lb(2))/2.
       !       end do
       !    end do
       ! end do

!give one cell a large offset
!        do i = 1,Nb
!        th = 2.*Pi*REAL(i-1)/REAL(Nb)
!        X(1,i+Nb*((3-1)*1)) = Lb/10*REAL(3-.5) + rad(2)*COS(th)
!        X(2,i+Nb*((3-1)*1)) = rad(2)*SIN(th) - 2 + 20/2
!        end do

        !dSow = Lb/REAL(Nw)
      !dSow = Lb(1)/(REAL(Nw))
      !!bottom wall
      ! do i = 1,Nw
      !    wX(1,i) = REAL(i-1)*dSow
      !    wX(2,i) = 0.0
      ! end do

      ! wDx(1,:) = dSow/Soe(2)
      ! wDx(2,:) = 0.000
      ! wD(:) = dSow/Soe(2)
      ! X(:,Ns+1:Np) = wX

            case(114)
    Tc = T0; Mc = M0
       if (Nm.ne.25) stop 'bad Nm'
       !do m = 1,18
          !do n = 1,1
             !call RANDOM_NUMBER(r)
             !do i = 1,Nb
                !th = 2.*Pi*REAL(i-1)/REAL(Nb)
                  !X(1,i+Nb*(n - 1 + (m-1)*1)) = Lb(1)/18*REAL(m-.5) + 1*rad(2)*COS(th)
                  !X(2,i+Nb*(n - 1+ (m-1)*1)) = rad(2)*SIN(th) - 0.8*(r-0.5) + Real(Lb(2))/8.
             !end do
          !end do
       !end do

       do m = 1,5
          do n = 1,5
             call RANDOM_NUMBER(r)
             r=0
             do i = 1,Nb
                th = 2.*Pi*REAL(i-1)/REAL(Nb)
                X(1,i+Nb*(n - 1 + (m-1)*5)) = Lb(1)/5.4*REAL(m-.5) + (rad(n + (m-1)*5)+0.1*COS(2.*th+2.*Pi*0))*COS(th) + 0.2*(r-0.5)
                X(2,i+Nb*(n - 1+ (m-1)*5)) = Lb(2)/5*REAL(n-.5) + (rad(n + (m-1)*5)+0.1*COS(2.*th+2.*Pi*0))*SIN(th) - 0.2*(r-0.5)
             end do
          end do
       end do


      dSow = Lb(1)/(REAL(Nw)/2.)
      !bottom wall
       do i = 1,Nw/2
          wX(1,i) = REAL(i-1)*dSow
          wX(2,i) = 0.0
       end do
       !top wall
       do i = Nw/2+1,Nw
          wX(1,i) = REAL(i-1-Nw/2)*dSow
          wX(2,i) = Real(Lb(2))/1.01  ! vert position of wall
       end do

       wDx(1,:) = dSow/Soe(2)
       wDx(2,:) = 0.000
       wD(:) = dSow/Soe(2)
       X(:,Ns+1:Np) = wX


    case(101)
       if (Nm.ne.0) stop 'bad Nm'

    case(102)

       if (Nm.ne.1) stop 'bad Nm'

       do i = 1,Nb
          th = 2.*Pi*REAL(i-1)/REAL(Nb)
          X(1,i) = Lb(1)/2 + ro*COS(th)
          X(2,i) = Lb(2)/4 + ro*SIN(th)
       end do

       dSow = Lb(1)/(REAL(Nw)/2.)
       do i = 1,Nw/2
          wX(1,i) = REAL(i-1)*dSow
          wX(2,i) = 0.001
       end do

       do i = Nw/2+1,Nw
          wX(1,i) = REAL(i-1-Nw/2)*dSow
          wX(2,i) = 0.5*Lb(2)
       end do

       wDx(1,:) = dSow/Soe(2)
       wDx(2,:) = 0.
       wD(:) = dSow/Soe(2)
       X(:,1:Np) = wX

    case(1)
       if (Nm.ne.1) stop 'bad Nm'
       do i = 1,Nb
          th = 2.*Pi*REAL(i-1)/REAL(Nb)
          X(1,i) = Lb(1)/2 + r1*COS(th) - 1.2*ro
          X(2,i) = Lb(2)/2 + r1*SIN(th)
       end do
    case(2)
       if (Nm.ne.2) stop 'bad Nm'
       do i = 1,Nb
          th = 2.*Pi*REAL(i-1)/REAL(Nb)
          X(1,i) = Lb(1)/4. +COS(th)
          X(2,i) = Lb(2)/4 + SIN(th)
          X(1,i+Nb) = 3.*Lb(1)/4. + COS(th)
          X(2,i+Nb) = 3.*Lb(2)/4. + SIN(th)
       end do
    case(-2)
       do i = 1,Nb
          Soe(i) = REAL(i-1)*len0(m)/REAL(Nb)
       end do
       do m = 1,Nm
          rad(m) = ro
          len0(m) = 2.*Pi*rad(m)
          dSo((m-1)*Nb+1:m*Nb) = len0(m)/REAL(Nb)
       end do
       if (Nm.ne.2) stop 'bad Nm'
       do i = 1,Nb
          th = 2.*Pi*REAL(i-1)/REAL(Nb)
          X(1,i) = Lb(1)/4. +ro*COS(th)
          X(2,i) = Lb(2)/2 +ro*SIN(th)+0.25*ro
          X(1,i+Nb) = 3.*Lb(1)/4. +ro*COS(th)
          X(2,i+Nb) = Lb(2)/2 +ro*SIN(th)-0.25*ro
       end do
    case(-3)
       do i = 1,Nb
          Soe(i) = REAL(i-1)*len0(m)/REAL(Nb)
       end do
       do m = 1,Nm
          call RANDOM_NUMBER(r)
          rad(m) = ro + (r1-ro)*r
          len0(m) = 2.*Pi*rad(m)
          dSo((m-1)*Nb+1:m*Nb) = len0(m)/REAL(Nb)
       end do
       i1 = INT(Lb(1)/(2.*ro))
       i2 = INT(Lb(2)/(2.*ro))
       print *,i1,i2
       do m = 1,i1
          do n = 1,i2
             if (ncnt.le.Nm) then
                ncnt = ncnt + 1
                call RANDOM_NUMBER(r)
                do i = 1,Nb
                   th = 2.*Pi*REAL(i-1)/REAL(Nb)
                   X(1,i+Nb*(ncnt-1)) = Lb(1)/REAL(i1)*REAL(m-.5) + (rad(ncnt))*COS(th) + 0.2*ro*(r-0.5)
                   X(2,i+Nb*(ncnt-1)) = Lb(2)/REAL(i2)*REAL(n-.5) + (rad(ncnt))*SIN(th) - 0.2*ro*(r-0.5) + Lb(2)/70.*REAL(m-1)
                end do
             end if
          end do
       end do
!       write(66,"(2E20.10)")X;stop
    case(-25)
       do i = 1,Nb
          Soe(i) = REAL(i-1)*len0(m)/REAL(Nb)
       end do
       do m = 1,Nm
          call RANDOM_NUMBER(r)
          rad(m) = ro + (r1-ro)*r
          len0(m) = 2.*Pi*rad(m)
          dSo((m-1)*Nb+1:m*Nb) = len0(m)/REAL(Nb)
       end do
       if (Nm.ne.25) stop 'bad Nm'
       do m = 1,5
          do n = 1,5
             call RANDOM_NUMBER(r)
             do i = 1,Nb
                th = 2.*Pi*REAL(i-1)/REAL(Nb)
                X(1,i+Nb*(n - 1 + (m-1)*5)) = Lb(1)/5.4*REAL(m-.5) + (rad(n + (m-1)*5))*COS(th) + 0.2*ro*(r-0.5)
                X(2,i+Nb*(n - 1+ (m-1)*5)) = Lb(2)/5.4*REAL(n-.5) + (rad(n + (m-1)*5))*SIN(th)&
                  & - 0.2*ro*(r-0.5) + Lb(2)/70.*REAL(m-1)
             end do
          end do
       end do
    case(25)
       !rad(Nm) = r1
       !Tc(Nm) = T1
       !Mc(Nm) = M1
    Tc = T0; Mc = M0
       if (Nm.ne.25) stop 'bad Nm'
       do m = 1,5
          do n = 1,5
             call RANDOM_NUMBER(r)
             do i = 1,Nb
                th = 2.*Pi*REAL(i-1)/REAL(Nb)
                X(1,i+Nb*(n - 1 + (m-1)*5)) = Lb(1)/5.4*REAL(m-.5) + COS(th) + 4.0*(r-0.5)
                X(2,i+Nb*(n - 1+ (m-1)*5)) = Lb(2)/5.4*REAL(n-.5) + SIN(th) - 0.2*(r-0.5)
             end do
          end do
       end do

      dSow = Lb(1)/REAL(Nw/2.)
      !bottom wall
       do i = 1,Nw/2
          wX(1,i) = REAL(i-1)*dSow
          wX(2,i) = 0.01
       end do

       !top wall
       do i = Nw/2+1,Nw
          wX(1,i) = REAL(i-1-Nw/2)*dSow
          wX(2,i) = Real(Lb(2))/1.03  ! vert position of wall
       end do

       wDx(1,:) = dSow/Soe(2)
       wDx(2,:) = 0.000
       wD(:) = dSow/Soe(2)
       X(:,Ns+1:Np) = wX

    case(100)
       rad = ro; Tc = T0; Mc = M0
       rad(Nm) = r1
       Tc(Nm) = T1
       Mc(Nm) = M1

       if (Nm.ne.100) stop 'bad Nm'
       do m = 1,10
          do n = 1,10
             call RANDOM_NUMBER(r)
             do i = 1,Nb
                th = 2.*Pi*REAL(i-1)/REAL(Nb)
                X(1,i+Nb*(n - 1 + (m-1)*10)) = Lb(1)/10.7*REAL(m-.5) + rad(n + (m-1)*10)*COS(th) - 0.1*(r-0.5)
                X(2,i+Nb*(n - 1 + (m-1)*10)) = Lb(2)/10.7*REAL(n-.5) + rad(n + (m-1)*10)*SIN(th) + 0.1*(r-0.5) + Lb(2)/300.*REAL(m)
             end do
          end do
       end do
       X(2,98*Nb:99*Nb-1) = X(2,98*Nb:99*Nb-1) - 0.35*ro
       X(1,89*Nb:90*Nb-1) = X(1,89*Nb:90*Nb-1) - 0.3*ro
       dSow = Lb(1)/REAL(Nw)
       do i = 1,Nw
          wX(1,i) = REAL(i-1)*dSow
          wX(2,i) = 0.
       end do
       wDx(1,:) = dSow/Soe(2)
       wDx(2,:) = 0.
       wD(:) = dSow/Soe(2)
       X(:,Ns+1:Np) = wX

    case(-900)
       do i = 1,Nb
          Soe(i) = REAL(i-1)*len0(m)/REAL(Nb)
       end do
       do m = 1,Nm
          call RANDOM_NUMBER(r)
          rad(m) = ro + (r1-ro)*r
          len0(m) = 2.*Pi*rad(m)
          dSo((m-1)*Nb+1:m*Nb) = len0(m)/REAL(Nb)
       end do
       if (Nm.ne.900) stop 'bad Nm'
       do m = 1,30
          do n = 1,30
             call RANDOM_NUMBER(r)
             do i = 1,Nb
                th = 2.*Pi*REAL(i-1)/REAL(Nb)
                X(1,i+Nb*(n - 1 + (m-1)*30)) = Lb(1)/30.*REAL(m-.5) + (rad(n + (m-1)*30))*COS(th) + 0.2*(r-0.5)
                X(2,i+Nb*(n - 1+ (m-1)*30)) = Lb(2)/30.*REAL(n-.5) + (rad(n + (m-1)*30))*SIN(th) - 0.2*(r-0.5) + Lb(2)/70.*REAL(m-1)
             end do
          end do
       end do
    case(-100)
       do i = 1,Nb
          Soe(i) = REAL(i-1)*len0(m)/REAL(Nb)
       end do
       do m = 1,Nm
          call RANDOM_NUMBER(r)
          rad(m) = ro + (r1-ro)*r
          len0(m) = 2.*Pi*rad(m)
          dSo((m-1)*Nb+1:m*Nb) = len0(m)/REAL(Nb)
       end do
       if (Nm.ne.100) stop 'bad Nm'
       do m = 1,10
          do n = 1,10
             call RANDOM_NUMBER(r)
             do i = 1,Nb
                th = 2.*Pi*REAL(i-1)/REAL(Nb)
                X(1,i+Nb*(n - 1 + (m-1)*10)) = Lb(1)/10.*REAL(m-.5) + (rad(n + (m-1)*10))*COS(th) + 0.2*ro*(r-0.5)
                X(2,i+Nb*(n - 1+ (m-1)*10)) = Lb(2)/10.*REAL(n-.5) + (rad(n + (m-1)*10))*SIN(th) &
                  & - 0.2*ro*(r-0.5) + Lb(2)/70.*REAL(m-1)
             end do
          end do
       end do
       !       write(77,"(2F20.10)")X; stop
    case(31)
       rad = ro; Tc = T0; Mc = M0
       rad(Nm-2) = r1
       Tc(Nm-2) = T1
       Mc(Nm-2) = M1

       if (Nm.ne.30) stop 'bad Nm'
       do m = 1,10
          do n = 1,3
             call RANDOM_NUMBER(r)
             do i = 1,Nb
                th = 2.*Pi*REAL(i-1)/REAL(Nb)
                X(1,i+Nb*(n - 1 + (m-1)*3)) = Lb(1)/11.5*REAL(m-.5) + rad(n + (m-1)*3)*COS(th) + 0.01*(r-0.5)
                X(2,i+Nb*(n - 1 + (m-1)*3)) = Lb(2)/5.5*REAL(n-.5) + rad(n + (m-1)*3)*SIN(th) - 0.01*(r-0.5) + Lb(2)/2 - 2*ro
                end do
          end do
       end do
       X(2,28*Nb+1:29*Nb) = X(2,28*Nb+1:29*Nb) - 0.1*ro
       X(2,29*Nb+1:30*Nb) = X(2,29*Nb+1:30*Nb) - 0.1*ro
       X(1,27*Nb+1:28*Nb) = X(1,27*Nb+1:28*Nb) + 1.7*ro
       dSow = Lb(1)/(REAL(Nw)/2.)
       do i = 1,Nw/2
          wX(1,i) = REAL(i-1)*dSow
          wX(2,i) = 0.
       end do
       do i = Nw/2+1,Nw
          wX(1,i) = REAL(i-1-Nw/2)*dSow
          wX(2,i) = 0.99*Lb(2)
       end do
       wDx(1,:) = dSow/Soe(2)
       wDx(2,:) = 0.
       wD(:) = dSow/Soe(2)
       X(:,Ns+1:Np) = wX


    case(32)
       rad = ro; Tc = T0; Mc = M0
       !rad(Nm-2) = r1
       !Tc(Nm-2) = T1
       !Mc(Nm-2) = M1

       if (Nm.ne.10) stop 'bad Nm'
       do m = 1,10
             call RANDOM_NUMBER(r)
             do i = 1,Nb
                th = 2.*Pi*REAL(i-1)/REAL(Nb)
                X(1,i+Nb*(1 - 1 + (m-1)*3)) = Lb(1)/11.5*REAL(m-.5) + rad(1 + (m-1)*3)*COS(th) + 0.01*(r-0.5)
                X(2,i+Nb*(1 - 1 + (m-1)*3)) = Lb(2)/5.5*REAL(1-.5) + rad(1 + (m-1)*3)*SIN(th) - 0.01*(r-0.5) + Lb(2)/2 - 2*ro
          end do
       end do
       !X(2,28*Nb+1:29*Nb) = X(2,28*Nb+1:29*Nb) - 0.1*ro
       !X(2,29*Nb+1:30*Nb) = X(2,29*Nb+1:30*Nb) - 0.1*ro
       !X(1,27*Nb+1:28*Nb) = X(1,27*Nb+1:28*Nb) + 1.7*ro
       dSow = Lb(1)/(REAL(Nw)/2.)
       do i = 1,Nw/2
          wX(1,i) = REAL(i-1)*dSow
          wX(2,i) = 0.
       end do
       do i = Nw/2+1,Nw
          wX(1,i) = REAL(i-1-Nw/2)*dSow
          wX(2,i) = 0.99*Lb(2)
       end do
       wDx(1,:) = dSow/Soe(2)
       wDx(2,:) = 0.
       wD(:) = dSow/Soe(2)
       X(:,Ns+1:Np) = wX

     case(127)
   rad=ro; Tc = T0; Mc = M0
       !rad(Nm) = r1
       if (Nm.ne.15) stop 'bad Nm'
       do m = 1,15
          do n = 1,1
             call RANDOM_NUMBER(r)
             do i = 1,Nb
                th = 2.*Pi*REAL(i-1)/REAL(Nb)
!                X(1,i+Nb*(n - 1 + (m-1)*3)) = Lb/5.4*REAL(m-.5) + rad(n + (m-1)*5)*COS(th) + 0.2*(r-0.5)
                  X(1,i+Nb*(n - 1 + (m-1)*1)) = Lb(1)/15*REAL(m-.5) + rad(2)*COS(th) + 0*(r-0.5)
!                X(2,i+Nb*(n - 1+ (m-1)*3)) = Lb/10*REAL(n-.5) + rad(n + (m-1)*5)*SIN(th) - 0.2*(r-0.5) + 0*Lb/70.*REAL(m-1)
!                  X(2,i+Nb*(n - 1+ (m-1)*5)) = Lb/10*REAL(n-.5) + rad(2)*SIN(th)+ 0*Lb/70.*REAL(m-1) + Lb/4 - 0.3*(r-0.5)
              X(2,i+Nb*(n - 1+ (m-1)*1)) = rad(2)*SIN(th) - 0.25*(r-0.5) + 0.5*Lb(2)/2
             end do
          end do
       end do

       dSow = Lb(1)/(REAL(Nw)/2.)
       do i = 1,Nw/2
          wX(1,i) = REAL(i-1)*dSow
          wX(2,i) = 0.
       end do
       do i = Nw/2+1,Nw
          wX(1,i) = REAL(i-1-Nw/2)*dSow
          wX(2,i) = 0.5*Lb(2)
       end do
       wDx(1,:) = dSow/Soe(2)
       wDx(2,:) = 0.
       wD(:) = dSow/Soe(2)
       X(:,Ns+1:Np) = wX

    case(0)
       rad = ro; Tc = T0; Mc = M0
       !rad(Nm) = r1
       !Tc(Nm) = T1
       !Mc(Nm) = M1
       open(1,file='./2DCell/D/restart.in',form='UNFORMATTED')
       read(1)Np_l,Nb_l,Nm_l,lt0  ;     print *,Np_l,Nb_l,Nm_l
       read(1)time0,rad,dSo,Soe,len0,Tc,Mc,wX,wDx,wD
       print *,Np_l,Nb_l,Nm_l
       if (Np_l.ne.Np.or.Nb_l.ne.Nb.or.Nm_l.ne.Nm) stop 'bad restart'
       read(1)X
       close(1)
    case(-1)
       rad = ro; Tc = T0; Mc = M0
       rad(Nm) = r1
       Tc(Nm) = T1
       Mc(Nm) = M1
       open(1,file='D/restart.in',form='UNFORMATTED')
       read(1)Np_l,Nb_l,Nm_l,lt0
       read(1)time0,rad,dSo,Soe,len0,Tc,Mc,wX,wDx,wD
       if (Np_l.ne.Np.or.Nb_l.ne.Nb.or.Nm_l.ne.Nm) stop 'bad restart'
       read(1)X
       close(1)
       !       luk = MAXLOC(rad)
       !       rad = ro
       Tc = T0; Mc = M0
       !       rad(luk(1)) = r1
       !       Tc(luk(1)) = T1
       !       Mc(luk(1)) = M1
       open(1,file='RESTART-NEW-TM')
       do i = 1,Nm
          write(1,"(I5,3F20.10)")i,rad(i),Tc(i),Mc(i)
       end do
       close(1)

    case default
       stop 'bad Nm'
    end select

  END SUBROUTINE ic

END MODULE icmod
