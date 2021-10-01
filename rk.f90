MODULE rkmod

  USE prms
  USE data
  USE stokes

  IMPLICIT none

  real                       :: fac_3,fac_6,fac_2

CONTAINS

  SUBROUTINE initRK

    fac_3 = 1./3.*dt
    fac_2 = .5*dt
    fac_6 = 1./6.*dt

  END SUBROUTINE initRK

  SUBROUTINE rk2(time,lt)

    real,dimension(2,Np)    :: dXac,dXcf,dXmc,Xc,UU
    real                    :: time_local
    real                    :: time
    integer                 :: lt

    character(50)           :: fn

! temporary
    integer :: i

        time_local = time
 !     print*, "call rhs"
    call rhs(X,dXac,time_local,lt)
 !    print*, "call meanu"
    call meanU(X,dXmc,time_local)
    Xc = X + 0.5*dt*(dXac+dXmc)

    !open(1,file='CellVel.xdir.rk',form='FORMATTED')
    !Do i = 1, Np
        !write(1,*) dXac(1,i)+dXmc(1,i)
    !End Do

    close(1)
    time_local = time + dt/2.
    call rhs(Xc,dXac,time_local,lt)
    call meanU(Xc,dXmc,time_local)
    X = X + dt*(dXac+dXmc)

    time = time + dt


    !time_local = time
    !call rhs(X,dXac,time_local,lt)
    !Xc = X + 0.5*dt*dXac
    !!call meanU(X,dXmc)
    !!Xc = X + 0.5*dt*(dXac+dXmc)

    !time_local = time + dt/2.
    !call rhs(Xc,dXac,time_local,lt)
    !X = X + dt*dXac
    !!call meanU(Xc,dXmc)
    !!X = X + dt*(dXac+dXmc)

!time = time + dt

        !open(1,file='CellVel.xdir.rk',form='FORMATTED')
    !Do i = 1, Np
        !write(1,*) UU(1,i)
    !End Do

    !close(1)



  END SUBROUTINE rk2

END MODULE rkmod

