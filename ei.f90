!PROGRAM expint
MODULE expint

  IMPLICIT none

!!$  do i = 1,1000
!!$     write(55,*)0.01*REAL(i),Ei(0.01*REAL(i))
!!$  end do

  integer, parameter   :: Nei = 10000
  real, parameter      :: Ei_xmax = 20.
  real                 :: Ei_dx
  real, dimension(Nei) :: Ei_tab,Ei_x
  
CONTAINS
  
  SUBROUTINE initEi

    real    :: x
    integer :: i

    Ei_dx = Ei_xmax/REAL(Nei)
    do i = 1,Nei
       Ei_x(i) = Ei_dx*REAL(i)
       Ei_tab(i) = Ei_ex( Ei_x(i) )
    end do

  END SUBROUTINE initEi

  SUBROUTINE testEi

    real    :: x
    integer :: i

    do i = 1,5*Nei
       x = REAL(i)*Ei_dx/5.
       write(88,*)x,Ei(x),Ei_ex(x)
    end do

  END SUBROUTINE testEi


  FUNCTION Ei(x) RESULT(y)
    real  :: x,y

    integer :: j

    j = MIN(INT(x/Ei_dx),Nei-1)
!    print *,x,j,Ei_dx,Ei_x(Nei)
    y = Ei_tab(j) + (Ei_tab(j+1)-Ei_tab(j))*(x-Ei_x(j))/Ei_dx
!    write(55,*) x

  END FUNCTION Ei

  FUNCTION Ei_ex(x) RESULT (y)
    real  :: x,y

    real, parameter :: a0 = -0.57721566
    real, parameter :: a1 =  0.99999193
    real, parameter :: a2 = -0.24991055
    real, parameter :: a3 =  0.05519968
    real, parameter :: a4 = -0.00976004
    real, parameter :: a5 =  0.00107857

    real, parameter :: AA1 = 2.334733  
    real, parameter :: AA2 = 0.250621
    real, parameter :: BB1 = 3.330657  
    real, parameter :: BB2 = 1.681534  

    if (x.le.1) then
       y = a0+x*(a1+x*(a2+x*(a3+x*(a4+x*a5)))) - LOG(x)
    else
       y = (AA2+x*(AA1+x))/(BB2+x*(BB1+x))/(x*EXP(x))
    end if

  END FUNCTION Ei_ex

!END PROGRAM expint
END MODULE expint
