MODULE data

  USE prms

  IMPLICIT none

  real, allocatable, dimension(:,:)    :: X      ! control points
  real, allocatable, dimension(:)      :: Soe    ! every point  ref coord
  real, allocatable, dimension(:)      :: dSo    ! dso of all ctrl points
  real, allocatable, dimension(:)      :: rad    ! initial radius
  real, allocatable, dimension(:)      :: Tc,Mc  ! elastiticy arrays
  integer, allocatable, dimension(:)   :: cell   ! cell number for any of Np points 
  real, allocatable, dimension(:,:)    :: wX     ! wall initial loc
  real, allocatable, dimension(:,:)    :: wDx    ! derivatives of wall 
  real, allocatable, dimension(:)      :: wD     ! stretch of wall
  real, allocatable, dimension(:)      :: len0   ! unstretched length


CONTAINS

  SUBROUTINE initdata

    allocate(X(2,Np),Soe(Nb),dSo(Np),rad(Nm),Tc(Nm),Mc(Nm))
    allocate(wX(2,Nw),wDx(2,Nw),wD(Nw))
    allocate(cell(Np),len0(Nm))

  END SUBROUTINE initdata

  SUBROUTINE finaldata

    deallocate(X,Soe,dSo,rad,Tc,Mc)
    deallocate(wX,wDx,wD,len0,cell)

  END SUBROUTINE finaldata

END MODULE data
