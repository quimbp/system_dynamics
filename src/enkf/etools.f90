! **************************************************************************
! ... etools.f90
! ... Quim Ballabrera, February 2013
! ...
! ... System Dynamics:
! ... ETOOLS
! ... Ensemble tools
! ...
! **************************************************************************


subroutine emean (X,m,n,xm)

implicit none

integer, intent(in)                       :: m,n
real(kind=8), dimension(m,n), intent(in)  :: X
real(kind=8), dimension(m), intent(out)   :: xm

integer i
real(kind=8) In

In = 1.0D0/DBLE(n)
do i=1,m
  xm(i) = SUM(X(i,:))*In
enddo

return
end subroutine emean
! ...
! ========================================================================
! ...
subroutine Hforward (x,n,iobs,p,y)

use sysdyn

! ...
! ... This version of Hforward works with the iobs(1:n) vector
! ... of 1 (observed component) and 0 (non-observed component).
! ... p = SUM(iobs) is the number of observations.
! ...

integer, intent(in)                     :: n,p
integer, dimension(p), intent(in)       :: iobs
real(kind=8), dimension(n), intent(in)  :: x
real(kind=8), dimension(p), intent(out) :: y

! ... Local variables:
! ...
integer i,ii

ii = 0
do i=1,n
  if (Iobs(i).gt.0) then
    ii = ii + 1
    y(ii) = x(i)
  endif
enddo
if (ii.NE.p) call stop_error (1,'Incompatible number of observations')

end subroutine Hforward

