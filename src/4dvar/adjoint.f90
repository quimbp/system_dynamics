! **************************************************************************
! ... adjoint.f90
! ... Quim Ballabrera, January 2018
! ...
! ... System Dynamics
! ... Adjoint methods applied to the Lorenz model
! ...
! **************************************************************************

module adjoint

use sysdyn

implicit none

integer, parameter                            :: NVE = lorenz_n
type(obs_type)                                :: obs

logical                                       :: verbose = .true.
real(dp)                                      :: cost_to
real(dp)                                      :: cost_tf
real(dp), dimension(NVE)                      :: r2

! ... Background covariance factor
! ... If zero, no background term is used
! ... 
real(dp)                                      :: alpha   = zero
real(dp), dimension(NVE)                      :: xb
real(dp), dimension(NVE,NVE)                  :: Bcov
real(dp), dimension(NVE,NVE)                  :: Binv

contains

! ...
! =====================================================================
! ...
subroutine rk4_adj (n,ax,t,h,axn,aderivs)

! ... Given value fort the adjoint variables x_adj and xn_adj,
! ... known at t, use the adjoint fourth-order Runge-Kutta method to 
! ... backpropagate the adjoint over an interval h and return the 
! ... incremented variable x_adj(1:n). 
!  ...
! ... The user supplies the subroutine aderivs(n,t,Ax,Adxdt), which 
! ... returns the adjoint of the derivatives dxdt at t.
! ... Based on rk4.f of the Numerical Recipes
! ...
! ... Fixed system reference in LORENZ_REF

integer, intent(in)                            :: n
real(dp), intent(in)                           :: t,h
real(dp), dimension(n), intent(inout)          :: ax(n)
real(dp), dimension(n), intent(inout)          :: axn(n)
external                                       :: aderivs

! ... Local variables:
! ...
real(dp), dimension(n)                         :: ak1,ak2,ak3,ak4
real(dp), dimension(n)                         :: ax1,ax2,ax3,ax4




! ... xn(:) = x(:) + (k1 + 2*k2 + 2*k3 + k4)/6
! ...
 ax(:) = ax(:) + axn(:)
ak1(:) = axn(:)/6d0
ak2(:) = axn(:)/3D0
ak3(:) = axn(:)/3D0
ak4(:) = axn(:)/6D0
axn(:) = zero

! ... x4 = x + k3
! ... call rsh(n,t+h,x4,k4)
! ... k4 = h*k4
! ...
ak4(:) = h*ak4(:)
ax4(:) = zero
call aderivs(n,t+h,ax4,ak4)
 ax(:) =  ax(:) + ax4(:)
ak3(:) = ak3(:) + ax4(:)
ax4(:) = zero

! ... x3 = x + 0.5*k2
! ... call rsh(n,t+h/2,x3,k3)
! ... k3 = h*k3
! ...
ak3(:) = h*ak3(:)
ax3(:) = zero
call aderivs(n,t+half*h,ax3,ak3)
 ax(:) =  ax(:) + ax3(:)
ak2(:) = ak2(:) + half*ax3(:)
ax3(:) = zero

! ... x2 = x + 0.5*k1
! ... call rhs(n,t+h/2,x2,k2)
! ... k2 = h*k2
! ...
ak2(:) = h*ak2(:)
ax2(:) = zero
call aderivs(n,t+half*h,ax2,ak2)
 ax(:) =  ax(:) + ax2(:)
ak1(:) = ak1(:) + half*ax2(:)
ax2(:) = zero

! ... x1 = x
! ... call rhs(n,t,x1,k1)
! ... k1 = h*k1
! ...
ak1(:) = h*ak1(:)
ax1(:) = zero
call aderivs(n,t,ax1,ak1)
 ax(:) = ax(:) + ax1(:)

return
end subroutine rk4_adj
! ...
! =====================================================================
! ...

end module adjoint

