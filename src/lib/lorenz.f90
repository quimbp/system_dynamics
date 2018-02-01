! *****************************************************************************
! ... lorenz.f90
! ... The 1963 Lorenz model.
! ... COSMO project
! ... Quim Ballabrera, March 2017
! *****************************************************************************

module lorenz

use types, only: dp

implicit none

private
public lorenz_n
public lorenz_dt
public lorenz_to
public lorenz_tf
public lorenz_sigma,lorenz_beta,lorenz_rho
public lorenz_ini
public lorenz_ref
public lorenz_rhs
public lorenz_ltm
public lorenz_adj

integer, parameter               :: lorenz_n     = 3
real(dp)                         :: lorenz_to    = 0.0_dp
real(dp)                         :: lorenz_tf    = 20.0_dp
real(dp)                         :: lorenz_dt    = 0.01_dp
real(dp)                         :: lorenz_sigma = 10.0_dp
real(dp)                         :: lorenz_beta  = 8.0_dp/3.0_dp
real(dp)                         :: lorenz_rho   = 28.0_dp

real(dp), dimension(3) :: lorenz_ini = [1.508870_dp,-1.531271_dp,25.46091_dp]
real(dp), dimension(3) :: lorenz_ref 

contains
! ...
! =============================================================================
! ...
subroutine lorenz_rhs (n,t,x,dxdt)

integer, intent(in)                   :: n
real(dp), intent(in)                  :: t
real(dp), dimension(n), intent(in)    :: x
real(dp), dimension(n), intent(out)   :: dxdt

dxdt(1) = lorenz_sigma*(x(2)-x(1))
dxdt(2) = x(1)*(lorenz_rho-x(3)) - x(2)
dxdt(3) = x(1)*x(2) - lorenz_beta*x(3)

return
end subroutine lorenz_rhs

! ...
! =============================================================================
! ...
subroutine lorenz_ltm (n,t,x,dxdt)

integer, intent(in)                   :: n
real(dp), intent(in)                  :: t
real(dp), dimension(n), intent(in)    :: x
real(dp), dimension(n), intent(out)   :: dxdt

dxdt(1) = lorenz_sigma*(x(2)-x(1))
dxdt(2) = x(1)*(lorenz_rho-lorenz_ref(3)) - lorenz_ref(1)*x(3) - x(2)
dxdt(3) = lorenz_ref(1)*x(2) + x(1)*lorenz_ref(2) - lorenz_beta*x(3)

return
end subroutine lorenz_ltm
! ...
! =============================================================================
! ...
subroutine lorenz_adj (n,t,x_adj,dxdt_adj)

integer, intent(in)                   :: n
real(dp), intent(in)                  :: t
real(dp), dimension(n), intent(inout) :: x_adj
real(dp), dimension(n), intent(inout) :: dxdt_adj

!dxdt(3) = lorenz_ref(1)*x(2) + x(1)*lorenz_ref(2) - lorenz_beta*x(3)
x_adj(1) = x_adj(1) + lorenz_ref(2)*dxdt_adj(3)
x_adj(2) = x_adj(2) + lorenz_ref(1)*dxdt_adj(3)
x_adj(3) = x_adj(3) -   lorenz_beta*dxdt_adj(3)
dxdt_adj(3) = 0.0_dp

!dxdt(2) = x(1)*(lorenz_rho-lorenz_ref(3)) - lorenz_ref(1)*x(3) - x(2)
x_adj(1) = x_adj(1) + (lorenz_rho-lorenz_ref(3))*dxdt_adj(2)
x_adj(2) = x_adj(2) -                            dxdt_adj(2)
x_adj(3) = x_adj(3) -              lorenz_ref(1)*dxdt_adj(2)
dxdt_adj(2) = 0.0_dp

!dxdt(1) = lorenz_sigma*(x(2)-x(1))
x_adj(1) = x_adj(1) - lorenz_sigma*dxdt_adj(1)
x_adj(2) = x_adj(2) + lorenz_sigma*dxdt_adj(1)
dxdt_adj(1) = 0.0_dp

return
end subroutine lorenz_adj

end module lorenz
