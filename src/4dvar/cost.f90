! **************************************************************************
! ... cost.f90
! ... Quim Ballabrera, January 2018
! ...
! ... System Dynamics
! ... Adjoint methods applied to the Lorenz model
! ...
! ... The FUNCTION JCOST calculates the cost function:
! ...
! ...            Jcost = 1/2 * Sum_i=1^p  [y^o(t_i) - M(xo,to,ti)]^2
! ... 
! ... as a function of the intial condition, xo.
! ...
! ... The subroutine runs the lorenz model and accounts for the 
! ... difference against observations.
! ...
! ... The SUBROUTINE DCOST calculates the derivative of the cost function
! ... by respect the initial conditions:
! ...
! ...                         DCOST(i) = d Jcost / d xo(i)
! ...
! ... by using the adjoint of the Lorenz model.
! ...
! **************************************************************************

real(dp) function Jcost (xo)

use sysdyn
use adjoint

implicit none

real(dp), dimension(NVE), intent(in)          :: xo

! ... Local variables:
! ...
integer i,step,nsteps,ncount
real(dp) to,tf,dt,t,innov,Jb
real(dp), dimension(NVE)                      :: x,xn,dd

! ... Model time interval:
! ...
dt = lorenz_dt
to = cost_to
tf = cost_tf
nsteps = nint((tf-to)/dt)

if (verbose) then
  write(*,*) 
  write(*,*) '== COST FUNCTION == '
  write(*,*) 'Initial time  = ', to
  write(*,*) 'Final   time  = ', tf
  write(*,*) 'Time interval = ', dt
  write(*,*) 'Model steps   = ', nsteps
endif

! ... Find the index:
! ...
do step=1,obs%Nobs
  if ((obs%time(step).ge.to).or. &
      (abs(obs%time(step)-to).lt.0.01*dt)) then
    obs%p = step
    exit
  endif
enddo

! ... Evaluate cost function:
! ...
Jcost = zero
ncount = 0

! ... Initial condition contribution
! ...
if (abs(obs%time(obs%p)-to).lt.0.01*dt) then
  if (verbose) write(*,*) 'Cost contribution, time = ', to
  do i=1,NVE
    if (obs%iobs(i).eq.1) then
      innov  = obs%y(i,obs%p) - xo(i)
      Jcost  = Jcost  + innov*innov/r2(i)
      ncount = ncount + 1
    endif
  enddo
  obs%p = obs%p + 1
endif

! ... Background field:
! ...
if (alpha.ge.1.d-8) then
  dd = xb(:) - xo(:)
  xn = MATMUL(Binv,dd)
  Jb = alpha*DOT_PRODUCT(dd,xn)
  if (verbose) write(*,*) 'Adding Background term', Jb
  Jcost = Jcost + Jb
endif

! ... Now run the model
! ...
x = xo
t = to
do step=1,nsteps

  call rk4 (lorenz_n,x,t,dt,xn,lorenz_rhs)
  x(:) = xn(:)
  t    = t + dt

  if (obs%p.le.obs%Nobs) then
    if (step.eq.obs%step(obs%p)) then
      !print*, 'Cost contribution, time = ', t
      do i=1,NVE
        if (obs%iobs(i).eq.1) then
          innov  = obs%y(i,obs%p) - x(i)
          Jcost  = Jcost  + innov*innov/r2(i)
          ncount = ncount + obs%iobs(i)
        endif
      enddo
      obs%p = obs%p + 1
    endif
  endif

enddo

if (verbose) then
  write(*,*) 'Number observ = ', ncount
endif

Jcost = half * Jcost

end function Jcost
! ...
! =================================================================================
! ...
subroutine dCost (xo,gradx)

use sysdyn
use adjoint

implicit none

real(dp), dimension(NVE), intent(in)          :: xo
real(dp), dimension(NVE), intent(out)         :: gradx

! ... Local variables:
! ...
integer iu,i,step,nsteps
real(dp) to,tf,dt,t
real(dp), dimension(NVE)                      :: x,xn,xr,dd
real(dp), dimension(NVE)                      :: gradJ
real(dp), dimension(NVE)                      :: p

! ... Model time interval:
! ...
dt = lorenz_dt
to = cost_to
tf = cost_tf
nsteps = nint((tf-to)/dt)

if (verbose) then
  write(*,*)
  write(*,*) 'Initial time  = ', to
  write(*,*) 'Final   time  = ', tf
  write(*,*) 'Time interval = ', dt
  write(*,*) 'Model steps   = ', nsteps
endif

! ... Find the obsevation pointer value (observation for the 
! ... interval end):
! ...
do step=1,obs%Nobs
  if ((obs%time(step).ge.tf).or. &
      (abs(obs%time(step)-tf).lt.0.01*dt)) then
    obs%p = step
    exit
  endif
enddo
if (obs%time(obs%p).gt.tf+0.01*dt) obs%p = obs%p - 1

! ... Run the forward model, saving each time step solution
! ... The Adjoint code is based on the linearization
! ... around each solution state
! ...
iu = unitfree()
open(iu,file='adj.tmp',status='unknown',form='unformatted', &
        access='direct',recl=dp*(NVE+1))

t = to
x = xo
!print*, 'saving step ', 0, t
write(iu,rec=1) t, x

do step=1,nsteps
  call rk4 (lorenz_n,x,t,dt,xn,lorenz_rhs)
  x = xn
  t = t + dt
  !print*, 'saving step ', step, t
  write(iu,rec=step+1) t, x
enddo
close(iu)

! ... Run the model backwards:
! ...
gradJ(:) = zero
p(:)     = zero

open(iu,file='adj.tmp',status='old',form='unformatted', &
        access='direct',recl=dp*(NVE+1))

read(iu,rec=nsteps+1) t, xr(:)
do step=nsteps,1,-1
  if (obs%p.gt.0) then
    if (abs(obs%time(obs%p)-t).lt.0.01*dt) then
      !print*, 'Assimilating observation ', obs%p, obs%time(obs%p)
      do i=1,NVE
        gradJ(i) = - obs%iobs(i)*(obs%y(i,obs%p)-xr(i))/r2(i)
      enddo
      obs%p = obs%p - 1
    endif
  endif
  read(iu,rec=step) t, lorenz_ref(:)
  call rk4_adj(lorenz_n,p,t,dt,gradJ,lorenz_adj)
enddo

! ... Initial time step
obs%p = 1
read(iu,rec=1) t, xr(:)
if (abs(obs%time(obs%p)-to).lt.0.01*dt) then
  if (verbose) print*, 'Assimilating observation initial time', to
  do i=1,NVE
    gradJ(i) = - obs%iobs(i)*(obs%y(i,obs%p)-xr(i))/r2(i)
    p(i) = p(i) + gradJ(i)
  enddo
endif

! ... Adding the Background term
! ...
if (alpha.ge.1D-8) then
  dd   = xb(:) - xr(:)
  p(:) = p(:) - alpha*MATMUL(Binv,dd)
endif

close(iu,status='DELETE')

gradx(:) = p(:)

end subroutine dCost
