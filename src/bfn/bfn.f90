! **************************************************************************
! ... nudging.f90
! ... Quim Ballabrera, February 2013
! ...
! ... System Dynamics
! ... Nudging applied to the Lorenz model
! ...
! **************************************************************************

program nudging

use sysdyn

implicit none

integer, parameter                          :: NVE = lorenz_n
type(obs_type)                              :: obs

logical default,Zexist

integer i,ii,step,nstep,iu,k,kk
real(dp) T,Jc,a,b,tt,apply

integer,  dimension(NVE)                    :: Iobs
real(dp), dimension(NVE)                    :: Xo,X,Xn,Xm
real(dp), dimension(NVE)                    :: Yo
character(len=180) nfile,ifile,ofile,rpath,rfile,bfile

! ... Lorenz model parameters
! ...
real(dp)                                    :: to,tf,dt
real(dp)                                    :: sigma
real(dp)                                    :: rho
real(dp)                                    :: beta

! ... Nudging parameters
! ...
real(dp) mux,muy,muz
real(dp) Tx,Ty,Tz

NAMELIST/namnudging/ifile,ofile,          &
                    rpath,rfile,          &
                    bfile,                &
                    mux,muy,muz,          &
                    Tx,Ty,Tz,             &
                    Iobs,                 &
                    sigma,rho,beta,       &
                    to,dt,tf

to    = lorenz_to
dt    = lorenz_dt
tf    = lorenz_tf
sigma = lorenz_sigma
rho   = lorenz_rho
beta  = lorenz_beta

nfile = ''
CALL getarg(1,nfile)
inquire(file=nfile,exist=Zexist)

if (len_trim(nfile).gt.0) then
  if (.not.Zexist) stop 'Namelist file not found'
endif

if (Zexist) then
  write(*,*) 'Reading ', trim(nfile)
  open(10,FILE=nfile,STATUS='old')
  rewind(10)
  read(10,namnudging)
  close(10)
  default = .false.
else
  default = .true.
endif

IF (default) THEN
  write(*,*) 'Default nudging parameters'

  ! ... Starting model value:
  ! ...
  ifile = 'nudging.ini'
  to    =  0.0D0
  Xo(:) = (/ -1.483,  6.517, 30.774 /)

  iu = unitfree()
  open(iu,file=trim(ifile),status='unknown')
  write(iu,*) to, Xo(:)
  close(iu)

  ! ... Observation file
  ofile = './obs_noise_02_sampl_010.dat'

  ! ... Output file
  rpath = './'
  rfile = './forward.dat'
  bfile = './backward.dat'

  ! ... Nudging coefficient:
  ! ...
  mux =  100D0
  muy =  100D0
  muz =  100D0
  Tx  =  0.10
  Ty  =  0.10
  Tz  =  0.10

  ! ... Model parameters:
  ! ...
  sigma = lorenz_sigma
  rho   = lorenz_rho
  beta  = lorenz_beta

  ! ... Model simulation time
  ! ...
  tf = 20.0D0
  dt = lorenz_dt

  ! ... What do we assimilate
  ! ...
  Iobs(:) = 1   ! All obs are assimilated
  !Iobs(:) = (/1, 1, 0/) 

  
  ! ... If file nudging.namelist does not exist, we
  ! ... create a copy that can be later modified ....
  ! ...
  inquire(file='bfn.namelist',exist=Zexist)
  if (.not.Zexist) then
    open(20,file='bfn.namelist',status='new')
    write(20,nml=namnudging)
    close(20)
  endif

endif

lorenz_sigma = sigma
lorenz_rho   = rho
lorenz_beta  = beta

! ... Read the intial solution
! ...
write(*,*) 'Reading initial solution from ', trim(ifile)
iu = unitfree()
open(iu,file=trim(ifile),status='old')
rewind(iu)
read(iu,*) tt, Xo(:)
close(iu)


! ... Read the observations
! ...
obs = obs_read (ofile)

write(*,*) 
write(*,*) 'Nobs                   = ', obs%Nobs
write(*,*) 'First time record      = ', obs%time(1)
write(*,*) 'Last time record       = ', obs%time(obs%Nobs)
write(*,*) 'Iobs                   = ', Iobs(:)
write(*,*) 'X nudging coefficient  = ', mux
write(*,*) 'Y nudging coefficient  = ', muy
write(*,*) 'Z nudging coefficient  = ', muz
write(*,*) 'Tx                     = ', Tx
write(*,*) 'Ty                     = ', Ty
write(*,*) 'Tz                     = ', Tz

! ... Fill the h vector:
! ...
obs%h(:) = .false.
where (Iobs(:).ge.1) obs%h = .true.

nstep = NINT((tf-to)/dt)

write(*,*) 
write(*,*) 'Model initial time = ', to
write(*,*) 'Model final time   = ', tf
write(*,*) 'Model time step    = ', dt
write(*,*) 'Number time steps  = ', nstep + 1
write(*,*) 'Model parameters'
write(*,*) 'sigma              = ', lorenz_sigma
write(*,*) 'rho                = ', lorenz_rho   
write(*,*) 'beta               = ', lorenz_beta  

! ... Fill the step vector:
! ...
obs%step(:) = -1
kk = 1
do i=1,nstep+1
 tt = to + (i-1)*dt
 if (abs(obs%time(kk)-tt).le.0.1*dt) then
   obs%step(kk) = i
   kk = kk + 1
   if (kk.gt.obs%Nobs) exit
 endif
enddo
do i=1,obs%Nobs
  print*, i, obs%time(i), obs%step(i)
enddo

! ... Make sure that the observation index
! ... points to the first observation
do i=1,obs%Nobs
  if (obs%time(i).ge.to) then
    obs%p = i
    exit
  endif
enddo
write(*,*)
write(*,*) 'First observation ', obs%p, obs%time(obs%p)
write(*,*)

! ... Nudging data assimilation:
! ...
write(*,*) 'Output file ', trim(filename_append(rpath,rfile))
open(77,file=trim(filename_append(rpath,rfile)),status='unknown')
rewind(77)

t = to
X = Xo
write(77,'(4F10.4)') t, X

do step=1,nstep

  ! ... Assimilation is performed all the time:
  ! ... dx/dt = Physics + mu*SUM (w(t'-t)*(y(t')-x(t))
  ! ...
  CALL rk4 (lorenz_n,X,t,dt,XN,rhwnud)

  X = XN
  t = t + dt
  write(77,'(4F10.4)') t, X
ENDDO
write(*,*) 'closing file ', trim(filename_append(rpath,rfile))
CLOSE(77)

! ... Backward simulation:
! ...
dt = -dt

write(*,*) 'Output file ', trim(filename_append(rpath,bfile))
open(77,file=trim(filename_append(rpath,bfile)),status='unknown')
rewind(77)

write(77,'(4F10.4)') t, X
do step=1,nstep

  ! ... Assimilation is performed all the time:
  ! ... dx/dt = Physics + mu*SUM (w(t'-t)*(y(t')-x(t))
  ! ...
  CALL rk4 (lorenz_n,X,t,dt,XN,back_rhwnud)

  X = XN
  t = t + dt
  write(77,'(4F10.4)') t, X
ENDDO
write(*,*) 'closing file ', trim(filename_append(rpath,bfile))
close(77)


open(77,file='bfn.out',status='unknown')
write(77,'(4F10.4)') t, X
close(77)



call stop_error(0,'Nudging Ok')

contains
! ...
! ==========================================================================
! ...
  subroutine rhnud (n,t,x,dxdt)

  implicit none

  integer, intent(in)     :: n
  real(dp), intent(in)    :: t
  real(dp), intent(in)    :: x(n)
  real(dp), intent(out)   :: dxdt(n)

  ! ... Local variables
  ! ...
  real(dp), dimension(n)  :: nr

  ! ... Lorenz model equations
  ! ...
  dxdt(1) = lorenz_sigma*(X(2)-X(1))  
  dxdt(2) = (lorenz_rho-X(3))*X(1) - X(2)
  dxdt(3) = X(1)*X(2) - lorenz_beta*X(3)

  ! ... Newtonian Relaxation:
  ! ... apply : is 0 or 1 depending if obs is available at that time
  ! ... obs%h : is True for the variables being assimilated, false otherwise
  ! ...
  nr(:) = zero
  if (obs%h(1)) nr(1) = apply*mux*(Yo(1)-X(1))
  if (obs%h(2)) nr(2) = apply*muy*(Yo(2)-X(2))
  if (obs%h(3)) nr(3) = apply*muz*(Yo(3)-X(3))

  dxdt(:) = dxdt(:) + nr(:)

  return
  end subroutine rhnud
! ...
! ==========================================================================
! ...
  subroutine rhwnud (n,t,x,dxdt)

  implicit none

  integer, intent(in)     :: n
  real(dp), intent(in)    :: t
  real(dp), intent(in)    :: x(n)
  real(dp), intent(out)   :: dxdt(n)
 
  ! ... Local variables
  ! ...
  real(dp) xsum,ysum,zsum
  real(dp), dimension(n)  :: nr

  ! ... Lorenz model equations
  ! ...
  dxdt(1) = sigma*(X(2)-X(1))  
  dxdt(2) = (rho-X(3))*X(1) - X(2)
  dxdt(3) = X(1)*X(2) - beta*X(3)

  ! ... Newtonian Relaxation:
  ! ... obs%h : is True for the variables being assimilated, false otherwise
  ! ...
  call weight_innov (t,x,xsum,ysum,zsum)

  nr(:) = zero
  if (obs%h(1)) nr(1) = mux*xsum
  if (obs%h(2)) nr(2) = muy*ysum
  if (obs%h(3)) nr(3) = muz*zsum

  dxdt(:) = dxdt(:) + nr(:)

  return
  end subroutine rhwnud
! ...
! ========================================================================
! ...
  subroutine weight_innov (Time,X,xsum,ysum,zsum)

  implicit none

  real(dp), intent(in)                   :: Time
  real(dp), dimension(NVE), intent(in)   :: X
  real(dp), intent(out)                  :: xsum,ysum,zsum

  integer i
  real(dp) TT,argx,argy,argz,wx,wy,wz

  xsum = 0.0D0
  ysum = 0.0D0
  zsum = 0.0D0
  do i=1,obs%Nobs
    TT   = ABS(Time - obs%time(i))
    argx = TT/Tx
    argy = TT/Ty
    argz = TT/Tz
    wx = 0D0
    wy = 0D0
    wz = 0D0
    IF (argx.LT.10) wx = exp(-argx*argx)
    IF (argy.LT.10) wy = exp(-argy*argy)
    IF (argz.LT.10) wz = exp(-argz*argz)
  
    if (Time.gt.obs%time(i)+0.001D0) then
      wx = 0D0
      wy = 0D0
      wz = 0D0
    endif
  
    !print*, Time, i, obs%time(i), wx
    xsum = xsum + wx*(obs%y(1,i)-X(1))
    ysum = ysum + wy*(obs%y(2,i)-X(2))
    zsum = zsum + wz*(obs%y(3,i)-X(3))
  enddo

  return
  end subroutine weight_innov
! ...
! ==========================================================================
! ...
  subroutine back_rhwnud (n,t,x,dxdt)

  implicit none

  integer, intent(in)     :: n
  real(dp), intent(in)    :: t
  real(dp), intent(in)    :: x(n)
  real(dp), intent(out)   :: dxdt(n)

  ! ... Local variables
  ! ...
  real(dp) xsum,ysum,zsum
  real(dp), dimension(n)  :: nr

  ! ... Lorenz model equations
  ! ...
  dxdt(1) = sigma*(X(2)-X(1))
  dxdt(2) = (rho-X(3))*X(1) - X(2)
  dxdt(3) = X(1)*X(2) - beta*X(3)

  ! ... Newtonian Relaxation:
  ! ... obs%h : is True for the variables being assimilated, false otherwise
  ! ...
  call back_weight_innov (t,x,xsum,ysum,zsum)

  nr(:) = zero
  if (obs%h(1)) nr(1) = mux*xsum
  if (obs%h(2)) nr(2) = muy*ysum
  if (obs%h(3)) nr(3) = muz*zsum

  dxdt(:) = dxdt(:) - nr(:)

  return
  end subroutine back_rhwnud
! ...
! ========================================================================
! ...
  subroutine back_weight_innov (Time,X,xsum,ysum,zsum)

  implicit none

  real(dp), intent(in)                   :: Time
  real(dp), dimension(NVE), intent(in)   :: X
  real(dp), intent(out)                  :: xsum,ysum,zsum

  integer i
  real(dp) TT,argx,argy,argz,wx,wy,wz

  xsum = 0.0D0
  ysum = 0.0D0
  zsum = 0.0D0
  do i=1,obs%Nobs
    TT   = ABS(Time - obs%time(i))
    argx = TT/Tx
    argy = TT/Ty
    argz = TT/Tz
    wx = 0D0
    wy = 0D0
    wz = 0D0
    IF (argx.LT.10) wx = exp(-argx*argx)
    IF (argy.LT.10) wy = exp(-argy*argy)
    IF (argz.LT.10) wz = exp(-argz*argz)
  
    !if (Time.lt.obs%time(i)+0.001D0) then
    !  wx = 0D0
    !  wy = 0D0
    !  wz = 0D0
    !endif
  
    !print*, Time, i, obs%time(i), wx
    xsum = xsum + wx*(obs%y(1,i)-X(1))
    ysum = ysum + wy*(obs%y(2,i)-X(2))
    zsum = zsum + wz*(obs%y(3,i)-X(3))
  enddo

  return
  end subroutine back_weight_innov

end program nudging

