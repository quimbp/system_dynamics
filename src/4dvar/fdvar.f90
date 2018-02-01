program fdvar

use sysdyn
use adjoint

implicit none

logical Zexist,default
integer iu,iter,step,nsteps,cycle,ncycles,i,j,nn
integer, dimension(NVE)                              :: Iobs
real(dp), dimension(NVE)                             :: xo,x,xn
real(dp), dimension(NVE)                             :: gradJ
real(dp) to,tf,time,Jo,Jmin,gtol,cdt
real(dp) sigma,rho,beta
real(dp) r2x,r2y,r2z
character(len=180) nfile,ifile,ofile,bfile,tfile,efile,afile,opath,mfile

real(dp) Jcost
external Jcost
external Dcost

namelist /NAMFDVAR/ifile,    &
                   ofile,    &
                   Iobs,     &
                   r2x,      &
                   r2y,      &
                   r2z,      &
                   bfile,    &
                   mfile,    &
                   alpha,    &
                   opath,    &
                   to,       &
                   tf,       &
                   ncycles,  &
                   afile,    &
                   tfile,    &
                   efile,    &
                   sigma,    &
                   rho,      &
                   beta  
               
call header()

! ... Default options
! ...
Iobs(:) = 1
ncycles = 10
opath   = './'
sigma = lorenz_sigma
rho   = lorenz_rho  
beta  = lorenz_beta 
r2x   = 2.0d0
r2y   = 2.0d0
r2z   = 2.0d0
bfile = ''
mfile = ''
alpha = zero

call getarg(1,nfile)
inquire(file=nfile,exist=Zexist)

if (len_trim(nfile).gt.0) then
  if (Zexist) then
    write(*,*) 'Reading ', trim(nfile)
    open(10,file=nfile,status='old')
    read(10,namfdvar)
    close(10)
    default = .false.
  else
    call stop_error(1,'Namelist file not found')
  endif
else
  default = .true.
endif

if (default) then
  ofile   = '../../work/obs_noise_02_sampl_010.dat'
  ifile   = '../../work/fdvar_ini.dat'
  afile   = '../../work/fdvar_opt.dat'
  tfile   = '../../work/fdvar_noise_02_sampl_010.dat'
  efile   = '../../work/fdvar_end.dat'
  bfile   = ''
  mfile   = ''
  to      = zero
  tf      = one

  ! ... If file enkf.namelist does not exist, we
  ! ... create a copy that can be later modified ....
  ! ...
  inquire(file='fdvar.namelist',exist=Zexist)
  if (.not.Zexist) then
    iu = unitfree()
    open(iu,file='fdvar.namelist',status='new')
    write(iu,nml=namfdvar)
    close(iu)
  endif
endif

lorenz_sigma = sigma
lorenz_rho   = rho  
lorenz_beta  = beta 

write(*,*)
write(*,*) 'Cost  initial time = ', to
write(*,*) 'Cost  final time   = ', tf
write(*,*) 'Model time step    = ', lorenz_dt
write(*,*) 'Model parameters'
write(*,*) 'sigma              = ', lorenz_sigma
write(*,*) 'rho                = ', lorenz_rho
write(*,*) 'beta               = ', lorenz_beta

! ... Read the observations
! ...
obs = obs_read(ofile)
Iobs(:) = 1
obs%iobs(:) = Iobs(:)

! ... Observation error:
! ...
r2(1) = r2x
r2(2) = r2y
r2(3) = r2z

! ... Background error covariance:
! ...
if (len_trim(bfile).gt.0) then
  write(*,*) 
  write(*,*) 'Reading error covariance ', trim(bfile)

  iu = unitfree()
  open(iu,file=bfile,status='old')
  read(iu,*) i, j
  if (i.ne.NVE) call stop_error(1,'Invalid covariance matrix')
  if (j.ne.NVE) call stop_error(1,'Invalid covariance matrix')
  do j=1,NVE
    read(iu,*) Bcov(:,j)
  enddo
  close(iu)

  ! ... Get the inverse of the error covariance matrix
  ! ...
  call matrinv(NVE,Bcov,Binv)

  ! ... Read the background field
  ! ...
  if (len_trim(mfile).eq.0) call stop_error(1,'No background field defined')

  write(*,*) 
  write(*,*) 'Reading background field ', trim(mfile)
  iu = unitfree()
  open(iu,file=mfile,status='old')
  read(iu,*) time, xb
  close(iu)
  write(*,*) 'xb = ', sngl(xb)
else
  alpha = zero
  Binv(:,:) = zero
endif



! ... Read the first guess
! ...
write(*,*) 
write(*,*) 'Reading first guess ', trim(ifile)

iu = unitfree()
open(iu,file=ifile,status='old')
read(iu,*) time, xo
close(iu)

write(*,*) 'xo = ', sngl(xo)



! ... Fill the step vector:
! ...
nn = nint((obs%time(obs%Nobs)-obs%time(1))/lorenz_dt)
call obs_map (obs,to,lorenz_dt,nn)
!do i=1,obs%Nobs
!  print*, i, obs%time(i), obs%step(i)
!enddo


! ... 
verbose = .false.
cost_to = to
cdt = (tf-to)/ncycles
write(*,*) 'to = ', to
write(*,*) 'tf = ', tf
write(*,*) 'dt = ', cdt

do cycle=1,ncycles

  write(*,*)
  write(*,*) 'Cycle = ', cycle

  cost_to = to
  cost_tf = to + cycle*cdt
  write(*,*) 'to            = ', cost_to
  write(*,*) 'tf            = ', cost_tf

  Jo = Jcost(xo)
  call Dcost(xo,gradJ)
  write(*,*) 'Initial state = ', sngl(xo)
  write(*,*) 'Inicial cost  = ', Jo
  write(*,*) 'Inicial grad  = ', sngl(gradJ)

  gtol = 1.0d-5
  x(:) = xo(:)
  call DFPMIN(x,lorenz_n,gtol,iter,Jmin,Jcost,DCost)

  write(*,*) 'optimal x = ', x
  write(*,*) 'optimal J = ', Jmin
  write(*,*) 'Num. iter = ', iter
  xo = x(:)
enddo

iu = unitfree()
write(*,*)
write(*,*) 'Writing minimization solution: ', trim(afile)
open(iu,file=trim(afile),status='unknown')
write(iu,'(4F10.4)') to, xo
close(iu)


! ... Opening output trajectory:
! ...
iu = unitfree()
write(*,*)
write(*,*) 'Opening output trajectory: ', trim(tfile)
open(iu,file=trim(tfile),status='unknown')
write(iu,'(4F10.4)') to, xo

! =============================================================
! ... Run the model
! ...

nsteps = (tf-to)/lorenz_dt

time = to
x(:) = xo(:)
do step=1,nsteps
  call rk4(lorenz_n,x,time,lorenz_dt,xn,lorenz_rhs)
  time = time + lorenz_dt
  x(:) = xn(:)
  write(iu,'(4F10.4)') time, x
enddo
close(iu)
! =============================================================

write(*,*) 'Opening file for final state: ', trim(efile)
open(iu,file=trim(efile),status='unknown')
write(iu,*) time, x
close(iu)

call stop_error(0,'FDVAR Ok')


end program fdvar
