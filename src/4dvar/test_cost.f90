program test_cost

use sysdyn
use adjoint

implicit none

real(dp), parameter                                  :: epsilon = 1.0D-4

logical Zexist,default
integer iu,iter,step,nsteps,cycle,ncycles,i,j
integer, dimension(NVE)                              :: Iobs
real(dp), dimension(NVE)                             :: xo,x,xn,yo
real(dp), dimension(NVE)                             :: gradJ,gradN
real(dp) to,tf,time,Jo,Jmin,gtol,Jp
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
  read(iu,*) to, xb
  close(iu)
  write(*,*) 'xb = ', sngl(xb)
else
  alpha     = zero
  Binv(:,:) = zero
  xb(:)     = zero
endif



! ... Read the first guess
! ...
write(*,*) 
write(*,*) 'Reading first guess ', trim(ifile)

iu = unitfree()
open(iu,file=ifile,status='old')
read(iu,*) to, xo
close(iu)

write(*,*) 'xo = ', sngl(xo)


yo = obs%y(:,1)
write(*,*) 'yo = ', yo

! ... Fill the step vector:
! ...
nsteps = nint((obs%time(obs%Nobs)-obs%time(1))/lorenz_dt)

call obs_map (obs,to,lorenz_dt,nsteps)
!do i=1,obs%Nobs
!  print*, i, obs%time(i), obs%step(i)
!enddo

cost_to = to
cost_tf = tf

write(*,*) 'Reference cost and gradient'
Jo = Jcost(xo)
call Dcost(xo,gradJ)


write(*,*)
write(*,*) 'Numerical gradient'
verbose = .false.
! ...
Jp = Jcost(xo+epsilon*[one,zero,zero])
gradN(1) = (Jp-Jo)/epsilon

Jp = Jcost(xo+epsilon*[zero,one,zero])
gradN(2) = (Jp-Jo)/epsilon

Jp = Jcost(xo+epsilon*[zero,zero,one])
gradN(3) = (Jp-Jo)/epsilon

write(*,*) 'tf              = ', cost_tf
write(*,*) 'Initial state   = ', sngl(xo)
write(*,*) 'Referen state   = ', sngl(xb)
write(*,*) 'alpha           = ', alpha
write(*,*) 'Initial cost    = ', Jo
write(*,*) 'Initial grad    = ', sngl(gradJ)
write(*,*) 'Numerical grad  = ', sngl(gradN)



call stop_error(0,'TEST_COST Ok')
end program test_cost
