! **************************************************************************
! ... enkf.f90
! ... Quim Ballabrera, February 2013
! ...
! ... System Dynamics
! ... Ensemble Kalman Filter applied to the Lorenz model
! ...
! **************************************************************************

program enkf

use sysdyn

implicit none

integer, parameter                          :: NVE = lorenz_n
type(obs_type)                              :: obs

logical default,Zexist,noassim

! ... Lorenz model
! ...
integer                                     :: nstep
real(dp)                                    :: sigma,rho,beta
real(dp)                                    :: dt
real(dp)                                    :: to
real(dp)                                    :: tf

! ... Obs
! ...
integer, dimension(NVE)                     :: Iobs

! ... ENKF
! ...
integer                                     :: rank,Np,Ncycles,cycle
real(dp)                                    :: forget,rseed,r2x,r2y,r2z
real(dp), dimension(1)                      :: noise
real(dp), dimension(NVE)                    :: Xbm,Xam,XR
real(dp), dimension(:), ALLOCATABLE         :: Ro,Yo,Hxa
real(dp), dimension(:,:), ALLOCATABLE       :: Xb,Xa

integer i,j,k,jj,ll,step,kobs,idum,iu
real(dp) te,xJ,tt
real(dp), dimension(NVE)                    :: Xn,Xs,Xe
real(dp), dimension(:,:), ALLOCATABLE       :: cov
character(LEN=180) nfile,ofile,efile,sfile,opath,dfile,word
character(LEN=180) bfile,afile,yfile

namelist/namenkf/rank,         &
                 efile,        &
                 sfile,        &
                 ofile,        &
                 Iobs,         &
                 opath,        &
                 dfile,        &
                 bfile,        &
                 afile,        &
                 yfile,        &
                 r2x,          &
                 r2y,          &
                 r2z,          &
                 forget,       &
                 rseed,        &
                 sigma,        &
                 rho,          &
                 beta,         &
                 dt,           &
                 idum


call header()

! ... Default initialization
! ...
idum  = -1
sigma = lorenz_sigma
rho   = lorenz_rho
beta  = lorenz_beta
dt    = lorenz_dt
to    = lorenz_to
tf    = lorenz_tf


nfile = ''
call getarg(1,nfile)
inquire(file=nfile,exist=Zexist)

if (len_trim(nfile).gt.0) then
  if (.not.Zexist) stop 'Namelist file not found'
endif

if (Zexist) then
  write(*,*) 'Reading ', trim(nfile)
  iu = unitfree()
  open(iu,file=nfile,status='old')
  read(iu,namenkf)
  close(iu)
  default = .false.
else
  default = .true.
endif

if (default) then
  write(*,*) 'Default parameters'

  opath = './'
  ofile = 'obs_noise_02_sampl_010.dat'
  sfile = 'lorenz_simulation_std.dat'
  efile = 'Eo.matrix'
  dfile = 'enkf.diags'
  bfile = 'xb.dat'
  afile = 'xa.dat'
  yfile = 'innov.dat'

  ! ... Enkf main parameters
  ! ...
  rank   = 10
  forget = 1.2000
  rseed  = 0.10D0
  r2x    = 2.0
  r2y    = 2.0
  r2z    = 2.0

  ! ... Model parameters:
  ! ...
  sigma = lorenz_sigma
  rho   = lorenz_rho
  beta  = lorenz_beta
  dt    = lorenz_dt

  ! ... What do we assimilate
  ! ...
  Iobs(:) = 1   ! All obs are assimilated

  ! ... If file enkf.namelist does not exist, we
  ! ... create a copy that can be later modified ....
  ! ...
  inquire(file='enkf.namelist',exist=Zexist)
  if (.not.Zexist) then
    open(20,file='enkf.namelist',status='new')
    write(20,nml=namenkf)
    close(20)
  endif

endif

! ... Update, if necessary, the model parameters:
! ...
lorenz_sigma = sigma
lorenz_rho   = rho
lorenz_beta  = beta


! ... Create output path
! ...
if ((trim(opath).EQ.'.').OR.         &
    (trim(opath).EQ.'./').OR.        &
    (len_trim(opath).EQ.0)) then
   ! ... Do nothing
   opath = '.'
else
  call SYSTEM ('mkdir -p '//trim(opath))
ENDIf

! ... Output Diagnostic file
! ...
open(72,file=filename_append(opath,dfile),STATUS='unknown')


! ... Read the observations
! ...
obs = obs_read(ofile)
Ncycles = obs%Nobs - 1         ! The first observation is not assimilated

if (SUM(Iobs).EQ.0) then
  noassim = .true.
else
  noassim = .false.
endif


! ... Random number generator seed
! ... if idum = -1, the default seed is used
! ...
if (idum.eq.-1) then
  write(*,*) 'Default random number generator seed'
else
  write(*,*) 'Resetting random number generator seed'
  call randseed(idum)
endif


! ... Observation error variance
! ...
XR = (/r2x,r2y,r2z/)

Np = SUM(Iobs)
allocate (Ro(Np))
allocate (Yo(Np))
allocate (Hxa(Np))
call Hforward (XR,NVE,Iobs,Np,Ro)

write(72,*) 
write(72,*) 'First time record = ', obs%time(1)
write(72,*) 'Last time record  = ', obs%time(obs%Nobs)
write(72,*) 'First time interv = ', obs%time(2)-obs%time(1)
write(72,*) 'Iobs  = ', Iobs(:)
write(72,*) 'Np    = ', Np

to = obs%time(1)
tf = obs%time(obs%Nobs)

nstep = NINT((tf-to)/dt)
write(72,*) 
write(72,*) 'Model dimension    = ', NVE
write(72,*) 'Model initial time = ', to
write(72,*) 'Model final time   = ', tf
write(72,*) 'Model time step    = ', dt
write(72,*) 'Number time steps  = ', nstep + 1
write(72,*) 'Num Assim cycles   = ', Ncycles


! ... Reading normalization
! ...
write(*,*) 
write(*,*) 'Reading standard deviation file :', trim(sfile)
iu = unitfree()
open(iu,file=sfile,status='old')
read(iu,*) tt, Xs
close(iu)

! ... Reading initial ensemble:
! ...
allocate (Xb(NVE,rank))

write(*,*) 'Reading initial ensemble: ', trim(efile)
iu = unitfree()
open(iu,file=efile,status='old')
read(iu,*) jj,ll
if (jj.ne.NVE) call stop_error(1,'Incompatible dimension')
if (ll.LT.rank) call stop_error(1,'ERROR: too few modes')
do k=1,rank
  read(iu,*) Xb(:,k)
enddo
close(iu)


! ... Mean from ensemble.
! ...
call emean (Xb,NVE,rank,Xbm)
write(*,*) 
write(*,*) 'Ensemble mean      = ', SNGL(Xbm)

! ... Ensemble anomalies
! ...
do i=1,NVE
  Xb(i,:) = Xb(i,:) - Xbm(i)        ! Remove the mean
enddo


write(*,*) 
write(*,*) 'Ensemble covariance '
allocate (cov(NVE,NVE))
do j=1,NVE
do i=j,NVE
  cov(i,j) = DOT_PRODUCT(Xb(i,:),Xb(j,:))/(rank-1.0D0)
  cov(j,i) = cov(i,j)
enddo
enddo
do i=1,NVE
  write(*,'(100F9.3)') cov(i,:)
enddo


kobs = 1

! ... Analysis
! ...
allocate (Xa(NVE,rank))
do i=1,NVE
  Xa(i,:) = Xbm(i) + Xb(i,:)
enddo
Xam(:) = Xbm(:)


open(34,file=filename_append(opath,bfile),status='unknown')
open(35,file=filename_append(opath,afile),status='unknown')
open(36,file=filename_append(opath,yfile),status='unknown')
write(34,'(4F10.4)') to, Xbm
write(35,'(4F10.4)') to, Xam

iu = 100
do i=1,rank
  iu = iu + 1
  write(word,'(T1,"/fort.",I3.3)') iu
  open(iu,file=trim(opath)//TRIM(word),STATUS='unknown')
enddo

  
do cycle=1,ncycles

  kobs = kobs + 1
  tf = obs%time(kobs)
  nstep = NINT((tf-to)/dt)

  ! ... What is being observed: obs%y has all the components
  ! ... with Hforward we extract the selected observations
  ! ... as indicated by the values of Iobs:
  ! ...
  call Hforward (obs%y(1:NVE,kobs),NVE,Iobs,Np,Yo)

  write(*,*)
  write(*,*) 'Cycle ', cycle, ' / ', Ncycles
  write(*,*)

  write(72,*) 'Cycle ', cycle, ' / ', Ncycles
  write(72,*) 'Obs pointer   : ', kobs
  write(72,*) 'to            : ', SNGL(to)
  write(72,*) 'tf (kobs)     : ', SNGL(tf)
  write(72,*) 'Yobs          : ', SNGL(Yo)
  write(72,*) 'nstep         : ', nstep
  write(72,*) 'rseed         : ', rseed


  ! ... We start with a model simulation
  ! ... 
  iu = 100
  do ll=1,rank
    iu = iu + 1
    te = to
    do i=1,NVE
      noise = randn()
      Xe(i) = Xa(i,ll) + rseed*Xs(i)*noise(1)
    enddo
    write(iu,'(4F10.4)') Te, Xe
    do step=1,nstep
      call rk4 (lorenz_n,Xe,Te,dt,Xn,lorenz_rhs)
      Xe = Xn
      te = te + dt
      write(iu,'(4F10.4)') Te, Xe
    enddo
    if (ABS(te-tf).GT.0.1*dt) call stop_error (1,'Lost time track')
    Xb(:,ll) = Xe(:)
  enddo

  ! ... Analysis step
  ! ...
  if (noassim) then
    Xa(:,:)  = Xb(:,:)
    call emean (Xb,NVE,rank,Xbm)
    Xam(:) = Xbm(:)
  else
    call etkf (72,forget,Xb,Xa,xbm,xam,NVE,rank,Iobs,Yo,Ro,Np,xJ)
  endif
  print*, 'xbm = ', SNGL(xbm)
  print*, 'xam = ', SNGL(xam)


  write(34,'(4F10.4)') tf, Xbm
  write(35,'(4F10.4)') tf, Xam

  call Hforward(xbm,NVE,Iobs,Np,Hxa)
  write(36,'(4F10.4)') tf, Yo-Hxa

  to = obs%time(kobs)

enddo


call stop_error (0,'EnKF Ok')
end program enkf
! ...
! ==========================================================================
! ...
