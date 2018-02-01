! **************************************************************************
! ... nudging.f90
! ... Quim Ballabrera, February 2013
! ...
! ... System Dynamics
! ...
! **************************************************************************

PROGRAM nudging

implicit none

use sysdyn

!INCLUDE 'lorenz.h90'
!INCLUDE 'obs.h90'
!INCLUDE 'nudging.h90'

integer, parameter                          :: NVE=3

logical default,Zexist,continuous

integer i,ii,step,ensemble,nensembles,iu,idum
real(dp) Tf,T,Jc,a,b,mumax, mumin,Tmin
real(dp), DIMENSION(NVE)                :: Xo,X,Xn,Xm,Xs
real(dp), DIMENSION(:,:), ALLOCATABLE   :: Yint
character(LEN=180) ofile,dfile,opath,efile,nfile

REAL(KIND=8) gasdev,xran1

NAMELIST/namnudging/continuous,default,idum,to,dt,tf,sigma,rho,beta,ofile,Iobs,dfile,opath

EXTERNAL rhnud,rhwnud
EXTERNAL gasdev,xran1


idum = -1
nensembles = 20
Tx = 0.25
Ty = 0.25
Tz = 0.25

call getarg(1,nfile)

inquire(file=nfile,exist=Zexist)

if (Zexist) then
  write(*,*) 'Reading namelist file: ', trim(nfile)
  open(10,file=nfile,status='old')
  rewind(10)
  read(10,namnudging)
  close(10)
else
  default = .true.
endif

IF (default) THEN
  opath = './'
  ofile = './obs_noise_02_sampl_10.dat'
  dfile = 'nudging_noise_02_sampl_010.dat'

  ! ... Model parameters:
  ! ...
  sigma = 10.0D0
  rho   = 28.0D0
  beta  = 8.0D0/3.0D0

  ! ... Initial conditions
  ! ...
  to =  0.0D0
  tf = 20.0D0
  dt = 0.01D0

  ! ... What do we assimilate
  ! ...
  Iobs(:) = 1   ! All obs are assimilated

  ! ... How do assimilate:
  ! ...
  continuous = .false.
ENDIF

! ... Outputh directory:
! ...
! ... Open output path
! ...
IF ((TRIM(opath).EQ.'.').OR.         &
    (TRIM(opath).EQ.'./').OR.        &
    (LEN_TRIM(opath).EQ.0)) THEN
   ! ... Do nothing
   opath = '.'
ELSE
  i = LEN_TRIM(opath)
  IF (opath(i:i).NE.'/') opath(i+1:i+1) = '/'
  WRITE(*,*) 'mkdir -p '//TRIM(opath)
  CALL SYSTEM ('mkdir -p '//TRIM(opath))
ENDIF


! ... Model mean and standard deviation
! ...
OPEN(10,file='../makeobs/true_simulation_ave.dat',status='old')
REWIND(10)
READ(10,*) T, Xm
CLOSE(10)

OPEN(10,file='../makeobs/true_simulation_std.dat',status='old')
REWIND(10)
READ(10,*) T, Xs
CLOSE(10)

WRITE(*,*) 
WRITE(*,*) 'Long term average'
WRITE(*,'(3F9.3)') Xm
WRITE(*,*) 'Standard deviation'
WRITE(*,'(3F9.3)') Xs

WRITE(*,*)
WRITE(*,*) 'Output diagnostic file : ', TRIM(opath)//TRIM(dfile)

OPEN(72,file=TRIM(opath)//TRIM(dfile),status='unknown')
REWIND(72)

! ... Read the observations
! ...
CALL read_obs(ofile)

WRITE(*,*) 
WRITE(*,*) 'obs file          = ', TRIM(ofile)
WRITE(*,*) 'Nobs              = ', Nobs
WRITE(*,*) 'First time record = ', Tobs(1)
WRITE(*,*) 'Last time record  = ', Tobs(Nobs)
WRITE(*,*) 'Iobs              = ', Iobs(:)
WRITE(*,*) 'continuous        = ', continuous

WRITE(72,*) 
WRITE(72,*) 'obs file          = ', TRIM(ofile)
WRITE(72,*) 'Nobs              = ', Nobs
WRITE(72,*) 'First time record = ', Tobs(1)
WRITE(72,*) 'Last time record  = ', Tobs(Nobs)
WRITE(72,*) 'Iobs              = ', Iobs(:)
WRITE(72,*) 'continuous        = ', continuous

To = Tobs(1)
Tf = Tobs(Nobs)

nstep = NINT((tf-to)/dt)

WRITE(*,*) 
WRITE(*,*) 'Model initial time = ', to
WRITE(*,*) 'Model final time   = ', tf
WRITE(*,*) 'Model time step    = ', dt
WRITE(*,*) 'Number time steps  = ', nstep + 1

WRITE(72,*) 
WRITE(72,*) 'Model initial time = ', to
WRITE(72,*) 'Model final time   = ', tf
WRITE(72,*) 'Model time step    = ', dt
WRITE(72,*) 'Number time steps  = ', nstep + 1

! ... Make the mapping between observations and the model time steps and the 
! ...
ii = 1
t  = to
DO i=1,nstep+1
  IF (ii.LE.Nobs) THEN
  IF (ABS(t-Tobs(ii)).LT.0.49*dt) THEN
    Pobs(ii) = i
    ii = ii + 1
  ENDIF
  ENDIF
  t = t + dt
ENDDO

! ... Interpolate the obs to each time step of the model
! ...
ALLOCATE (Yint(NVE,nstep+1))

t = to
DO step=1,nstep+1
  IF (step.EQ.1) THEN
    Yint(:,1) = Yobs(:,1)
  ELSE IF (step.EQ.nstep+1) THEN
    Yint(:,nstep+1) = Yobs(:,Nobs)
  ELSE
    DO ii=1,Nobs-1
      IF (Tobs(ii).LE.t.AND.Tobs(ii+1).GT.t) EXIT
    ENDDO
! ... Y = a*t + b
! ... Yo(ii) = a*Tobs(ii) + b
! ... Yo(ii+1) = a*Tobs(ii+1) + b
! ... Yo(ii+1)-Yo(ii) = a*(Tobs(ii+1)-Tobs(ii))
    DO i=1,NVE
      a = (Yobs(i,ii+1)-Yobs(i,ii))/(Tobs(ii+1)-Tobs(ii))
      b = Yobs(i,ii)-a*Tobs(ii)
      Yint(i,step) = a*t + b
    ENDDO
  ENDIF
  t = t + dt
ENDDO

WRITE(*,*) 
WRITE(*,*) 'Idum         = ', idum
WRITE(*,*) 'nensembles   = ', nensembles

WRITE(72,*) 
WRITE(72,*) 'Idum         = ', idum
WRITE(72,*) 'nensembles   = ', nensembles

mumax = 1.0D0/dt
mumin = 1.0D0/(3.0d0*dt)

WRITE(*,*) 
WRITE(*,*) 'min mu = ', mumin, mumin*dt
WRITE(*,*) 'max mu = ', mumax, mumax*dt

WRITE(72,*) 
WRITE(72,*) 'min mu = ', mumin, mumin*dt
WRITE(72,*) 'max mu = ', mumax, mumax*dt

OPEN(45,file=TRIM(opath)//'initial.condition',status='unknown')
REWIND(45)
WRITE(45,*) NVE, nensembles

OPEN(46,file=TRIM(opath)//'nudging.coefs',status='unknown')
REWIND(46)

OPEN(47,file=TRIM(opath)//'scale.coefs',status='unknown')
REWIND(47)

Tmin = 1*dt

DO ensemble=1,nensembles

  DO i=1,NVE
    Xo(i) = Xm(i) + Xs(i)*gasdev(idum)
  ENDDO
  WRITE(45,*) ensemble*1.0
  WRITE(45,*) Xo(:)

  mux = (mumax-mumin)*xran1(idum) + mumin
  muy = (mumax-mumin)*xran1(idum) + mumin
  muz = (mumax-mumin)*xran1(idum) + mumin

  Tx = (0.25-Tmin)*xran1(idum) + Tmin
  Ty = (0.25-Tmin)*xran1(idum) + Tmin
  Tz = (0.25-Tmin)*xran1(idum) + Tmin

  Tx = 0.10D0; Ty = 0.10D0; Tz = 0.10D0

  WRITE(72,*)
  WRITE(72,*) 'ensemble, unit         = ', ensemble, 100+ensemble-1
  WRITE(72,*) 'Initial conditions     = ', Xo
  WRITE(72,*) 'X nudging coefficient  = ', mux
  WRITE(72,*) 'Y nudging coefficient  = ', muy
  WRITE(72,*) 'Z nudging coefficient  = ', muz
  WRITE(72,*) 'X decorrelation scale  = ', Tx
  WRITE(72,*) 'Y decorrelation scale  = ', Ty
  WRITE(72,*) 'Z decorrelation scale  = ', Tz
  WRITE(46,'(3F9.3)') mux,muy,muz
  WRITE(47,'(3F7.3)') Tx,Ty,Tz

! ... Nudging data assimilation:
! ...
  iu = 100 + ensemble - 1
  WRITE(efile,'(T1,"fort.",i3.3)') iu
  OPEN(iu,file=TRIM(opath)//TRIM(efile),STATUS='unknown')
  REWIND(iu)

  t = To
  X = Xo
  WRITE(iu,'(4F10.4)') t, X

  ii = 2
  DO step=1,nstep

    applyx = 0
    applyy = 0
    applyz = 0
    Yo(:) = Yint(:,step+1)

    IF (continuous) THEN
      applyx = 1
      applyy = 1
      applyz = 1
    ELSE
      IF (step.EQ.Pobs(ii)) THEN
        applyx = 1
        applyy = 1
        applyz = 1
        ii = ii + 1
      ENDIF
    ENDIF

    CALL rk4 (X,NVE,T,dt,XN,rhwnud)
    X = XN
    t = t + dt
    WRITE(iu,'(4F10.4)') t, X

  ENDDO
  CLOSE(iu)

ENDDO

CLOSE(45)

STOP 'Ok'
END PROGRAM nudging
! ...
! ==========================================================================
! ...
SUBROUTINE rhnud (t,x,dxdt)

IMPLICIT NONE

INCLUDE 'lorenz.h90'
INCLUDE 'nudging.h90'
INCLUDE 'obs.h90'

REAL(KIND=8), INTENT(in)    :: t
REAL(KIND=8), INTENT(in)    :: x(NVE)
REAL(KIND=8), INTENT(out)   :: dxdt(NVE)

dxdt(1) = sigma*(X(2)-X(1))  
dxdt(2) = (rho-X(3))*X(1) - X(2)
dxdt(3) = X(1)*X(2) - beta*X(3)

IF (Iobs(1).GE.1) dxdt(1) = dxdt(1) + applyx*mux*(Yo(1)-X(1))
IF (Iobs(2).GE.1) dxdt(2) = dxdt(2) + applyy*muy*(Yo(2)-X(2))
IF (Iobs(3).GE.1) dxdt(3) = dxdt(3) + applyz*muz*(Yo(3)-X(3))

RETURN
END SUBROUTINE rhnud
! ...
! ==========================================================================
! ...
SUBROUTINE rhwnud (t,x,dxdt)

IMPLICIT NONE

INCLUDE 'lorenz.h90'
INCLUDE 'nudging.h90'
INCLUDE 'obs.h90'

REAL(KIND=8), INTENT(in)    :: t
REAL(KIND=8), INTENT(in)    :: x(NVE)
REAL(KIND=8), INTENT(out)   :: dxdt(NVE)

REAL(KIND=8) xsum,ysum,zsum

CALL weight_obs (t,x,xsum,ysum,zsum)

dxdt(1) = sigma*(X(2)-X(1))
dxdt(2) = (rho-X(3))*X(1) - X(2)
dxdt(3) = X(1)*X(2) - beta*X(3)

IF (Iobs(1).GE.1) dxdt(1) = dxdt(1) + mux*xsum
IF (Iobs(2).GE.1) dxdt(2) = dxdt(2) + muy*ysum
IF (Iobs(3).GE.1) dxdt(3) = dxdt(3) + muz*zsum

RETURN
END SUBROUTINE rhwnud
! ...
! ==========================================================================
! ...
SUBROUTINE weight_obs (Time,X,xsum,ysum,zsum)

IMPLICIT NONE

INCLUDE 'lorenz.h90'
INCLUDE 'nudging.h90'
INCLUDE 'obs.h90'

REAL(KIND=8), INTENT(in)                   :: Time
REAL(KIND=8), DIMENSION(NVE), INTENT(in)   :: X
REAL(KIND=8), INTENT(out)                  :: xsum,ysum,zsum

INTEGER i
REAL(KIND=8) TT,argx,argy,argz,wx,wy,wz

xsum = 0.0D0
ysum = 0.0D0
zsum = 0.0D0
DO i=1,Nobs
  TT   = Time - Tobs(i)
  argx = TT/Tx
  argy = TT/Ty
  argz = TT/Tz
  argx = argx*argx
  argy = argy*argy
  argz = argz*argz
  wx   = exp(-argx)
  wy   = exp(-argy)
  wz   = exp(-argz)
  xsum = xsum + wx*(Yobs(1,i)-X(1))
  ysum = ysum + wy*(Yobs(2,i)-X(2))
  zsum = zsum + wz*(Yobs(3,i)-X(3))
ENDDO

RETURN
END SUBROUTINE weight_obs
