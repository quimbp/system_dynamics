! **************************************************************************
! ... makeobs.f90
! ... Quim Ballabrera, February 2013
! ...
! ... System dynamics
! ...
! **************************************************************************

program makeobs

use sysdyn

logical                               :: default = .false.
logical                               :: Zexist

real(dp), dimension(:), allocatable   :: Y,Xm,Xs,Xe,Xo
real(dp), dimension(:), allocatable   :: T
real(dp), dimension(:,:), allocatable :: X

integer nl,no,step,idum,n,iu
real(dp) dto,var,t1,t2,noise(1)
character(len=180) ifile,ofile,nfile,line

namelist /namobs/ idum,ifile,ofile,var,dto

call header()

nfile = ''
call getarg(1,nfile)
inquire(file=nfile,exist=Zexist)
if (Zexist) then
  write(*,*) 'Reading ', trim(nfile)
  open(10,file=nfile,status='old')
  rewind(10)
  read(10,namobs)
  close(10)
  default = .false.
else
  write(*,*) 'No valid namelist specified.'
  write(*,*) 'Default makeobs parameters'
  default = .true.
endif

if (default) then
  write(*,*) 'Using default parameters ...'
  idum  = -1
  ifile = './lorenz_simulation.dat'
  ofile = 'lorenz_obs.dat'
  var   = 2.00
  dto   = 0.25
endif

iu = unitfree()
write(*,*) 'Reading lorenz simulation from : ', TRIM(ifile)
open(iu,file=ifile,status='old')

nl = numlines(iu)
IF (nl.LE.1) STOP 'MAKEOBS ERROR: not enough model outputs'

rewind(iu)
read(iu,'(A)') line
n = numwords(line) - 1     ! First column is time

write(*,*)
write(*,*) 'n         = ', n
write(*,*) 'idum      = ', idum
write(*,*) 'ifile     = ', TRIM(ifile)
write(*,*) 'ofile     = ', TRIM(ofile)
write(*,*) 'var       = ', var
write(*,*) 'dto       = ', dto
write(*,*)

! ... Random number generator seed
! ... if idum = -1, the default seed is used
! ...
if (idum.eq.-1) then
  write(*,*) 'Default random number generator seed'
else
  write(*,*) 'Resetting random number generator seed'
  call randseed(idum)
endif

allocate(Xo(n))
allocate(Xm(n))
allocate(Xs(n))
allocate(Xe(n))
allocate(Y(n))

rewind(iu)
read(iu,*) t1,Xo
read(iu,*) t2,Xo
dt = t2 - t1
no = nint(dto/dt)

write(*,*) 
write(*,*) 'Model entries             : ', nl
write(*,*) 'Model time-step           : ', dt
write(*,*) 'Observation time-step     : ', dto
write(*,*) 'Observation sampling rate : ', no

! ... Trajectory sampling
! ...
allocate (T(nl))
allocatE (X(n,nl))

rewind(iu)
do i=1,nl
  read(iu,*,err=1000) T(i), X(:,i)
enddo

write(*,*) 'Initial time              : ', T(1)
write(*,*) 'Final time                : ', T(nl)

! ... Standard deviation
! ...
do i=1,n
  CALL avevar (X(i,:),Xm(i),Xs(i))
  Xs(i) = maxval((/Xs(i),0.0D0/))
  Xs(i) = SQRT(Xs(i))
enddo
write(*,*) 
write(*,*) 'Mean value and Standard deviation'
write(*,'(T2,A,T5,F9.3,T15,F9.3)') 'X',Xm(1),Xs(1)
write(*,'(T2,A,T5,F9.3,T15,F9.3)') 'Y',Xm(2),Xs(2)
write(*,'(T2,A,T5,F9.3,T15,F9.3)') 'Z',Xm(3),Xs(3)


open(9,file='lorenz_simulation_ave.dat',status='unknown')
write(9,'(4F10.4)') 0.0, Xm(:)
close(9)

open(9,file='lorenz_simulation_std.dat',status='unknown')
write(9,'(4F10.4)') 0.0, Xs(:)
close(9)


Xe = sqrt(var)

write(*,*) 
write(*,*) 'Error variance ', var
write(*,*) 'Standard deviation of error'
write(*,'(T2,A,T5,F9.3)') 'X',Xe(1)
write(*,'(T2,A,T5,F9.3)') 'Y',Xe(2)
write(*,'(T2,A,T5,F9.3)') 'Z',Xe(3)

open(20,file=ofile,status='unknown')
do step=1,nl
  if (MOD(step,no).EQ.1) THEN
    do i=1,n
      noise = randn()
      Y(i) = X(i,step) + Xe(i)*noise(1)
    enddo
    write(20,'(4F10.4)') T(step), Y(:)
  endif
enddo

call stop_error(0,'makeobs program Ok')

1000 call stop_error(1,'Error while reading the trajectory')

end program makeobs
