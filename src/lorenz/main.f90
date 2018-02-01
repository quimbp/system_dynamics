use sysdyn

implicit none

logical                               :: fini = .false.
logical                               :: fpar = .false.
logical                               :: ftra = .false.
logical                               :: fend = .false.
logical                               :: plot = .false.

character(len=180)                    :: ifile = 'lorenz.ini'
character(len=180)                    :: pfile = 'lorenz.params'
character(len=180)                    :: tfile = 'lorenz.dat'
character(len=180)                    :: efile = 'lorenz.end'

logical                               :: fto = .false.
logical                               :: ftf = .false.
logical                               :: fdt = .false.

real(dp)                              :: to = zero
real(dp)                              :: tf = 20.0_dp
real(dp)                              :: dt = 0.01_dp

integer                               :: na,ierr,step,nsteps
real(dp), dimension(lorenz_n)         :: x,xo,xn
real(dp)                              :: time

namelist /PARAMS/ lorenz_sigma, lorenz_beta, lorenz_rho

call header()

! ... Lineargs
! ...
call lineargs_ini(na)

call argstr('-pa',fpar,pfile)
call argstr('-i',fini,ifile)
call argstr('-tr',ftra,tfile)
call argstr('-e',fend,efile)
call argdbl('-to',fto,to)
call argdbl('-tf',ftf,tf)
call argdbl('-te',ftf,tf)
call argdbl('-dt',fdt,dt)
call argflg('-pl',plot)

! ... Model parameters
! ... It tries to open the file lorenz.params
! ... If it fails, the parameters are the default ones from sysdyn
! ...
write(*,*)
open(1,file=pfile,status='old',iostat=ierr)
if (ierr.eq.0) then
  write(*,*) 'Reading lorenz model parameters from file: ', trim(pfile)
  rewind(1)
  read(1,nml=PARAMS)
  close(1)
else
  write(*,*) 'Using default lorenz parameters'
endif

write(*,*) 'sigma = ', lorenz_sigma
write(*,*) 'beta  = ', lorenz_beta 
write(*,*) 'rho   = ', lorenz_rho  
write(*,*)

! ... Model initial conditions
! ... It tries to open the file lorenz.ini
! ... If it fails, the parameters are the default ones from sysdyn
! ...
open(1,file=ifile,status='old',iostat=ierr)
if (ierr.eq.0) then
  write(*,*) 'Reading initial condition from file: ', trim(ifile)
  read(1,*) time, xo(:)
  if (.not.fto) to = time
else
  write(*,*) 'Using default initial condition'
  xo(:) = lorenz_ini(:)
endif
 
write(*,*) 'xo = ', xo 
write(*,*)

nsteps = abs((tf - to)/dt)
write(*,*) 
write(*,*) 'Initial time  = ', to
write(*,*) 'Final   time  = ', tf
write(*,*) 'Time interval = ', dt
write(*,*) 'Model steps   = ', nsteps

! ... Opening output trajectory:
! ...
write(*,*) 
write(*,*) 'Opening output trajectory: ', trim(tfile)
open(10,file=trim(tfile),status='unknown')
write(10,'(4F10.4)') to, xo

! =============================================================
! ... Run the model
! ...

time = to
x(:) = xo(:)
do step=1,nsteps
  call rk4(lorenz_n,x,time,dt,xn,lorenz_rhs)
  time = time + dt
  x(:) = xn(:)
  write(10,'(4F10.4)') time, x
enddo
close(10)
! =============================================================

write(*,*) 'Opening file for final state: ', trim(efile)
open(10,file=trim(efile),status='unknown')
write(10,*) time, x
close(10)

! ... Call python to plot the results if option -plot has been used
! ...
if (plot) call system('python3 lorenz.py '//trim(tfile))

call stop_error(0,'lorenz model Ok')
end

