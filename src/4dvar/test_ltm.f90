! **************************************************************************
! ... test_ltm.f90
! ... Quim Ballabrera, January 2018
! ...
! ... System Dynamics
! ... Adjoint methods applied to the Lorenz model
! ...
! ... Test the Tangent Linear Model (the derivative of the Linear Model)
! ... by calculating the numerical local derivative of the Nonlinear model
! ...
! ...                 L dx  = [M(xr+alpha*dx) -M(xr)]/alpha
! ... So,
! ...          [ (M(xr+alpha*dx)-M(xr) ] / [ alpha*L dx ] -> 1
! ...
! **************************************************************************

program test_ltm

use sysdyn

implicit none

integer, parameter                            :: NVE = lorenz_n

logical Zexist

integer iu,i,iter,rr,record
real(dp) alpha,tr
real(dp), dimension(NVE)                      :: xr,yr,dx,dy
real(dp), dimension(NVE)                      :: xp,yp,dyp

character(len=180) nfile,tfile,dfile

NAMELIST /namltm/ tfile,          &
                  record,         &
                  dfile

call header()

! ... Default initialization
! ...
tfile  = '../../work/lorenz_simulation.dat'
record = 101
dx     = [0.1D0,0.2D0,0.5D0]


call getarg(1,nfile)
inquire(file=nfile,exist=Zexist)

if ((len_trim(nfile).gt.0).and. &
    (.not.Zexist))    call stop_error(1,'Namelist file not found')


if (Zexist) then
  write(*,*) 'Reading ', trim(nfile)
  iu = unitfree()
  open(iu,file=nfile,status='old')
  rewind(iu)
  read(iu,namltm)
  close(iu)

  ! ... Read the perturbation:
  ! ...
  open(iu,file=tfile,status='old')
  read(iu,*) tr, dx
  close(iu)
  
endif


iu = unitfree()
open(iu,file=tfile,status='old')
do rr=1,record
  read(iu,*,end=2000,err=2000) tr, xr
enddo

! ... Save the reference for the LTM code
! ... LORENZ_REF(1:LORENZ_N) is a vector in the LORENZ module
! ...
lorenz_ref(:) = xr(:)


write(*,*) 
write(*,*) 'Reference state : '
write(*,'(100(F12.7))') SNGL(xr)

write(*,*) 
write(*,*) 'Perturbation vector : '
write(*,'(100(F12.8))') SNGL(dx)


call lorenz_rhs(NVE,tr,xr,yr)     ! Reference Run

write(*,*)
write(*,*) 'Positive perturbations:'
write(*,*) '------------------------------------------------------'

alpha = one
do iter=1,5
  xp(:) = xr(:) + alpha*dx(:)
  call lorenz_rhs(NVE,tr,xp,yp)
  dyp(:) = yp(:) - yr(:)

  call lorenz_ltm(NVE,tr,dx,dy)
  dy = alpha*dy

  write(*,*)
  write(*,'(T2,"alpha, a = ",F9.6)') alpha
  write(*,'(T2,"M(x+a*dx)-M(x)               = ",8F10.6)') dyp(:)
  write(*,'(T2,"a*LTM(dx)                    = ",8F10.6)') dy(:)
  write(*,'(T2,"(M(x+a*dx)-M(x))/(a*LTM(dx)) = ",8F10.6)') (dyp(i)/dy(i),i=1,NVE)

  alpha = 0.1D0*alpha
enddo


write(*,*)
write(*,*) 'Negative perturbations:'
write(*,*) '------------------------------------------------------'

alpha = one
do iter=1,5
  xp(:) = xr(:) - alpha*dx(:)
  call lorenz_rhs(NVE,tr,xp,yp)
  dyp(:) = yp(:) - yr(:)

  call lorenz_ltm(NVE,tr,dx,dy)
  dy = alpha*dy

  write(*,*)
  write(*,'(T2,"alpha, a = ",F9.6)') alpha
  write(*,'(T2,"M(x-a*dx)-M(x)               = ",8F10.6)') dyp(:)
  write(*,'(T2,"a*LTM(dx)                    = ",8F10.6)') dy(:)
  write(*,'(T2,"(M(x-a*dx)-M(x))/(a*LTM(dx)) = ",8F10.6)') (dyp(i)/dy(i),i=1,NVE)

  alpha = 0.1D0*alpha
enddo




! ... Success run
! ...
call stop_error(0,'TEST_LTM Ok')

! ... Faulted run
! ...
2000 call stop_error(1,'TEST _LTM: Error reading input file')

end program test_ltm

