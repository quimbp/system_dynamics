! **************************************************************************
! ... test_adj.f90
! ... Quim Ballabrera, January 2018
! ...
! ... System Dynamics
! ... Adjoint methods applied to the Lorenz model
! ...
! ... Use the adjoint operator definition
! ...                         < A x, p > = < x, A*T p>
! ... where
! ... A = Linear Tangent Model
! ... x = Perturbation
! ... p = Adjoint variable
! ...
! **************************************************************************

program test_adj

use sysdyn

implicit none

integer, parameter                            :: NVE = lorenz_n

logical Zexist

integer iu,iter,rr,record
real(dp) alpha,tr,prod1,prod2
real(dp), dimension(NVE)                      :: xr,dx,Ldx
real(dp), dimension(NVE)                      :: p        ! Adjoint variable

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


write(*,*)
write(*,*) 'Adjoint validation:'
write(*,*) '------------------------------------------------------'

alpha = one
do iter=1,5

  ! ... The reference field is in vector LORENZ_REF(:)
  ! ...
  call lorenz_ltm(NVE,tr,alpha*dx,Ldx)
  prod1 = DOT_PRODUCT(Ldx,Ldx)

  p(:) = zero
  call lorenz_adj(NVE,tr,p,Ldx)
  prod2 = alpha*DOT_PRODUCT(p,dx)

  write(*,*)
  write(*,'(T2,"alpha, a = ",F9.6)') alpha
  write(*,'(T2,"< LTM(a*dx), LTM(a*dx) >  = ",F10.6)') prod1
  write(*,'(T2,"< ADJ(LTM(dx)), a*dx >    = ",F10.6)') prod2

  alpha = 0.1D0*alpha
enddo


! ... Success run
! ...
call stop_error(0,'TEST_ADJ Ok')

! ... Faulted run
! ...
2000 call stop_error(1,'TEST_ADJ: Error reading input file')

end program test_adj

