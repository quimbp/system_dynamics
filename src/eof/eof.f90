! **************************************************************************
! ... eof.f90
! ... Quim Ballabrera, March 2013
! ...
! ... System Dynamics
! ...
! **************************************************************************

program eof

use sysdyn

implicit none

logical Zexist
integer i,j,k,Nv,Nt,step,rank,iu

real(dp) xsum,xsum2
real(dp), dimension(:), allocatable        :: Zm,Zs
real(dp), dimension(:,:), allocatable      :: Z,cov
real(dp), dimension(:), allocatable        :: Lambda,ev,pc
real(dp), dimension(:,:), allocatable      :: U
real(dp), dimension(:), allocatable        :: T

character(len=180) nfile,tfile,cfile,efile,pfile,mfile,sfile,xfile
character(len=180) rfile,line

namelist/nameof/tfile,     &
                cfile,     &
                efile,     &
                pfile,     &
                mfile,     &
                sfile,     &
                xfile

call header()

tfile = 'lorenz_simulation.dat'
mfile = 'lorenz.mean'
sfile = 'lorenz.stdv'
cfile = 'lorenz.cov'
rfile = 'lorenz.cor'
efile = 'lorenz.eof'
xfile = 'lorenz.spm'
pfile = 'lorenz.pc'

nfile = ''
CALL getarg(1,nfile)
inquire(file=nfile,exist=Zexist)

if (len_trim(nfile).gt.0) then
  if (Zexist) then
    write(*,*) 'Reading ', trim(nfile)
    iu = unitfree()
    open(iu,FILE=nfile,STATUS='old')
    rewind(iu)
    read(iu,nameof)
    close(iu)
  else
    call stop_error (1,'Namelist file not found')
  endif
endif


! ... Open input trajectory
! ...
iu = unitfree()
write(*,*) 'Opening file ', trim(tfile)
open(iu,file=tfile,status='old')
rewind(iu)
Nt = numlines(iu)

rewind(iu)
read(iu,'(A)') line
Nv = numwords(line) - 1

allocate (Z(Nv,Nt))
allocate (T(Nt))
allocate (Zs(Nv))
allocate (cov(Nv,Nv))

! ... Load the data:
! ...
rewind(iu)
do step=1,Nt
  READ(iu,*) T(step),Z(:,step)
enddo
close(iu)

! ... Long term mean:
! ...
allocate (Zm(Nv))
Zm(:) = 0.0D0
do step=1,Nt
  Zm(:) = Zm(:) + Z(:,step)
enddo
Zm(:) = Zm(:) / Nt

write(*,*) 'Saving long term mean ', trim(mfile)
open(iu,file=mfile,status='unknown')
rewind(iu)
write(iu,*) 0.0, SNGL(Zm)
close(iu)

! ... Anomalies:
! ...
do i=1,Nv
  Z(i,:) = Z(i,:) - Zm(i)
enddo

! ... Calculate the covariance matrix
! ...
do j=1,Nv
do i=j,Nv
  xsum = 0.0D0
  do k=1,Nt
    xsum = xsum + Z(i,k)*Z(j,k)
  enddo
  cov(i,j) = xsum / (Nt-1.0D0)
  cov(j,i) = cov(i,j)
enddo
enddo

write(*,*) 'Covariance : '
do i=1,Nv
  write(*,'(100F9.3)') cov(i,:)
enddo

write(*,*) 'Saving covariance matrix ', trim(cfile)
iu = unitfree()
open(iu,file=cfile,status='unknown')
rewind(iu)
write(iu,*) Nv, Nv
do j=1,Nv
  write(iu,*) cov(:,j)
enddo
close(iu)


! ... Standard deviation
! ...
do i=1,Nv
  Zs(i) = SQRT(DOT_PRODUCT(Z(i,:),Z(i,:))/(Nt-1.0D0))
enddo

write(*,*) 'Saving standard deviation ', trim(sfile)
open(iu,file=sfile,status='unknown')
rewind(iu)
write(iu,*) 0.0D0, SNGL(Zs)
close(iu)


! ... Scale the anomalies:
! ...
do i=1,Nv
  Z(i,:) = Z(i,:)/Zs(i)
enddo

! ... Calculate the correlation matrix
! ...
do j=1,Nv
do i=j,Nv
  xsum = 0.0D0
  do k=1,Nt
    xsum = xsum + Z(i,k)*Z(j,k)
  enddo
  cov(i,j) = xsum / (Nt-1.0D0)
  cov(j,i) = cov(i,j)
enddo
enddo

write(*,*)
write(*,*) 'Mean : '
write(*,*) SNGL(Zm)
write(*,*) 'Standard deviation : '
write(*,*) SNGL(Zs)
write(*,*) 'Correlation : '
do i=1,Nv
  write(*,'(100F9.3)') cov(i,:)
enddo

! ... Save the covariance matrix
! ...
write(*,*) 'Saving correlation matrix ', trim(rfile)
iu = unitfree()
open(iu,file=rfile,status='unknown')
rewind(iu)
write(iu,*) Nv, Nv
do j=1,Nv
  write(iu,*) cov(:,j)
enddo
close(iu)


! ... Just to be safe: print the trace of the covariance matrix:
! ...
xsum = 0.0D0
do i=1,Nv
do k=1,Nt
  xsum = xsum + Z(i,k)*Z(i,k)
enddo
enddo
write(*,*) 'Trace coavariance matrix (Z Z^T): ', xsum/(Nt-1.0D0)

! ... EOFs:
! ...
rank = MIN(nv,nt-1)
write(*,*) 'Rank: ', rank

allocate (U(Nv,rank))
allocate (Lambda(rank))

call get_the_eofs (Z,Nv,Nt,rank,U,Lambda)

! ... The covariance matrix seen by the EOFs: Cov = U Lambda U^T
! ...
xsum = 0.0D0
do i=1,Nv
do k=1,rank
  xsum = xsum + Lambda(k)*U(i,k)*U(i,k)
enddo
enddo
write(*,*) 'Trace coavariance matrix (U Lambda U^T): ', xsum

! ... Explained variance:
! ...
allocate (ev(rank))
xsum  = SUM(Lambda(:))
xsum2 = 0.0D0
do k=1,rank
  xsum2 = xsum2 + 100.0D0*Lambda(k)/xsum
  ev(k) = xsum2
enddo

write(*,*) 'Saving scpectrum file ', trim(xfile)
iu = unitfree()
open (iu,file=xfile,status='unknown')
rewind(iu)
do k=1,rank
  write(iu,'(i5,F18.5,2F9.3)') k, Lambda(k), 100.0*Lambda(k)/xsum, ev(k)
enddo
close(iu)

write(*,*) 'Saving the EOFs matrix ', trim(efile)
iu = unitfree()
open(iu,file=efile,status='unknown')
rewind(iu)
write(iu,*) Nv, rank
do k=1,rank
  write(iu,*) U(:,k)
enddo
close(iu)

! ... Principal components:
! ...
allocate (pc(rank))
write(*,*) 'Saving the PC matrix ', trim(pfile)
iu = unitfree()
open(iu,file=pfile,status='unknown')
rewind(iu)

do step=1,Nt
  do k=1,rank
    xsum = 0.0D0
    do i=1,Nv
      xsum = xsum + U(i,k)*Z(i,step)
    enddo
    pc(k) = xsum
  enddo
  write(iu,'(4F10.4)') T(step), pc(:)
enddo
close(iu)


call stop_error (0,'EOF Ok')
end program eof
! ...
! =====================================================================
! ... 
SUBROUTINE get_the_eofs (Z,Nv,Nt,rank,U,Lambda)

use sysdyn

implicit none

integer Nv,Nt,rank
real(dp), dimension(Nv,Nt), intent(in)       :: Z
real(dp), dimension(Nv,rank), intent(out)    :: U
real(dp), dimension(rank), intent(out)       :: Lambda

integer Nm,k
real(dp), dimension(:), allocatable          :: DD
real(dp), dimension(:,:), allocatable        :: UU,VV


if (Nv.GE.Nt) then
  write(*,*) 'Singular Value Decomposition of Z ...'

  Nm = Nt
  allocate (UU(Nv,Nt))
  allocate (VV(Nt,Nt))
  allocate (DD(Nt))

  UU(:,:) = Z(:,:)
  call svdcmp (UU,Nv,Nt,Nv,Nt,DD,VV)

else
  write(*,*) 'Singular Value Decomposition of Z^T ...'

  Nm = Nv
  allocate (UU(Nv,Nv))
  allocate (VV(Nt,Nv))
  allocate (DD(Nv))

  VV = TRANSPOSE(Z)
  call svdcmp (VV,Nt,Nv,Nt,Nv,DD,UU)

endif

DD(:) = DD(:)**2
call eigsort (DD,Nm,UU,Nv)

IF (rank.GT.Nm) call stop_error(1,'ERROR: rank to large !!!!!')

! ... Output arrays:
! ...
do k=1,rank
  Lambda(k) = DD(k)/(Nt-1.0D0)
  U(:,k)    = UU(:,k)
enddo

deallocate (UU)
deallocate (VV)
deallocate (DD)

return
end
! ...
! ====================================================================
! ...
