use sysdyn

real(dp), parameter                         :: dt=0.01D0

integer i,no,n,NVE
real(dp), dimension(:), allocatable     :: to,t
real(dp), dimension(:,:), allocatable   :: xo,xa
character(len=180) line
character(len=180) ifile,ofile

call header()

ifile = ''
ofile = ''
call getarg(1,ifile)
if (len_trim(ifile).eq.0) call stop_error(1,'Empty input filename')

call getarg(2,ofile)
if (len_trim(ofile).eq.0) ofile = 'out.dat'

write(*,*) 'Reading file ', TRIM(ifile)
open(10,file=ifile,status='old')
no = numlines(10)
write(*,*) 'Number of lines : ', no

rewind(10)
read(10,'(A)') line
NVE = numwords(line) - 1
write(*,*) 'NVE = ', NVE


allocate (to(no))
allocate (xo(NVE,no))
rewind(10)
do i=1,no
  read(10,*) to(i), xo(:,i)
enddo
close(10)

n = NINT((to(no)-to(1))/dt + 1)

allocate (t(n))
allocate (xa(NVE,n))

t(1) = to(1)
do i=2,n
  t(i) = t(i-1) + dt
enddo

! ... Interpolation using Akima splines
! ... In librery libsysdin.a
! ...
do i=1,NVE
   xa(i,:) = akima(to,xo(i,:),t)
enddo

write(*,*) 'Saving interpolation in ', TRIM(ofile)
open(10,file=ofile,status='unknown')
rewind(10)
do i=1,n
  write(10,'(4F10.4)') t(i), xa(:,i)
enddo
close(10)


call stop_error(0,'interpolation Ok')
end

