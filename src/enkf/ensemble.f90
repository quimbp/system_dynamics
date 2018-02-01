! **************************************************************************
! ... ensemble.f90
! ... Quim Ballabrera, January 2018
! ...
! ... System Dynamics
! ... Initial Ensemble. From the various alternatives, we use the 
! ... simplest one: Random sampling of the true trajectory of the 
! ... lorenz simulation. Other approaches are:
! ...   Random values
! ...   Empirical Orthogonal Functions
! ...   Singular vectors
! ...   Breeding vectors
! ... Some authors have shown that initialization methods matter little
! ... in some applications: Manusson et al. (2009)
! ... Thus, the intialization here is done as follows:
! ... (1) We read the trajectory of the Lorenz model,
! ...     we randomly take a series of random states separated more than
! ...     0.25 time units (25 time steps), which is the shortest decorrelation
! ...     scales of the model.
! ... (2) Random distribution around the long term mean, where the amplitude
! ...     of the perturbation is given by the standard deviation.
! ...
! **************************************************************************

program ensemble

use sysdyn

implicit none

integer, parameter                          :: NVE = lorenz_n

logical valid,Zexist,default

integer                                     :: idum     = -1
integer                                     :: rank     = 10
integer                                     :: method   = 1 ! The only coded one
real(dp)                                    :: rseed    = zero

integer i,j,k,kk,iu,nt,ko,kf
integer, dimension(:), allocatable          :: ind,indx

real(dp) tt
real(dp), dimension(NVE)                    :: xm,xs,Xbm
real(dp), dimension(:,:), allocatable       :: Xb,cov
character(LEN=180) nfile,tfile,efile,sfile,mfile

namelist/namens/rank,       &
                method,     &
                efile,      &
                tfile,      &
                mfile,      &
                sfile,      &
                rseed,      &
                rank,       &
                idum

call header()

! ... Default initialization
! ...
idum  = -1


nfile = ''
call getarg(1,nfile)
inquire(file=nfile,exist=Zexist)

if (len_trim(nfile).gt.0) then
  if (.not.Zexist) stop 'Namelist file not found'
endif

if (Zexist) then
  write(*,*) 'Reading ', trim(nfile)
  open(10,file=nfile,status='old')
  read(10,namens)
  close(10)
  default = .false.
else
  default = .true.
endif

if (default) then
  write(*,*) 'Default parameters'

  method = 1

  tfile = 'lorenz_simulation.dat'
  sfile = ''
  mfile = ''
  efile = 'Eo.matrix'

  rank   = 10
  rseed  = 0.

  ! ... If file enkf.namelist does not exist, we
  ! ... create a copy that can be later modified ....
  ! ...
  inquire(file='ensemble.namelist',exist=Zexist)
  if (.not.Zexist) then
    open(20,file='ensemble.namelist',status='new')
    write(20,nml=namens)
    close(20)
  endif

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


allocate (Xb(NVE,rank))

! ... Switch between initialization methods
! ... In the first version, only method 1 has been coded.
! ...
select case (method)

  case (1)
    iu = unitfree()
    write(*,*)
    write(*,*) 'Method 1'
    write(*,*) 'Sampling the trajectory of the model'
    write(*,*) 'Opening file ', trim(tfile)
    open(iu,file=tfile,status='old')
    nt = numlines(iu)
  
    allocate (ind(rank))
    allocate (indx(rank))
    ! ... Unconstrained sample
    ! ...
    ! ind = random_sample(rank,1,nt)
  
    ! ... Constrained sample
    ! ... To avoid taking two states too close:
    ! ...
    ind(1:1) = random_sample(1,1,nt)
    do k=2,rank
      valid = .false.
      do while (.not.valid)
        ind(k:k) = random_sample(1,1,nt)
        valid = .true.
        do kk=1,k-1
          if (abs(ind(kk)-ind(k)).le.25) valid = .false.
        enddo
      enddo
    enddo
 
    ! ... Sort the index: the order of the elements is in indx:
    ! ... ind(indx(1)) < ind(indx(2)) < ...
    ! ...
    call indexx(dble(ind),indx)
    write(*,'(T2,"Sampled states : ",1000(I6))') (ind(indx(i)),i=1,rank)

    rewind(iu)
    ko = 1
    do i=1,rank
      kf = ind(indx(i))
      do k=ko,kf
        read(iu,*) tt, Xb(:,i)
      enddo
      ko = kf + 1
    enddo
    close(iu)

  case (2)
    iu = unitfree()
    write(*,*)
    write(*,*) 'Method 2'
    write(*,*) 'Random ensemble around the long term mean'
    write(*,*) 'Opening mean file ', trim(mfile)
    open(iu,file=mfile,status='old')
    read(iu,*,err=2000,end=2000) tt, xm(:)
    close(iu)

    write(*,*) 'Opening std file ', trim(sfile)
    open(iu,file=sfile,status='old')
    read(iu,*,err=2000,end=2000) tt, xs(:)
    close(iu)

    do j=1,rank
      do i=1,NVE
        Xb(i,j) = xm(i) + xs(i)*randn()
      enddo
    enddo


  case default
    call stop_error(1,'Unknown ensemble initialization mode')

end select


! ... Save the ensemble
! ...
write(*,*) 
write(*,*) 'Saving ensemble in file ', trim(efile)
write(*,*) 'Ensemble'

iu = unitfree()
open(iu,file=efile,status='unknown')
rewind(iu)
write(iu,*) NVE, rank
do i=1,rank
  write(*,'(1000F9.4)') Xb(:,i)
  write(iu,*) Xb(:,i)
enddo
close(iu)


!! ... Reading normalization
!! ...
!write(*,*) 
!write(*,*) 'Reading standard deviation file :', trim(sfile)
!open(10,file=sfile,status='old')
!read(10,*) tt, Xs
!close(10)
!
! ... Reading initial ensemble:
! ...
!

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

! ... Ensemble covariance
! ...
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

call stop_error (0,'Ensemble Ok')


! ... Faulted run
! ...
2000 call stop_error (1,'Error reading input file')

end program ensemble
! ...
! ==========================================================================
! ...
