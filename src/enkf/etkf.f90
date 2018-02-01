! **************************************************************************
! ... etkf.f90
! ... Quim Ballabrera, February 2013
! ...
! ... System Dynamics:
! ... ETKF equations
! ... The code assumes diagonal observation error covariance matrix
! ...
! **************************************************************************


subroutine etkf (iu,forget,X,Xa,xf,xam,n,rank,iobs,yobs,R,Nobs,xJ)

use sysdyn

! ...    n       :: dimension of the system
! ...    rank    :: ensemble size
! ...    Nobs    :: number of observations

implicit none

integer, intent(in)                            :: iu,n,rank,Nobs
real(kind=8), intent(in)                       :: forget
real(kind=8), dimension(n,rank), intent(inout) :: X,Xa
real(kind=8), dimension(n), intent(out)        :: xf,xam
integer, dimension(Nobs), intent(in)           :: iobs
real(kind=8), dimension(Nobs), intent(in)      :: yobs,R
real(kind=8), intent(out)                      :: xJ

logical verb
integer i,j,k,ii
real(kind=8) xlami,xsum,xJ1,xJ2
real(kind=8), dimension(:), ALLOCATABLE        :: wm,ws
real(kind=8), dimension(Nobs)                  :: Hxf
real(kind=8), dimension(Nobs)                  :: innovation
real(kind=8), dimension(:), ALLOCATABLE        :: vwork,vwork2
real(kind=8), dimension(:,:), ALLOCATABLE      :: rwork,Wa
real(kind=8), dimension(Nobs,rank)             :: HE
real(kind=8), dimension(rank,Nobs)             :: HER
real(kind=8), dimension(:,:), ALLOCATABLE      :: aLambda,U

if (rank.LE.1) call stop_error(1,'ETKF ERROR: rank <= 1')

xlami = (DBLE(rank)-1.0D0)/forget

! ... First thing: evaluate HX, Center it, and calculate anomalies:
! ...
write(iu,*) 'etkf    '
write(iu,*) 'Using Nobs    : ', Nobs
write(iu,*) 'iobs(:)       : ', iobs(:)
write(iu,*) 'Using forget  : ', forget


! ................................................................
! ................................................................
! ... Global operations: Does not matter local or global analysis

! ... H(xr+dx)
! ...
do k=1,rank
  call Hforward (X(:,k),n,iobs,Nobs,HE(:,k))
enddo

! ... H(xr) = Ensemble average of HE
! ...
call emean (HE,Nobs,rank,Hxf)

! ... H'*dx = H(xr+dx) - H(xr)
! ...
do i=1,Nobs
  HE(i,:) = HE(i,:) - Hxf(i)
enddo

! ... Innovation vector:
! ...
innovation(:) = yobs(:) - Hxf(:)

! ... Second thing: Evaluate xf, and the ensemble anomalies:
! ...
call emean (X,n,rank,xf)

do i=1,n
  X(i,:) = X(i,:) - xf(i)
enddo

! ... End Global operations
! ................................................................
! ................................................................



! ... Third thing: Evaluate HER = (HE)^T R^{-1}
! ...
do j=1,Nobs
  HER(:,j) = HE(j,:)/R(j)
enddo

! ... Fourth thing: Evaluate W = HER*HE
! ... and its eigenvalue factorization:
! ...
allocate (rwork(rank,rank))
rwork = MATMUL(HER,HE)
do i=1,rank
  rwork(i,i) = xlami + rwork(i,i)
enddo

allocate (U(rank,rank))
allocate (wm(rank))
allocate (ws(rank))
verb = .true.
call eigencov (rank,rwork,U,wm,ii,verb)
if (ii.LT.rank) call stop_error(1,'LOST THE rank')

! ... Root square of the inverse ...
! ...
do i=1,rank
  ws(i) = 1.0D0 / SQRT(wm(i))
enddo

! ... Fourth thing: Evaluate aLambda = [ ((rank-1)/forget)*I + HER*HE]^-1
! ... What we do is to evaluate is square root: 
! ... Wa = SQRT(rank-1)*U Ws, aLambda = Wa * Wa^T
! ...
allocate (Wa(rank,rank))
allocate (aLambda(rank,rank))

do j=1,rank
  Wa(:,j) = U(:,j)*ws(j)      ! Root square of aLambda !!
enddo

do j=1,rank
do i=j,rank
  xsum = 0.0D0
  do k=1,rank
    xsum = xsum + Wa(i,k)*Wa(j,k)
  enddo
  aLambda(i,j) = xsum
  aLambda(j,i) = xsum
enddo
enddo

! ... Central analysis:
! ... Project the innovation vector on HER:
! ... wm is the amplitude of the modes of the central analysis:
! ...
allocate (vwork(rank))
vwork = MATMUL(HER,innovation)
wm    = MATMUL(aLambda,vwork)

! ... Cost function xJ:
! ... The first term does not depend on the amount of information projected
! ... onto the analysis subspace. The second term removes the part corrected
! ... by the analysis subspace !!!!!
! ... 
xsum = 0.0D0
do i=1,rank
  xsum = xsum + innovation(i)**2/R(i)
enddo
xJ1 = xsum

allocate (vwork2(rank))
do i=1,rank
  xsum = 0.0D0
  do j=1,rank
    xsum = xsum + Wa(j,i)*vwork(k)
  enddo
  vwork2(i) = xsum
enddo
xJ2 = - DOT_PRODUCT(vwork2,vwork2)

xJ = xJ1 + xJ2

! ... Add the central analysis (vector) to each column of Wa (matrix):
! ...
do k=1,rank
  Wa(:,k) = SQRT(rank-1.0D0)*Wa(:,k) + wm(:)
enddo

! ... The ensemble of analysis perturbations:
! ...
Xa = MATMUL(X,Wa)

! ... The ensemble of analysis states:
! ...
do k=1,rank
  Xa(:,k) = xf(:) + Xa(:,k)
enddo

call emean (Xa,n,rank,xam)


! ... Deallocate the local working arrays
! ...
deallocate (rwork)
deallocate (vwork)
deallocate (vwork2)
deallocate (wm)
deallocate (ws)
deallocate (aLambda)
deallocate (U)
deallocate (Wa)

return
end subroutine etkf
! ...
! ========================================================================
! ...
