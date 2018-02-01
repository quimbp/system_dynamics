! ****************************************************************************
! ... Minimization utilitites
! ... System Dynamics
! ... Quim Ballabrera, January 2018
! ... Source: Numerical Recipes
! ...
! ... subroutine DFPMIN
! ... subroutine LNSRCH
! ****************************************************************************

module minimization

use types, only: dp

implicit none

contains

SUBROUTINE dfpmin(p,n,gtol,iter,fret,func,dfunc)

integer, parameter                         :: ITMAX = 200
real(dp), parameter                        :: STPMX = 100.0D0
real(dp), parameter                        :: EPS   = 3.0D-8
real(dp), parameter                        :: TOLX  = 4.0D0*EPS


integer, intent(in)                        :: n
real(dp), dimension(n), intent(inout)      :: p
integer, intent(out)                       :: iter
real(dp), intent(in)                       :: gtol
real(dp), intent(out)                      :: fret

real(dp)                                   :: func
external                                   :: dfunc,func

! ... Local variables
! ...
integer i,its,j
logical check

real(dp) den,fac,fad,fae,fp,stpmax,xsum,sumdg,sumxi,temp,test
real(dp), dimension(n)      :: dg,g,hdg,pnew,xi
real(dp), dimension(n,n)    :: hessin

fp=func(p)
call dfunc(p,g)

xsum=0.d0
hessin(:,:) = 0.0d0
do i=1,n
  hessin(i,i) = 1.d0
  xi(i)       = -g(i)
  xsum         = xsum + p(i)**2
enddo
stpmax=STPMX*max(sqrt(xsum), dble(n))

do its=1,ITMAX
  iter=its
  call lnsrch(n,p,fp,g,xi,pnew,fret,stpmax,check,func)
  fp=fret

  xi(:) = pnew(:) - p(:)
  p(:)  = pnew(:)

  test=0.d0
  do i=1,n
    temp = abs(xi(i))/max(abs(p(i)),1.d0)
    if (temp.gt.test) test = temp
  enddo
  if (test.lt.TOLX) return

  dg(:) = g(:)
  call dfunc(p,g)
  test = 0.d0
  den  = max(fret,1.d0)
  do i=1,n
    temp = abs(g(i))*max(abs(p(i)),1.d0)/den
    if (temp.gt.test) test = temp
  enddo
  if(test.lt.gtol) return

  dg(:) = g(:) - dg(:)
  do i=1,n
    hdg(i) = DOT_PRODUCT(hessin(i,:),dg(:))
  enddo

  fac   = 0.d0
  fae   = 0.d0
  sumdg = 0.d0
  sumxi = 0.d0
  do i=1,n
    fac   = fac   + dg(i)*xi(i)
    fae   = fae   + dg(i)*hdg(i)
    sumdg = sumdg + dg(i)**2
    sumxi = sumxi + xi(i)**2
  enddo

  if(fac**2.gt.EPS*sumdg*sumxi)then
    fac=1.d0/fac
    fad=1.d0/fae
    do i=1,n
      dg(i)=fac*xi(i)-fad*hdg(i)
    enddo
    do i=1,n
      do j=1,n
        hessin(i,j) = hessin(i,j)       &
                    + fac*xi(i)*xi(j)   &
                    - fad*hdg(i)*hdg(j) &
                    + fae*dg(i)*dg(j)
      enddo
    enddo
  endif

  do i=1,n
    xi(i) = -DOT_PRODUCT(hessin(i,:),g(:))
  enddo

enddo
stop  'Error: too many iterations in dfpmin'

end subroutine dfpmin
!C  (C) Copr. 1986-92 Numerical Recipes Software *5sV1.
! ...
! ====================================================================================
! ...
subroutine lnsrch(n,xold,fold,g,p,x,f,stpmax,check,func)

real(dp), parameter                    :: ALF  = 1.0D-4
real(dp), parameter                    :: TOLX = 1.0D-7

integer, intent(in)                    :: n
real(dp), intent(in)                   :: fold
real(dp), intent(in)                   :: stpmax
real(dp), dimension(n), intent(in)     :: xold
real(dp), dimension(n), intent(in)     :: g
real(dp), dimension(n), intent(inout)  :: p
real(dp), dimension(n), intent(out)    :: x
real(dp), intent(out)                  :: f
logical, intent(out)                   :: check

real(dp)                               :: func
external                               :: func

! ... Local variables
! ...
integer i
real(dp) a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2
real(dp) slope,xsum,temp,test,tmplam

check=.false.
xsum = sqrt(DOT_PRODUCT(p,p))
if (xsum.gt.stpmax) then
  p(:) = p(:)*stpmax/xsum
endif
slope = DOT_PRODUCT(g,p)
test=0.d0
do i=1,n
  temp=abs(p(i))/max(abs(xold(i)),1.d0)
  if(temp.gt.test)test=temp
enddo
alamin=TOLX/test
alam=1.d0
do
  x(:) = xold(:) + alam*p(:)
  f = func(x)
  if(alam.lt.alamin)then
    x(:) = xold(:)
    check=.true.
    return
  else if(f.le.fold+ALF*alam*slope)then
    return
  else
    if(alam.eq.1.d0)then
      tmplam=-slope/(2.d0*(f-fold-slope))
    else
      rhs1=f-fold-alam*slope
      rhs2=f2-fold2-alam2*slope
      a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
      b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
      if(a.eq.0.d0)then
        tmplam=-slope/(2.d0*b)
      else
        disc=b*b-3.d0*a*slope
        tmplam=(-b+sqrt(disc))/(3.d0*a)
      endif
      if(tmplam.gt..5d0*alam)tmplam=.5d0*alam
    endif
  endif
  alam2=alam
  f2=f
  fold2=fold
  alam=max(tmplam,.1d0*alam)
enddo
end subroutine lnsrch
!C  (C) Copr. 1986-92 Numerical Recipes Software *5sV1.
! ...
! ====================================================================================
! ...

end module minimization
