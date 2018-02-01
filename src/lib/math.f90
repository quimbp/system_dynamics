! ****************************************************************************
! ... Mathematical utilitites
! ... COSMO Project
! ... Quim Ballabrera, March 2017
! ... About norm2: This has become a standard function after Fortran 2008.
! ****************************************************************************

module math

use types, only: dp
use constants, only: zero,one,nan
use utils

implicit none
private
public identity,akima,dakima,arange,mean,indexx,lookup_table, &
       avevar,rand,randn,randseed,eigencov,random_sample,svdcmp, &
       eigsort,matrinv
!public norm2

interface akima
  module procedure akimas,akimav
end interface akima

interface rand
  module procedure rand1,randm
end interface rand

interface randn
  module procedure randn1,randnm
end interface randn

contains
! ...
! =============================================================================
! ...
function akimas (x,f,xx) result(ff)

real(dp), dimension(:), intent(in)  :: x,f
real(dp), intent(in)                :: xx
real(dp)                            :: ff

! ... Local variables
! ...
real(dp) df(size(x))

df = dakima(x,f)
ff = evlak(xx,x,f,df)

end function akimas
! ...
! =============================================================================
! ...
function akimav (x,f,xx) result(ff)

real(dp), dimension(:), intent(in)  :: x,f
real(dp), dimension(:), intent(in)  :: xx
real(dp), dimension(size(xx))       :: ff

! ... Local variables
! ...
integer i,n,nn
real(dp) df(size(x))

n = size(x)
nn = size(xx)

df = dakima(x,f)
DO i=1,nn
  ff(i) = evlak(xx(i),x,f,df)
ENDDO

end function akimav
! ...
! =============================================================================
! ...
real(dp) function evlak(xc,x,f,df) 
! ... Given a montone increasing array of x values (Note: routine does not 
! ... check this), samples of values of the function f and 1st derivative  
! ... df at the xs, interpolates by cubic polynomial.  The array df can be 
! ... found by calling subroutine akima.                                   
!                                                                       
! ... If y falls outside the interval, the polynomial on the               
! ... nearest interval is used to extrapolate.                             

real(dp), intent(in)                  :: xc
real(dp), dimension(:), INTENT(in)    :: x,f,df

! ... Local variables
! ...
integer k,n
integer init         ! Search for proper interval initiated at previous call
save init
data init/1/
real(dp) dx,t,s

n = size(x)

! ... Locate sample interval containing xc:  after xc lies in [x(init),      
! ... x(init+1)), unless xc lies outside [x(1), x(n)] when the interval     
! ... is the intervals containing the apprpriate end point.                

init = min(init, n)
if (xc.gt.x(init)) then
  do k=init,n
    if (x(k).gt.xc) then
      init = k-1
      goto 1300
    endif
  enddo
  init = n-1
else
  do k=init,1,-1
    if (x(k).le.xc) then
      init = k
      goto 1300
    endif
  enddo
  init = 1
endif

1300 continue
dx = x(init+1) - x(init)

! ...  Evaluate the cubic interpolator                                      
! ...
t = (xc - x(init))/dx
s = 1.0D0 - t
evlak = s**2*((1.0D0 + 2.0D0*t)*f(init)   + t*dx*df(init))   +     &
        t**2*((1.0D0 + 2.0D0*s)*f(init+1) - s*dx*df(init+1))
end function evlak
! ...
! =============================================================================
! ...
function dakima (x,f) result(df)
! ... Given the array of samples of a function f and sample points x,      
! ... assumed to be monotone, generates slopes for an interpolating rule    
! ... according to Akima's algorithm (Lancaster and Salkauskas,            
! ... Curve and surface fitting, 1986 academic press, p 82).       

real(dp), dimension(:), intent(in)  :: x,f
real(dp), dimension(size(x))        :: df

! ... Local variables
! ...
integer                             :: i,n
real(dp)                            :: Sn,eps,D1,D2
real(dp), DIMENSION(4)              :: S

n = size(x)
eps  = 1.0e-6*ABS(x(n) - x(1))

S(1) = (f(2) - f(1))/(x(2) - x(1))
S(2) = S(1)
S(3) = S(1)
S(4) = (f(3) - f(2))/(x(3) - x(2))
Sn   = (f(n) - f(n-1))/(x(n) - x(n-1))

do i=1, n
  D1 = abs(S(2) - S(1))
  D2 = abs(S(4) - S(3))
  df(i) = (D2*S(2) + D1*S(3))/(D1 + D2 + eps)

  S(1) = S(2)
  S(2) = S(3)
  S(3) = S(4)
  if (i+3.le.n) then
    S(4) = (f(i+3) - f(i+2))/(x(i+3)-x(i+2))
  else
    S(4) = Sn
  endif
enddo

if (n.eq.2) return

! ... If 3 or more points use gradient from a parabola for 1st & last df   
! ...
df(1) = qakima(f(2)-f(1),x(2)-x(1),f(3)-f(1),x(3)-x(1))
df(n) = qakima(f(n-1)-f(n),x(n-1)-x(n),f(n-2)-f(n),x(n-2)-x(n))

contains
  real(dp) pure function qakima(u1,x1,u2,x2)
  real(dp), intent(in)       :: u1,x1,u2,x2
  qakima = (u1/x1**2-u2/x2**2)/(1.0_dp/x1-1.0_dp/x2)
  end function qakima

end function dakima
! ...
! =============================================================================
! ...
function arange(xo,xf,n)
! ... Returns the n-dimensional vector (xo,xo+dx,xo+2dx,...,xf)

integer, intent(in)     :: n
real(dp), intent(in)    :: xo,xf
real(dp), dimension(n)  :: arange

! ... Local variables
! ...
integer i
real(dp) dx

dx = (xf-xo)/(n-1)
do i=1,n
  arange(i) = (i-1)*dx + xo
enddo

end function arange
! ...
! =============================================================================
! ...
function identity(n)
! ... Returns the (n x n) identity matrix

integer, intent(in)          :: n
real(dp), dimension(n,n)     :: identity

! ... Local variables
! ...
integer i

identity = zero
do i=1,n
  identity(i,i) = one
enddo

end function identity
! ...
! =============================================================================
! ...
real(dp) function mean(A,W)
! ... Calculates the weighted mean = 1/N * Sum W(i)*A(i)
! ... Weights are optional.

real(dp), dimension(:), intent(in)     :: A
real(dp), dimension(:), optional       :: W

integer n
real(dp) Sw

mean = nan
n = size(A)
if (n.eq.0) return

if (present(W)) then
  Sw   = sum(W)
  mean = dot_product(W,A)/Sw
else
  mean = sum(A)/N
endif

end function mean
! ...
! =============================================================================
! ...
subroutine avevar (A,ave,var)

real(dp), dimension(:), intent(in)   :: A
real(dp), intent(out)                :: ave,var

! ... Local variables:
! ...
integer j,n
real(dp) s,ep

ave = nan
var = nan
n = size(A)
if (n.eq.0) return

ave = sum(A(1:n))/n

var = 0.0d0
ep  = 0.0d0
do j=1,n
  s   = A(j) - ave
  ep  = ep + s
  var = var + s*s
enddo
var = (var-ep**2/n)/(n-1)

return
end
!  (C) Copr. 1986-92 Numerical Recipes Software *5sV1.
! ...
! =============================================================================
! ...
!real(dp) function norm2(A)
! ... Calculates the weighted mean = 1/N * Sum W(i)*A(i)
! ... Weights are optional.

!real(dp), dimension(:), intent(in)     :: A

!norm2 = sqrt(dot_product(A,A))

!end function norm2
! ...
! =============================================================================
! ...
! ... Subroutine indexx from Numerical Recipes in Fortran.
! ... The Art of Scientific Computing. Press et al., 1992.
! ... 
      subroutine indexx(arr,indx)

      real(dp), dimension(:), intent(in)          :: arr
      integer, dimension(size(arr)), intent(out)  :: indx
      
      ! ... Local varibles:
      ! ...
      integer, parameter                          :: M=7
      INTEGER                                     :: n,i,indxt,ir,itemp,j,jstack,k,l
      INTEGER, dimension(size(arr))               :: istack
      real(dp)                                    :: a            

      n = size(arr)

      do j=1,n
        indx(j)=j
      enddo

      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do i=j-1,1,-1
            if(arr(indx(i)).le.a) goto 2
            indx(i+1)=indx(i)
          enddo
          i=0
2         indx(i+1)=indxt
        enddo
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      end subroutine indexx
!  (C) Copr. 1986-92 Numerical Recipes Software *5sV1.
! ...
! =============================================================================
! ...
  real(dp) function lookup_table (x,f,xo) result(fo)

  real(dp), dimension(:), intent(in)        :: x
  real(dp), dimension(:), intent(in)        :: f
  real(dp), intent(in)                      :: xo

  integer i1,i2,n

  n = size(x)

  if (x(1).lt.x(n)) then
    ! ... Ordered in ascending order
    if (xo.le.x(1)) then
      fo = f(1)
      return
    else if (xo.ge.x(n)) then
      fo = f(n)
    else
      i1 = locate(x(1:n),xo)
      i2 = i1 + 1
      fo = f(i1) + (xo-x(i1))/(x(i2)-x(i1))*(f(i2)-f(i1))
    endif
  else
    ! ... Ordered in descending order
    if (xo.ge.x(1)) then
      fo = f(1)
      return
    else if (xo.le.x(n)) then
      fo = f(n)
    else
      i1 = locate(x(1:n),xo)
      i2 = i1 + 1
      fo = f(i1) + (xo-x(i1))/(x(i2)-x(i1))*(f(i2)-f(i1))
    endif
  endif

  end function lookup_table
! ...
! =============================================================================
! ...
function rand1 ()

real(dp)                                   :: rand1

call random_number(rand1)  ! GNU RANDOM GENERATOR
 
end function rand1
! ...
! =============================================================================
! ...
function randm (m) result(ff)

integer, intent(in)                        :: m
real(dp), dimension(m)                     :: ff  

call random_number(ff)  ! GNU RANDOM GENERATOR
 
end function randm
! ...
! ========================================================================
! ...
function randn1 () 

real(dp)                                   :: randn1

! ... Local variables
! ...
real(dp), dimension(1)                     :: ff

ff = randnm(1)
randn1 = ff(1)

end function randn1
! ...
! ========================================================================
! ...
function randnm (m) result(ff)

integer, intent(in)                        :: m
real(dp), dimension(m)                     :: ff

! ... Local variables
! ...
integer i
reaL(dp), parameter                    :: s  =  0.449871D0
real(dp), parameter                    :: t  = -0.386595D0
real(dp), parameter                    :: a  =  0.19600D0
real(dp), parameter                    :: b  =  0.25472D0
real(dp), parameter                    :: r1 =  0.27597D0
real(dp), parameter                    :: r2 =  0.27846D0
real(dp) u,v,x,y,q

do i=1,m

  do 
    call random_number(u)  ! GNU RANDOM GENERATOR
    call random_number(v)  ! GNU RANDOM GENERATOR
    v = 1.7156D0 * (v - 0.5D0)

    ! ... Evaluate the quadratic form
    ! ...
    x = u - s
    y = abs(v) - t
    q = x**2 + y*(a*y - b*x)
 
    if (q .lt. r1) exit
    if (q .gt. r2) cycle
    if (v**2 .LT. -4D0*LOG(u)*u**2) exiT
  enddo
  ff(i) = v/u

enddo
 
end function randnm
! ...
! =============================================================================
! ...
subroutine randseed(idum)

integer, intent(in)                    :: idum

integer                                :: n
integer, dimension(:), allocatable     :: rseed

! ... Random seed
! ...
call random_seed(size=n)
allocate(rseed(n))
rseed(:) = idum
call random_seed(put=rseed)


end subroutine randseed
! ...
! =============================================================================
! ...
function random_sample (N,io,if) result(ff)

! ... Samples randomly the segment [io,if].
! ... N samples. The values of the sample are in random_sample(1:N).
! ...
integer, intent(in)                   :: N
integer, intent(in), OPTIONAL         :: io,if
integer, dimension(N)                 :: ff

integer i,il,lio,lif
real(dp) nn(1)

if (PRESENT(io)) then
  lio = io
else
  lio = 1
endif

if (PRESENT(if)) then
  lif = if
else
  lif = N
endif

il = lif - lio + 1
do i=1,N
  ff(i) = INT( io+il*rand1() )
enddo

RETURN
end function random_sample
! ...
! =============================================================================
! ...
! ***************************************************************************
! ... eigencov.f90
! ... Eigenvectors and eigenvalues of a covariance (square, symmetric) matrix.
! ...
! ... Input:
! ... N = Dimension of the P(N,N) matrix
! ... P = Input covariance matrix
! ... Output
! ... U = Eigenvectors (N,N)
! ... D = Eigenvalues (N)
! ... rank = integer with the number of non-zero eigenvalues.
! ...
! ... Uses routines tred2, tqli, and eigsort of the Numerical Recipes.
! ===========================================================================
! ...
subroutine eigencov (N,P,U,D,rank,verbose)

implicit none

INTEGER N,rank
REAL(dp) P(N,N),U(N,N)
REAL(dp) D(N)
logical verbose

INTEGER i,k
REAL(dp) xsum
REAL(dp) E(N)

if (verbose) then
  xsum = 0.0D0
  DO i=1,n
    xsum = xsum + P(i,i)
  ENDDO
  WRITE(6,*) 'Trace covariance matrix : ', xsum
endif

U = P

CALL tred2 (U,n,n,D,E)
CALL tqli (D,E,n,n,U)

DO i=1,n
  IF(D(i).LT.0.0) D(i) = 0.0D0
ENDDO

CALL eigsort(D,n,U,n)

rank = COUNT(D.NE.0)

if (verbose) then
  xsum = 0.0D0
  DO i=1,n
  DO k=1,rank
    xsum = xsum + D(k)*U(i,k)*U(i,k)
  ENDDO
  ENDDO
  WRITE(6,*) '>> TEST: Trace covariance matrix : ', xsum
  WRITE(6,*) '>> using rank: ', rank
endif

RETURN
END
! ...
! ========================================================================
! ...
      SUBROUTINE tred2(a,n,np,d,e)
      INTEGER n,np
      DOUBLE PRECISION a(np,np),d(np),e(np)
      INTEGER i,j,k,l
      DOUBLE PRECISION f,g,h,hh,scale
      do 18 i=n,2,-1
        l=i-1
        h=0.d0
        scale=0.d0
        if(l.gt.1)then
          do 11 k=1,l
            scale=scale+abs(a(i,k))
11        continue
          if(scale.eq.0.d0)then
            e(i)=a(i,l)
          else
            do 12 k=1,l
              a(i,k)=a(i,k)/scale
              h=h+a(i,k)**2
12          continue
            f=a(i,l)
            g=-sign(sqrt(h),f)
            e(i)=scale*g
            h=h-f*g
            a(i,l)=f-g
            f=0.d0
            do 15 j=1,l
!     Omit following line if finding only eigenvalues
              a(j,i)=a(i,j)/h
              g=0.d0
              do 13 k=1,j
                g=g+a(j,k)*a(i,k)
13            continue
              do 14 k=j+1,l
                g=g+a(k,j)*a(i,k)
14            continue
              e(j)=g/h
              f=f+e(j)*a(i,j)
15          continue
            hh=f/(h+h)
            do 17 j=1,l
              f=a(i,j)
              g=e(j)-hh*f
              e(j)=g
              do 16 k=1,j
                a(j,k)=a(j,k)-f*e(k)-g*a(i,k)
16            continue
17          continue
          endif
        else
          e(i)=a(i,l)
        endif
        d(i)=h
18    continue
!     Omit following line if finding only eigenvalues.
      d(1)=0.d0
      e(1)=0.d0
      do 24 i=1,n
!     Delete lines from here ...
        l=i-1
        if(d(i).ne.0.d0)then
          do 22 j=1,l
            g=0.d0
            do 19 k=1,l
              g=g+a(i,k)*a(k,j)
19          continue
            do 21 k=1,l
              a(k,j)=a(k,j)-g*a(k,i)
21          continue
22        continue
        endif
!     ... to here when finding only eigenvalues.
        d(i)=a(i,i)
!     Also delete lines from here ...
        a(i,i)=1.d0
        do 23 j=1,l
          a(i,j)=0.d0
          a(j,i)=0.d0
23      continue
!     ... to here when finding only eigenvalues.
24    continue
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software *5sV1.
! ...
! ========================================================================
! ...
      SUBROUTINE tqli(d,e,n,np,z)
      INTEGER n,np
      DOUBLE PRECISION d(np),e(np),z(np,np)
!CU    USES pythag
      INTEGER i,iter,k,l,m
      DOUBLE PRECISION b,c,dd,f,g,p,r,s
      do 11 i=2,n
        e(i-1)=e(i)
11    continue
      e(n)=0.d0
      do 15 l=1,n
        iter=0
1       do 12 m=l,n-1
          dd=abs(d(m))+abs(d(m+1))
          if (abs(e(m))+dd.eq.dd) goto 2
12      continue
        m=n
2       if(m.ne.l)then
          if(iter.eq.40) STOP 'too many iterations in tqli'
          iter=iter+1
          g=(d(l+1)-d(l))/(2.d0*e(l))
          r=pythag(g,1.d0)
          g=d(m)-d(l)+e(l)/(g+sign(r,g))
          s=1.d0
          c=1.d0
          p=0.d0
          do 14 i=m-1,l,-1
            f=s*e(i)
            b=c*e(i)
            r=pythag(f,g)
            e(i+1)=r
            if(r.eq.0.d0)then
              d(i+1)=d(i+1)-p
              e(m)=0.d0
              goto 1
            endif
            s=f/r
            c=g/r
            g=d(i+1)-p
            r=(d(i)-g)*s+2.d0*c*b
            p=s*r
            d(i+1)=g+p
            g=c*r-b
!     Omit lines from here ...
            do 13 k=1,n
              f=z(k,i+1)
              z(k,i+1)=s*z(k,i)+c*f
              z(k,i)=c*z(k,i)-s*f
13          continue
!     ... to here when finding only eigenvalues.
14        continue
          d(l)=d(l)-p
          e(l)=g
          e(m)=0.d0
          goto 1
        endif
15    continue
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software *5sV1.
! ...
! =====================================================================
! ...
subroutine eigsort(d,r,v,n)

implicit none

integer, intent(in)                       :: n,r
real(dp), dimension(r), intent(inout)     :: d
real(dp), dimension(n,r), intent(inout)   :: v

! ... Local variables
! ...
integer i,j,k
real(dp) p

do i=1,r-1
  k=i
  p=d(i)
  do j=i+1,r
    if (d(j).ge.p) then
      k=j
      p=d(j)
    endif
  enddo
  if (k.ne.i) then
    d(k)=d(i)
    d(i)=p
    do j=1,n
      p=v(j,i)
      v(j,i)=v(j,k)
      v(j,k)=p
    enddo
  endif
enddo

return
end subroutine eigsort
! ...
! =============================================================================
! ...
real(dp) pure function pythag(a,b)

implicit none

real(dp), intent(in)                    :: a,b

! ... Local variables
! ...
real(dp) absa,absb

absa=abs(a)
absb=abs(b)
if(absa.gt.absb)then
  pythag=absa*sqrt(1.d0+(absb/absa)**2)
else
  if(absb.eq.0.d0)then
    pythag=0.d0
  else
    pythag=absb*sqrt(1.d0+(absa/absb)**2)
  endif
endif
return
end function pythag
!  (C) Copr. 1986-92 Numerical Recipes Software *5sV1.
! ...
! =============================================================================
! ...
      SUBROUTINE svdcmp(a,m,n,mp,np,w,v)

      implicit none

      integer m,mp,n,np
      real(dp) a(mp,np),v(np,np),w(np)
!CU    USES pythag
      integer i,its,j,jj,k,l,nm
      real(dp) anorm,c,f,g,h,s,scale,x,y,z,rv1(n)

      g=0.0d0
      scale=0.0d0
      anorm=0.0d0
      do 25 i=1,n
        l=i+1
        rv1(i)=scale*g
        g=0.0d0
        s=0.0d0
        scale=0.0d0
        if(i.le.m)then
          do 11 k=i,m
            scale=scale+abs(a(k,i))
11        continue
          if(scale.ne.0.0d0)then
            do 12 k=i,m
              a(k,i)=a(k,i)/scale
              s=s+a(k,i)*a(k,i)
12          continue
            f=a(i,i)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,i)=f-g
            do 15 j=l,n
              s=0.0d0
              do 13 k=i,m
                s=s+a(k,i)*a(k,j)
13            continue
              f=s/h
              do 14 k=i,m
                a(k,j)=a(k,j)+f*a(k,i)
14            continue
15          continue
            do 16 k=i,m
              a(k,i)=scale*a(k,i)
16          continue
          endif
        endif
        w(i)=scale *g
        g=0.0d0
        s=0.0d0
        scale=0.0d0
        if((i.le.m).and.(i.ne.n))then
          do 17 k=l,n
            scale=scale+abs(a(i,k))
17        continue
          if(scale.ne.0.0d0)then
            do 18 k=l,n
              a(i,k)=a(i,k)/scale
              s=s+a(i,k)*a(i,k)
18          continue
            f=a(i,l)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,l)=f-g
            do 19 k=l,n
              rv1(k)=a(i,k)/h
19          continue
            do 23 j=l,m
              s=0.0d0
              do 21 k=l,n
                s=s+a(j,k)*a(i,k)
21            continue
              do 22 k=l,n
                a(j,k)=a(j,k)+s*rv1(k)
22            continue
23          continue
            do 24 k=l,n
              a(i,k)=scale*a(i,k)
24          continue
          endif
        endif
        anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
25    continue
      do 32 i=n,1,-1
        if(i.lt.n)then
          if(g.ne.0.0d0)then
            do 26 j=l,n
              v(j,i)=(a(i,j)/a(i,l))/g
26          continue
            do 29 j=l,n
              s=0.0d0
              do 27 k=l,n
                s=s+a(i,k)*v(k,j)
27            continue
              do 28 k=l,n
                v(k,j)=v(k,j)+s*v(k,i)
28            continue
29          continue
          endif
          do 31 j=l,n
            v(i,j)=0.0d0
            v(j,i)=0.0d0
31        continue
        endif
        v(i,i)=1.0d0
        g=rv1(i)
        l=i
32    continue
      do 39 i=min(m,n),1,-1
        l=i+1
        g=w(i)
        do 33 j=l,n
          a(i,j)=0.0d0
33      continue
        if(g.ne.0.0d0)then
          g=1.0d0/g
          do 36 j=l,n
            s=0.0d0
            do 34 k=l,m
              s=s+a(k,i)*a(k,j)
34          continue
            f=(s/a(i,i))*g
            do 35 k=i,m
              a(k,j)=a(k,j)+f*a(k,i)
35          continue
36        continue
          do 37 j=i,m
            a(j,i)=a(j,i)*g
37        continue
        else
          do 38 j= i,m
            a(j,i)=0.0d0
38        continue
        endif
        a(i,i)=a(i,i)+1.0d0
39    continue
      do 49 k=n,1,-1
        do 48 its=1,30
          do 41 l=k,1,-1
            nm=l-1
            if((abs(rv1(l))+anorm).eq.anorm)  goto 2
            if((abs(w(nm))+anorm).eq.anorm)  goto 1
41        continue
1         c=0.0d0
          s=1.0d0
          do 43 i=l,k
            f=s*rv1(i)
            rv1(i)=c*rv1(i)
            if((abs(f)+anorm).eq.anorm) goto 2
            g=w(i)
            h=pythag(f,g)
            w(i)=h
            h=1.0d0/h
            c= (g*h)
            s=-(f*h)
            do 42 j=1,m
              y=a(j,nm)
              z=a(j,i)
              a(j,nm)=(y*c)+(z*s)
              a(j,i)=-(y*s)+(z*c)
42          continue
43        continue
2         z=w(k)
          if(l.eq.k)then
            if(z.lt.0.0d0)then
              w(k)=-z
              do 44 j=1,n
                v(j,k)=-v(j,k)
44            continue
            endif
            goto 3
          endif
          if(its.eq.30) stop  'no convergence in svdcmp'
          x=w(l)
          nm=k-1
          y=w(nm)
          g=rv1(nm)
          h=rv1(k)
          f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0d0*h*y)
          g=pythag(f,1.0d0)
          f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
          c=1.0d0
          s=1.0d0
          do 47 j=l,nm
            i=j+1
            g=rv1(i)
            y=w(i)
            h=s*g
            g=c*g
            z=pythag(f,h)
            rv1(j)=z
            c=f/z
            s=h/z
            f= (x*c)+(g*s)
            g=-(x*s)+(g*c)
            h=y*s
            y=y*c
            do 45 jj=1,n
              x=v(jj,j)
              z=v(jj,i)
              v(jj,j)= (x*c)+(z*s)
              v(jj,i)=-(x*s)+(z*c)
45          continue
            z=pythag(f,h)
            w(j)=z
            if(z.ne.0.0d0)then
              z=1.0d0/z
              c=f*z
              s=h*z
            endif
            f= (c*g)+(s*y)
            x=-(s*g)+(c*y)
            do 46 jj=1,m
              y=a(jj,j)
              z=a(jj,i)
              a(jj,j)= (y*c)+(z*s)
              a(jj,i)=-(y*s)+(z*c)
46          continue
47        continue
          rv1(l)=0.0d0
          rv1(k)=f
          w(k)=x
48      continue
3       continue
49    continue
      return
      end subroutine svdcmp
!C  (C) Copr. 1986-92 Numerical Recipes Software *5sV1.
! ...
! =============================================================================
! ...
subroutine matrinv (N,A,Y)

implicit none

integer, intent(in)                        :: N
REAL(dp), dimension(N,N), intent(in)       :: A
REAL(dp), dimension(N,N), intent(out)      :: Y

integer i,j
integer, DIMENSION (N)                     :: indx
REAL(dp) D
REAL(dp), DIMENSION (N,N)                  :: AA

! ... Saving matrix
! ...
AA = A

Y(:,:) = 0.0D0
DO i=1,N
  Y(i,i) = 1.0D0
ENDDO

call ludcmp (AA,N,N,INDX,D)

do j=1,N
  calL lubksb (AA,N,N,INDX,Y(1,j))
enddo

RETURN
END SUBROUTINE matrinv
! ...
! ================================================================================
! ...
subroutine ludcmp(a,n,np,indx,d)

implicit none

real(dp), parameter                           :: TINY=1.0D-20

integer, intent(in)                           :: n,np
real(dp), dimension(np,np), intent(inout)     :: a
integer, dimension(n), intent(out)            :: indx
real(dp), intent(out)                         :: d

! ... Local variables
! ...
integer i,imax,j,k
real(dp) aamax,dum,sum,vv(n)

d = 1.d0
do i=1,n
  aamax = 0.d0
  do j=1,n
    if (abs(a(i,j)).gt.aamax) aamax = abs(a(i,j))
  enddo
  if (aamax.eq.0.d0) stop  'singular matrix in ludcmp'
  vv(i) = 1.d0/aamax
enddo

do j=1,n
  do i=1,j-1
    sum = a(i,j)
    do k=1,i-1
      sum = sum - a(i,k)*a(k,j)
    enddo
    a(i,j) = sum
  enddo
  aamax = 0.d0
  do i=j,n
    sum = a(i,j)
    do k=1,j-1
      sum = sum-a(i,k)*a(k,j)
    enddo
    a(i,j) = sum
    dum = vv(i)*abs(sum)
    if (dum.ge.aamax) then
      imax  = i
      aamax = dum
    endif
  enddo
  if (j.ne.imax)then
    do k=1,n
      dum       = a(imax,k)
      a(imax,k) = a(j,k)
      a(j,k)    = dum
    enddo
    d = -d
    vv(imax) = vv(j)
  endif
  indx(j) = imax
  if(a(j,j).eq.0.d0) a(j,j) = TINY
  if(j.ne.n)then
    dum = 1.d0/a(j,j)
    do i=j+1,n
      a(i,j) = a(i,j)*dum
    enddo
  endif
enddo

return
end subroutine ludcmp
!C  (C) Copr. 1986-92 Numerical Recipes Software *5sV1.
! ...
! ================================================================================
! ...
subroutine lubksb(a,n,np,indx,b)

implicit none

integer, intent(in)                    :: n,np
integer, dimension(n), intent(in)      :: indx
real(dp), dimension(np,np), intent(in) :: a
real(dp), dimension(n), intent(out)    :: b

! ... Local variables
! ...
integer i,ii,j,ll
real(dp) sum

ii=0
do i=1,n
  ll=indx(i)
  sum=b(ll)
  b(ll)=b(i)
  if (ii.ne.0)then
    do j=ii,i-1
      sum=sum-a(i,j)*b(j)
    enddo
  else if (sum.ne.0.d0) then
    ii=i
  endif
  b(i)=sum
enddo
do i=n,1,-1
  sum=b(i)
  do j=i+1,n
    sum=sum-a(i,j)*b(j)
  enddo
  b(i)=sum/a(i,i)
enddo
return
end subroutine lubksb
!  (C) Copr. 1986-92 Numerical Recipes Software *5sV1.

end module math
