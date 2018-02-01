! *****************************************************************************
! ... obs.f90
! ... System Dynamics
! ... Quim Ballabrera, June 2009
! *****************************************************************************

module observations

use types, only: dp
use constants, only: nan

implicit none

type obs_type
  integer                           :: n              ! Number of obs
  integer                           :: Nobs           ! Number time records
  integer                           :: p = 1          ! Observation to be used
  real(dp)                          :: missing = nan  ! Optional missing value
  real(dp), dimension(:), pointer   :: time           ! Observational time
  real(dp), dimension(:,:), pointer :: y              ! Observed values
  integer, dimension(:), pointer    :: step           ! Associated model step
  integer, dimension(:), pointer    :: ind            ! Observed Model-index 
  integer, dimension(:), pointer    :: iobs           ! 1/0 observation flag multiplier
end type obs_type

contains
! ...
! =============================================================================
! ...
function obs_read (filename) result (obs)

use utils, only: unitfree,numlines,numwords

implicit none

character(len=*), intent(in)          :: filename
type(obs_type)                        :: obs

integer iu,i
character(len=180) line

write(*,*)
write(*,*) 'Reading observation file : ', trim(filename)

iu = unitfree()
open(iu,file=filename,status='old')
obs%Nobs = numlines(iu)

rewind(iu)
read(iu,'(A)') line
obs%n = numwords(line) - 1    ! First column is time

write(*,*) 'Number of columns    = ', obs%n + 1
write(*,*) 'Number of lines      = ', obs%Nobs

allocate(obs%time(obs%Nobs))
allocate(obs%y(obs%n,obs%Nobs))
allocate(obs%step(obs%Nobs))
allocate(obs%iobs(obs%n))
obs%iobs(:) = 1                ! By default we assimilate everything

rewind(iu)
do i=1,obs%Nobs
  read(iu,*,err=100,end=100) obs%time(i), obs%y(:,i)
enddo
close(iu)
return

100 continue
write(*,*) 'Skipping lines. Nobs = ', i - 1
obs%Nobs = i - 1
close(iu)
return


end function obs_read
! ...
! =============================================================================
! ...
subroutine obs_map (obs,to,dt,nsteps)

implicit none

type(obs_type), intent(inout)         :: obs
real(dp), intent(in)                  :: to,dt
integer, intent(in)                   :: nsteps

! ... Local variables
! ...
integer i,kk
real(dp) tt

obs%step(:) = -1
kk = 1
do i=1,nsteps+1
 tt = to + (i-1)*dt
 if (abs(obs%time(kk)-tt).le.0.1*dt) then
   obs%step(kk) = i - 1
   kk = kk + 1
   if (kk.gt.obs%Nobs) exit
 endif
enddo

end subroutine obs_map
! ...
! =============================================================================
! ...
end module observations
