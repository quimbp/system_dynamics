! SYSDYN module
! Quim Ballabrera, Novemeber 2017

module sysdyn

  use types, only: sp,dp
  use constants
  use utils
  use lineargs
  use help
  use lorenz
  use math
  use runge_kutta
  use dates
  use observations
  use minimization

END MODULE sysdyn
