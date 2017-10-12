!=======================================================================
! DEFINITIONS.F90
! C. P. Batty & D. A. Hubber - 8/12/2006
! Definition of precision kind variables
!=======================================================================
MODULE DEFINITIONS

  integer, parameter :: DP = selected_real_kind(p=15) ! double precision
  integer, parameter :: SP = selected_real_kind(p=6)  ! single precision

#ifdef DOUBLE_PRECISION
  integer, parameter :: PR = DP                       ! particle precision
#else
  integer, parameter :: PR = SP                       ! default = single
#endif

  integer, parameter :: ILP = 4                       ! Integer long precision

END MODULE DEFINITIONS
!=======================================================================
