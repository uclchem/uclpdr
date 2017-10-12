!=======================================================================
!
!  Calculate the photon escape probability (beta) using the expression
!  of de Jong, Dalgarno & Chu (1975, ApJ, 199, 69) taking into account
!  the escape probability along each ray to the PDR surface.
!
!-----------------------------------------------------------------------
FUNCTION ESCAPE_PROBABILITY(NRAYS,TAU) RESULT(BETA)

   USE HEALPIX_TYPES
   IMPLICIT NONE

   INTEGER(KIND=I4B), INTENT(IN) :: NRAYS
   REAL(KIND=DP),     INTENT(IN) :: TAU(0:NRAYS-1)

   INTEGER(KIND=I4B) :: K
   REAL(KIND=DP)     :: BETA,BETA_RAY(0:NRAYS-1)

!  Initialize the escape probability values along each ray
   BETA_RAY=0.0D0

!  Calculate the escape probability along each ray
   DO K=0,NRAYS-1 ! Loop over rays

!!$!     Prevent exploding beta values caused by strong masing (tau << 0)
!!$!     Impose tau = -5 and calculate the escape probability accordingly
!!$      IF(TAU(K).LT.-5.0D0) THEN
!!$         BETA_RAY(K)=(1.0D0-EXP(5.0D0))/(-5.0D0)

!     Limit the escape probability to unity for masing transitions (tau <= 0)
      IF(TAU(K).LE.0) THEN
         BETA_RAY(K)=1.0D0

!     Prevent floating point overflow caused by very low opacity (tau < 1E-8)
      ELSE IF(ABS(TAU(K)).LT.1.0D-8) THEN
         BETA_RAY(K)=1.0D0

!     For all other cases use the standard escape probability formalism
      ELSE
         BETA_RAY(K)=(1.0D0-EXP(-TAU(K)))/TAU(K)
      END IF

   END DO ! End of loop over rays

!  The total escape probability must be divided by the number of rays to
!  account for the fraction of the total solid angle covered by each ray
!  (assuming that each ray covers the same fraction of the total 4π sr).
!  In the case of only 1 ray (i.e., semi-infinite slab geometry) the ray
!  subtends a solid angle of 2π sr, since the photons escape through the
!  hemisphere in the outward direction, so its escape probability should
!  be divided by two.
   BETA=SUM(BETA_RAY)/REAL(MIN(NRAYS,2),KIND=DP)

   RETURN
END FUNCTION ESCAPE_PROBABILITY
!=======================================================================
