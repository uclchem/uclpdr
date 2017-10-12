!=======================================================================
!
!  Calculate the level populations at LTE for the given species.
!
!-----------------------------------------------------------------------
SUBROUTINE CALCULATE_LTE_POPULATIONS(NLEVEL,ENERGY,WEIGHT,DENSITY, &
                                   & TEMPERATURE,POPULATION)

   USE HEALPIX_TYPES
   IMPLICIT NONE

   INTEGER(KIND=I4B), INTENT(IN)  :: NLEVEL
   REAL(KIND=DP),     INTENT(IN)  :: ENERGY(:),WEIGHT(:)
   REAL(KIND=DP),     INTENT(IN)  :: DENSITY,TEMPERATURE
   REAL(KIND=DP),     INTENT(OUT) :: POPULATION(:)

   INTEGER(KIND=I4B) :: ILEVEL
   REAL(KIND=DP) :: PARTITION_FUNCTION

!  Initialize the level populations
   POPULATION=0.0D0

   PARTITION_FUNCTION=0.0D0
   DO ILEVEL=1,NLEVEL
      POPULATION(ILEVEL)=WEIGHT(ILEVEL)*EXP(-ENERGY(ILEVEL)/(KB*TEMPERATURE))
      PARTITION_FUNCTION=PARTITION_FUNCTION+POPULATION(ILEVEL)
   END DO
   POPULATION=POPULATION*DENSITY/PARTITION_FUNCTION

!  Check that the sum of the level populations matches the total density to within 0.1%
   IF(ABS(DENSITY-SUM(POPULATION))/DENSITY.GT.1.0D-3) THEN
      WRITE(6,"('ERROR! Sum of LTE level populations differs from the total density by',F4.1,'%')") &
         & 1.0D2*ABS(SUM(POPULATION)-DENSITY)/DENSITY
      STOP
   END IF

   RETURN
END SUBROUTINE CALCULATE_LTE_POPULATIONS
!=======================================================================
