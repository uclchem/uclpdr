!=======================================================================
!
!  Check for convergence in the population densities of all levels of
!  each coolant of each particle. Set the relevant convergence flags.
!  Calculate the percentage of particles that have converged for each
!  coolant.
!
!-----------------------------------------------------------------------
SUBROUTINE CHECK_POPULATION_CONVERGENCE(NPART,NCOOL,COOLANT,PARTICLE,PERCENTAGE_CONVERGED)

   USE HEALPIX_TYPES
   USE COOLANT_MODULE
   USE PARTICLE_MODULE
   USE MAIN_MODULE, ONLY : POPULATION_LIMIT,POPULATION_CONVERGENCE_CRITERION

   IMPLICIT NONE

   INTEGER(KIND=I4B),   INTENT(IN)    :: NPART,NCOOL
   TYPE(COOLANT_TYPE),  INTENT(INOUT) :: COOLANT(:)
   TYPE(PARTICLE_TYPE), INTENT(INOUT) :: PARTICLE(:)
   REAL(KIND=DP),       INTENT(OUT)   :: PERCENTAGE_CONVERGED(:)

   INTEGER(KIND=I4B) :: I,N,P
   REAL(KIND=DP) :: RELATIVE_CHANGE

   DO P=1,NPART ! Loop over particles
      PARTICLE(P)%COOLANT%CONVERGED=.TRUE.
      DO N=1,NCOOL ! Loop over coolants
         DO I=1,COOLANT(N)%NLEVEL ! Loop over levels

!           Skip this level if its population density is below the cut-off
            IF(PARTICLE(P)%COOLANT(N)%POPULATION(I).LT.POPULATION_LIMIT*PARTICLE(P)%COOLANT(N)%DENSITY) CYCLE

!           Skip this level if its population density has not changed
            IF(PARTICLE(P)%COOLANT(N)%POPULATION(I).EQ.PARTICLE(P)%COOLANT(N)%PREVIOUS_POPULATION(I)) CYCLE

!           Calculate the relative change in population density between this iteration and the previous
            RELATIVE_CHANGE=ABS(PARTICLE(P)%COOLANT(N)%POPULATION(I)-PARTICLE(P)%COOLANT(N)%PREVIOUS_POPULATION(I)) &
                          & *2/(PARTICLE(P)%COOLANT(N)%POPULATION(I)+PARTICLE(P)%COOLANT(N)%PREVIOUS_POPULATION(I))

!           If the relative change is greater than the criterion for convergence, set the flag to false
            IF(RELATIVE_CHANGE.GT.POPULATION_CONVERGENCE_CRITERION) THEN
               PARTICLE(P)%COOLANT(N)%CONVERGED=.FALSE.
               EXIT
            END IF

         END DO ! End of loop over levels
      END DO ! End of loop over coolants
   END DO ! End of loop over particles

!  Calculate the percentage of particles whose populations have converged for each coolant
   PERCENTAGE_CONVERGED=0
   DO N=1,NCOOL
      DO P=1,NPART
         IF(PARTICLE(P)%COOLANT(N)%CONVERGED) PERCENTAGE_CONVERGED(N)=PERCENTAGE_CONVERGED(N)+1
      END DO
      PERCENTAGE_CONVERGED(N)=PERCENTAGE_CONVERGED(N)/REAL(NPART,KIND=DP)*100
      IF(INT(PERCENTAGE_CONVERGED(N)).GE.100) COOLANT(N)%CONVERGED=.TRUE.
   END DO

   RETURN
END SUBROUTINE CHECK_POPULATION_CONVERGENCE
!=======================================================================
