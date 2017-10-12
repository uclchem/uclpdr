!=======================================================================
!
!  Check for convergence in the relative abundances of all species
!  of each particle. Set the relevant convergence flags. Calculate
!  the percentage of particles that have converged.
!
!-----------------------------------------------------------------------
SUBROUTINE CHECK_CHEMISTRY_CONVERGENCE(NPART,NSPEC,PARTICLE,PERCENTAGE_CONVERGED)

   USE HEALPIX_TYPES
   USE PARTICLE_MODULE
   USE MAIN_MODULE, ONLY : ABUNDANCE_LIMIT,CHEMISTRY_CONVERGENCE_CRITERION

   IMPLICIT NONE

   INTEGER(KIND=I4B),   INTENT(IN)    :: NPART,NSPEC
   TYPE(PARTICLE_TYPE), INTENT(INOUT) :: PARTICLE(:)
   REAL(KIND=DP),       INTENT(OUT)   :: PERCENTAGE_CONVERGED

   INTEGER(KIND=I4B) :: I,P
   REAL(KIND=DP) :: RELATIVE_CHANGE

   PARTICLE%CHEMISTRY_CONVERGED=.TRUE.
   DO P=1,NPART ! Loop over particles
      DO I=1,NSPEC ! Loop over species

!        Skip this species if its abundance is below the cut-off
         IF(PARTICLE(P)%ABUNDANCE(I).LT.ABUNDANCE_LIMIT) CYCLE

!        Skip this species if its abundance has not changed
         IF(PARTICLE(P)%ABUNDANCE(I).EQ.PARTICLE(P)%PREVIOUS_ABUNDANCE(I)) CYCLE

!        Calculate the relative change in abundance between this iteration and the previous
         RELATIVE_CHANGE=ABS(PARTICLE(P)%ABUNDANCE(I)-PARTICLE(P)%PREVIOUS_ABUNDANCE(I)) &
                       & *2/(PARTICLE(P)%ABUNDANCE(I)+PARTICLE(P)%PREVIOUS_ABUNDANCE(I))

!        If the relative change is greater than the criterion for convergence, set the flag to false
         IF(RELATIVE_CHANGE.GT.CHEMISTRY_CONVERGENCE_CRITERION) THEN
            PARTICLE(P)%CHEMISTRY_CONVERGED=.FALSE.
            EXIT
         END IF

      END DO ! End of loop over species
   END DO ! End of loop over particles

!  Calculate the percentage of particles whose abundances have converged
   PERCENTAGE_CONVERGED=COUNT(PARTICLE%CHEMISTRY_CONVERGED)/REAL(NPART,KIND=DP)*100

   RETURN
END SUBROUTINE CHECK_CHEMISTRY_CONVERGENCE
!=======================================================================
