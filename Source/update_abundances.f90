!=======================================================================
!
!  Determine the list of particles whose abundances need to be updated
!  on the current iteration and load their properties into dummy arrays
!  to be passed to the chemistry routine. Reset their abundances to the
!  appropriate starting values before updating them.
!
!-----------------------------------------------------------------------
SUBROUTINE UPDATE_ABUNDANCES(TEMPERATURE_ITERATION,CHEMISTRY_ITERATION,NPART,PARTICLE,INITIAL_ABUNDANCE, &
                           & NCHEM,DUMMY_INDEX,DUMMY_ABUNDANCE,DUMMY_RATE,DUMMY_DENSITY,DUMMY_TEMPERATURE)

   USE HEALPIX_TYPES
   USE PARTICLE_MODULE

   IMPLICIT NONE

   INTEGER(KIND=I4B),   INTENT(IN)  :: TEMPERATURE_ITERATION,CHEMISTRY_ITERATION,NPART
   TYPE(PARTICLE_TYPE), INTENT(INOUT)  :: PARTICLE(*)
   REAL(KIND=DP),       INTENT(IN)  :: INITIAL_ABUNDANCE(:)
   INTEGER(KIND=I4B),   INTENT(OUT) :: NCHEM,DUMMY_INDEX(:)
   REAL(KIND=DP),       INTENT(OUT) :: DUMMY_ABUNDANCE(:,:),DUMMY_RATE(:,:)
   REAL(KIND=DP),       INTENT(OUT) :: DUMMY_DENSITY(:),DUMMY_TEMPERATURE(:)

   INTEGER(KIND=I4B) :: I,P
   REAL(KIND=DP)     :: TEMPERATURE_DIFFERENCE,TEMPERATURE_DIFFERENCE_LIMIT

   TEMPERATURE_DIFFERENCE_LIMIT = 0.1D0 ! Specify a temperature difference limit of 0.1 K

!  Initialize the list of particles to be updated and the associated dummy arrays
   DUMMY_INDEX = 0
   DUMMY_DENSITY = 0.0D0
   DUMMY_TEMPERATURE = 0.0D0
   DUMMY_ABUNDANCE = 0.0D0
   DUMMY_RATE = 0.0D0

!  Determine the list of particles whose abundances are to be updated
   I = 0 ! Set the list index counter to zero

!  Calculate abundances for all non-ionized particles on the first iteration
   IF(TEMPERATURE_ITERATION.EQ.1 .AND. CHEMISTRY_ITERATION.EQ.1) THEN
      DO P=1,NPART ! Loop over particles
         IF(PARTICLE(P)%TYPE.EQ.1) CYCLE ! Exclude particles in the ionized HII region
         I = I+1
         DUMMY_INDEX(I) = P
      END DO ! End of loop over particles

!  Update only the PDR particles (whose abundances may depend on shielding-dependent photorates)
   ELSE IF(TEMPERATURE_ITERATION.EQ.1 .AND. CHEMISTRY_ITERATION.GT.1) THEN
      DO P=1,NPART ! Loop over particles
         IF(PARTICLE(P)%TYPE.EQ.1) CYCLE ! Exclude particles in the ionized HII region
         IF(PARTICLE(P)%TYPE.EQ.3) CYCLE ! Exclude particles in the dark cloud region
         I = I+1
         DUMMY_INDEX(I) = P
      END DO ! End of loop over particles

!  Update all non-ionized particles whose temperatures have changed by more than the specified limit
   ELSE IF(TEMPERATURE_ITERATION.GT.1 .AND. CHEMISTRY_ITERATION.EQ.1) THEN
      DO P=1,NPART ! Loop over particles
         IF(PARTICLE(P)%TYPE.EQ.1) CYCLE ! Exclude particles in the ionized HII region

!        Calculate the absolute difference between the current and previous temperatures
         TEMPERATURE_DIFFERENCE = ABS(PARTICLE(P)%GAS_TEMPERATURE-PARTICLE(P)%PREVIOUS_TEMPERATURE)
         IF(TEMPERATURE_DIFFERENCE.LT.TEMPERATURE_DIFFERENCE_LIMIT) CYCLE

         I = I+1
         DUMMY_INDEX(I) = P
      END DO ! End of loop over particles

!  Update only the PDR particles whose temperatures have changed by more than the specified limit
   ELSE IF(TEMPERATURE_ITERATION.GT.1 .AND. CHEMISTRY_ITERATION.GT.1) THEN
      DO P=1,NPART ! Loop over particles
         IF(PARTICLE(P)%TYPE.EQ.1) CYCLE ! Exclude particles in the ionized HII region
         IF(PARTICLE(P)%TYPE.EQ.3) CYCLE ! Exclude particles in the dark cloud region

!        Calculate the absolute difference between the current and previous temperatures
         TEMPERATURE_DIFFERENCE = ABS(PARTICLE(P)%GAS_TEMPERATURE-PARTICLE(P)%PREVIOUS_TEMPERATURE)
         IF(TEMPERATURE_DIFFERENCE.LT.TEMPERATURE_DIFFERENCE_LIMIT) CYCLE

         I = I+1
         DUMMY_INDEX(I) = P
      END DO ! End of loop over particles
   END IF
   NCHEM = I ! Set the total number of particles in the list

!  Store the previously-calculated abundances for comparison
   DO P=1,NPART
      PARTICLE(P)%PREVIOUS_ABUNDANCE = PARTICLE(P)%ABUNDANCE
   END DO

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP    SHARED(NCHEM,PARTICLE,INITIAL_ABUNDANCE) &
!$OMP    SHARED(DUMMY_INDEX,DUMMY_DENSITY,DUMMY_TEMPERATURE) &
!$OMP    SHARED(DUMMY_ABUNDANCE,DUMMY_RATE) &
!$OMP    PRIVATE(P)
!  Load the relevant particle properties into dummy arrays to pass to the chemistry routine
   DO I=1,NCHEM

!     Include only the particles that will be updated on this iteration
      P = DUMMY_INDEX(I)

!     Set the abundances to their appropriate starting values
      DUMMY_ABUNDANCE(:,I) = INITIAL_ABUNDANCE

      DUMMY_RATE(:,I) = PARTICLE(P)%RATE
      DUMMY_DENSITY(I) = PARTICLE(P)%GAS_DENSITY
      DUMMY_TEMPERATURE(I) = PARTICLE(P)%GAS_TEMPERATURE
   END DO
!$OMP END PARALLEL DO

   RETURN
END SUBROUTINE UPDATE_ABUNDANCES
!=======================================================================
