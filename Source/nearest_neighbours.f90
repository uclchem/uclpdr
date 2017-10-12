!=======================================================================
!
!  Calculate the population densities of all levels of a given coolant
!  for the specified particle by taking the nearest neighbour weighted
!  interpolation of the converged population densities of neighbouring
!  particles.
!
!-----------------------------------------------------------------------
SUBROUTINE USE_NEARBY_POPULATIONS(NPART,PARTICLE_ID,COOLANT_ID,PARTICLE)

   USE HEALPIX_TYPES
   USE PARTICLE_MODULE

   INTEGER(KIND=I4B),   INTENT(IN)    :: NPART,PARTICLE_ID,COOLANT_ID
   TYPE(PARTICLE_TYPE), INTENT(INOUT) :: PARTICLE(*)

   INTEGER(KIND=I4B) :: I,L,M,N,P,NPOINT
   REAL(KIND=DP) :: STEP_SIZE,STEP_TOTAL

   L = 8 ! Distance weighting factor (a.k.a. power parameter)

   P = PARTICLE_ID
   N = COOLANT_ID

!  Store the previously-calculated population densities for comparison
   PARTICLE(P)%COOLANT(N)%PREVIOUS_POPULATION = PARTICLE(P)%COOLANT(N)%POPULATION
   PARTICLE(P)%COOLANT(N)%POPULATION = 0.0D0

   NPOINT = 0
   STEP_TOTAL = 0.0D0

   DO M=1,NPART ! Loop over particles

      IF(.NOT.PARTICLE(M)%COOLANT(N)%CONVERGED) CYCLE ! Skip this particle if its population densities have not converged

!     Calculate the distance of this particle from the particle of interest
      STEP_SIZE = SQRT((PARTICLE(M)%COORDINATES(1)-PARTICLE(P)%COORDINATES(1))**2 &
                   & + (PARTICLE(M)%COORDINATES(2)-PARTICLE(P)%COORDINATES(2))**2 &
                   & + (PARTICLE(M)%COORDINATES(3)-PARTICLE(P)%COORDINATES(3))**2)
      IF(STEP_SIZE.EQ.0) CYCLE

!     Add the contribution from this particle, inversely weighted by its distance from the particle of interest
      PARTICLE(P)%COOLANT(N)%POPULATION = PARTICLE(P)%COOLANT(N)%POPULATION + PARTICLE(M)%COOLANT(N)%POPULATION/STEP_SIZE**L
      STEP_TOTAL = STEP_TOTAL + 1.0D0/STEP_SIZE**L
      NPOINT = NPOINT + 1

   END DO ! End of loop over particles

!  Calculate the population densities by inverse distance weighted (IDW) interpolation
   IF(NPOINT.GT.0) THEN
      PARTICLE(P)%COOLANT(N)%POPULATION = PARTICLE(P)%COOLANT(N)%POPULATION/STEP_TOTAL
   END IF

   RETURN
END SUBROUTINE USE_NEARBY_POPULATIONS
!=======================================================================

!=======================================================================
!
!  Calculate new temperatures for particles that have not achieved the
!  required thermal balance criterion using nearest neighbour weighted
!  interpolation of the temperatures from nearby converged particles.
!
!-----------------------------------------------------------------------
SUBROUTINE USE_NEARBY_TEMPERATURES(NPART,PARTICLE)

   USE HEALPIX_TYPES
   USE PARTICLE_MODULE
   USE MAIN_MODULE, ONLY : FCRIT,TMIN,TMAX

   INTEGER(KIND=I4B),   INTENT(IN)    :: NPART
   TYPE(PARTICLE_TYPE), INTENT(INOUT) :: PARTICLE(*)

   INTEGER(KIND=I4B) :: I,L,M,N,P,NPOINT
   REAL(KIND=DP) :: STEP_SIZE,STEP_TOTAL
   REAL(KIND=DP) :: RELATIVE_DIFFERENCE

   L = 8 ! Distance weighting factor (a.k.a. power parameter)

   DO P=1,NPART ! Loop over particles

!     Calculate the relative difference between the total heating and total cooling rates
      RELATIVE_DIFFERENCE=2.0D0*ABS(SUM(PARTICLE(P)%HEATING_RATE)-SUM(PARTICLE(P)%COOLING_RATE))/ &
                              & ABS(SUM(PARTICLE(P)%HEATING_RATE)+SUM(PARTICLE(P)%COOLING_RATE))

!     Skip particles that have reached the thermal balance convergence criterion
!     or that have temperatures that fall outside the allowed min/max limits
      IF(RELATIVE_DIFFERENCE.LE.FCRIT .OR. &
       & PARTICLE(P)%GAS_TEMPERATURE.LE.TMIN .OR. &
       & PARTICLE(P)%GAS_TEMPERATURE.GE.TMAX) CYCLE

      PARTICLE(P)%GAS_TEMPERATURE = 0.0D0
      NPOINT = 0
      STEP_TOTAL = 0.0D0

      DO M=1,NPART ! Loop over particles

         IF(PARTICLE(M)%PREVIOUS_DIFFERENCE.GT.FCRIT) CYCLE ! Skip this particle if its temperature has not converged

!        Calculate the distance of this particle from the particle of interest
         STEP_SIZE = SQRT((PARTICLE(M)%COORDINATES(1)-PARTICLE(P)%COORDINATES(1))**2 &
                      & + (PARTICLE(M)%COORDINATES(2)-PARTICLE(P)%COORDINATES(2))**2 &
                      & + (PARTICLE(M)%COORDINATES(3)-PARTICLE(P)%COORDINATES(3))**2)
         IF(STEP_SIZE.EQ.0) CYCLE

!        Add the contribution from this particle, inversely weighted by its distance from the particle of interest
         PARTICLE(P)%GAS_TEMPERATURE = PARTICLE(P)%GAS_TEMPERATURE + PARTICLE(M)%GAS_TEMPERATURE/STEP_SIZE**L
         STEP_TOTAL = STEP_TOTAL + 1.0D0/STEP_SIZE**L
         NPOINT = NPOINT + 1

      END DO ! End of loop over particles

!     Calculate the temperature by inverse distance weighted (IDW) interpolation
      IF(NPOINT.GT.0) THEN
         PARTICLE(P)%GAS_TEMPERATURE = PARTICLE(P)%GAS_TEMPERATURE/STEP_TOTAL
      END IF

!     Set the previously-calculated gas temperature to the same value
      PARTICLE(P)%PREVIOUS_TEMPERATURE = PARTICLE(P)%GAS_TEMPERATURE

   END DO ! End of loop over particles

   RETURN
END SUBROUTINE USE_NEARBY_TEMPERATURES
!=======================================================================
