!=======================================================================
!
!  Construct the matrix of transition rates for a given coolant species
!  using the escape probability formalism to determine the local field.
!
!-----------------------------------------------------------------------
SUBROUTINE CONSTRUCT_TRANSITION_MATRIX(NRAYS,N,COOLANT,PARTICLE,TRANSITION_MATRIX)

   USE HEALPIX_TYPES
   USE COOLANT_MODULE
   USE PARTICLE_MODULE
   USE GLOBAL_MODULE, ONLY : T_CMB
   USE FUNCTIONS_MODULE
   USE SUBROUTINES_MODULE

   IMPLICIT NONE

   INTEGER(KIND=I4B),   INTENT(IN)    :: NRAYS
   INTEGER(KIND=I4B),   INTENT(IN)    :: N
   TYPE(COOLANT_TYPE),  INTENT(IN)    :: COOLANT
   TYPE(PARTICLE_TYPE), INTENT(INOUT) :: PARTICLE
   REAL(KIND=DP),       INTENT(OUT)   :: TRANSITION_MATRIX(:,:)

   INTEGER(KIND=I4B) :: I,J,K
   REAL(KIND=DP) :: TPOP,RHO_GRAIN,DUST_EMISSIVITY
   REAL(KIND=DP) :: S_ij,B_ij,B_ij_CMB,B_ij_DUST,BETA_ij
   REAL(KIND=DP), ALLOCATABLE :: RADIATION_FIELD(:,:)
   REAL(KIND=DP), ALLOCATABLE :: COLLISIONAL_RATE(:,:)

#ifdef USE_ALI
   REAL(KIND=DP) :: S_ij_PREVIOUS,LAMBDA_ij
#endif

!  Allocate and initialize the mean integrated radiation field
   ALLOCATE(RADIATION_FIELD(1:COOLANT%NLEVEL,1:COOLANT%NLEVEL))
   RADIATION_FIELD=0.0D0

!  Initialize the local emissivities
   PARTICLE%COOLANT(N)%EMISSIVITY=0.0D0

   DO I=1,COOLANT%NLEVEL
      DO J=1,COOLANT%NLEVEL
         IF(J.GE.I) EXIT

!        Calculate the background radiation field including contributions
!        from CMB blackbody emission and dust modified blackbody emission
         IF(COOLANT%FREQUENCY(I,J).LT.1.0D15) THEN
            B_ij_CMB=(2*HP*COOLANT%FREQUENCY(I,J)**3)/(C**2)/(EXP(HP*COOLANT%FREQUENCY(I,J)/(KB*T_CMB))-1.0D0)
#ifdef DUST
            RHO_GRAIN=2.0D0 ! Grain mass density (g cm^-3)
            DUST_EMISSIVITY=(RHO_GRAIN*PARTICLE%DUST_DENSITY)*(0.01*(1.3*COOLANT%FREQUENCY(I,J)/3.0D11))
            B_ij_DUST=(2*HP*COOLANT%FREQUENCY(I,J)**3)/(C**2)/(EXP(HP*COOLANT%FREQUENCY(I,J)/(KB*PARTICLE%DUST_TEMPERATURE))-1.D0)*DUST_EMISSIVITY
#else
            B_ij_DUST=0.0D0
#endif
            B_ij=B_ij_CMB+B_ij_DUST
         ELSE
            B_ij=0.0D0
         END IF

!        If the population n_i is zero then the source function is zero and the
!        mean integrated radiation field is just the background radiation field
         IF(PARTICLE%COOLANT(N)%POPULATION(I).EQ.0) THEN
            S_ij=0.0D0
            BETA_ij=1.0D0

!        If the difference between n_i.g_j and n_j.g_i is vanishingly small, then
!        calculate the source function directly and set the escape probability to 1
         ELSE IF(ABS(PARTICLE%COOLANT(N)%POPULATION(I)*COOLANT%WEIGHT(J)-PARTICLE%COOLANT(N)%POPULATION(J)*COOLANT%WEIGHT(I)).EQ.0) THEN
            S_ij=HP*COOLANT%FREQUENCY(I,J)*PARTICLE%COOLANT(N)%POPULATION(I)*COOLANT%A_COEFF(I,J)/(4.0*PI)
            BETA_ij=1.0D0

         ELSE
!           Calculate the source function
            S_ij=(2*HP*COOLANT%FREQUENCY(I,J)**3)/(C**2) &
              & /((PARTICLE%COOLANT(N)%POPULATION(J)*COOLANT%WEIGHT(I)) &
              &  /(PARTICLE%COOLANT(N)%POPULATION(I)*COOLANT%WEIGHT(J))-1.0D0)

!           Calculate the escape probability
            BETA_ij=ESCAPE_PROBABILITY(NRAYS,PARTICLE%COOLANT(N)%OPACITY(0:NRAYS-1,I,J))

         END IF

!        Calculate the local emissivity (erg cm^-3 s^-1) for the transition line
         IF(PARTICLE%COOLANT(N)%POPULATION(I).GT.0) THEN
            PARTICLE%COOLANT(N)%EMISSIVITY(I,J)=PARTICLE%COOLANT(N)%POPULATION(I)*COOLANT%A_COEFF(I,J) &
                                             & *HP*COOLANT%FREQUENCY(I,J)*BETA_ij*(S_ij-B_ij)/S_ij
!!$            IF(N.EQ.1) WRITE(90,"(2I3,10ES12.5)") I,J,PARTICLE%AV(0),PARTICLE%COOLANT(N)%POPULATION(I),BETA_ij,S_ij,B_ij,PARTICLE%COOLANT(N)%EMISSIVITY(I,J)
         END IF

!        Calculate the mean integrated radiation field <J_ij>
         RADIATION_FIELD(I,J)=(1.0D0-BETA_ij)*S_ij+BETA_ij*B_ij
         RADIATION_FIELD(J,I)=RADIATION_FIELD(I,J)

#ifdef USE_ALI
!        Calculate the source function in the same manner for
!        the populations determined on the previous iteration
         IF(PARTICLE%COOLANT(N)%PREVIOUS_POPULATION(I).EQ.0) THEN
            S_ij_PREVIOUS=0.0D0
         ELSE IF(ABS(PARTICLE%COOLANT(N)%PREVIOUS_POPULATION(I)*COOLANT%WEIGHT(J)-PARTICLE%COOLANT(N)%PREVIOUS_POPULATION(J)*COOLANT%WEIGHT(I)).EQ.0) THEN
            S_ij_PREVIOUS=HP*COOLANT%FREQUENCY(I,J)*PARTICLE%COOLANT(N)%PREVIOUS_POPULATION(I)*COOLANT%A_COEFF(I,J)/(4.0*PI)
         ELSE
            S_ij_PREVIOUS=(2*HP*COOLANT%FREQUENCY(I,J)**3)/(C**2) &
              & /((PARTICLE%COOLANT(N)%PREVIOUS_POPULATION(J)*COOLANT%WEIGHT(I)) &
              &  /(PARTICLE%COOLANT(N)%PREVIOUS_POPULATION(I)*COOLANT%WEIGHT(J))-1.0D0)
         END IF

!        Use the Accelerated Lambda Iteration method to speed up convergence
!        by amplifying the incremental difference of the new source function
         LAMBDA_ij=PARTICLE%COOLANT(N)%LAMBDA(I,J)
         RADIATION_FIELD(I,J)=(1.0D0-BETA_ij)*(LAMBDA_ij*S_ij+(1.0D0-LAMBDA_ij)*S_ij_PREVIOUS)+BETA_ij*B_ij
         RADIATION_FIELD(J,I)=RADIATION_FIELD(I,J)
#endif

      END DO ! End of loop over levels (j)
   END DO ! End of loop over levels (i)

!  Allocate and calculate the collisional rates
   ALLOCATE(COLLISIONAL_RATE(1:COOLANT%NLEVEL,1:COOLANT%NLEVEL))
   CALL CALCULATE_COLLISIONAL_RATES(COOLANT,PARTICLE%GAS_DENSITY,PARTICLE%GAS_TEMPERATURE,PARTICLE%ABUNDANCE,COLLISIONAL_RATE)

!  Construct the transition matrix: R_ij = A_ij + B_ij.<J> + C_ij
   DO I=1,COOLANT%NLEVEL
      DO J=1,COOLANT%NLEVEL
         TRANSITION_MATRIX(I,J)=COOLANT%A_COEFF(I,J) &
                            & + COOLANT%B_COEFF(I,J)*RADIATION_FIELD(I,J) &
                            & + COLLISIONAL_RATE(I,J)
      END DO
   END DO

   DEALLOCATE(RADIATION_FIELD,COLLISIONAL_RATE)

   RETURN
END SUBROUTINE CONSTRUCT_TRANSITION_MATRIX
!=======================================================================
