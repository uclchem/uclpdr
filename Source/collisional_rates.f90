!=======================================================================
!
!  Calculate the total collisional rates (s^-1) for a given coolant at
!  the specified temperature by summing the individual rates from each
!  of the available collision partners by linear interpolation between
!  the available rate coefficients at specific temperature values.
!
!-----------------------------------------------------------------------
SUBROUTINE CALCULATE_COLLISIONAL_RATES(COOLANT,DENSITY,TEMPERATURE, &
                                     & ABUNDANCE,COLLISIONAL_RATE)

   USE HEALPIX_TYPES
   USE COOLANT_MODULE
   USE GLOBAL_MODULE, ONLY : nH,nHx,nH2,nHe,nelect
   USE FUNCTIONS_MODULE

   IMPLICIT NONE

   TYPE(COOLANT_TYPE), INTENT(IN)  :: COOLANT
   REAL(KIND=DP),      INTENT(IN)  :: DENSITY,TEMPERATURE
   REAL(KIND=DP),      INTENT(IN)  :: ABUNDANCE(:)
   REAL(KIND=DP),      INTENT(OUT) :: COLLISIONAL_RATE(:,:)

   INTEGER(KIND=I4B) :: I,J,K,KLO,KHI,PARTNER_ID
   REAL(KIND=DP) :: PARA_FRACTION,ORTHO_FRACTION
   REAL(KIND=DP) :: STEP,C_COEFF

!  Initialize the collisional rates
   COLLISIONAL_RATE=0.0D0

!  Calculate the H2 ortho/para ratio at equilibrium for the specified
!  temperature and the resulting fractions of H2 in para & ortho form
   PARA_FRACTION=1.0D0/(1.0D0+ORTHO_PARA_RATIO(TEMPERATURE))
   ORTHO_FRACTION=1.0D0-PARA_FRACTION

   DO PARTNER_ID=1,7 ! Loop over collision partners

!     Skip the collision partner if no rates are available
      IF(COOLANT%TEMPERATURE(PARTNER_ID,1).EQ.0.0D0) CYCLE

!     Determine the two nearest temperature values
!     present within the list of collisional rates
      KLO=0; KHI=0
      DO K=1,COOLANT%NTEMP ! Loop over temperatures
         IF(COOLANT%TEMPERATURE(PARTNER_ID,K).GT.TEMPERATURE) THEN
            KLO=K-1
            KHI=K
            EXIT
         ELSE IF(COOLANT%TEMPERATURE(PARTNER_ID,K).EQ.0.0D0) THEN
            KLO=K-1
            KHI=K-1
            EXIT
         END IF
      END DO

!     If the required temperature is above or below the range of available
!     temperature values then use the highest or lowest value in the range
      IF(KHI.EQ.0) THEN
         KLO=COOLANT%NTEMP
         KHI=COOLANT%NTEMP
      ELSE IF(KHI.EQ.1) THEN
         KLO=1
         KHI=1
      END IF

!     Calculate the "distance" between the two temperature
!     values, to be used in the linear interpolation below
      IF(KLO.EQ.KHI) THEN
         STEP=0.0D0
      ELSE
         STEP=(TEMPERATURE-COOLANT%TEMPERATURE(PARTNER_ID,KLO)) &
           & /(COOLANT%TEMPERATURE(PARTNER_ID,KHI)-COOLANT%TEMPERATURE(PARTNER_ID,KLO))
      END IF

!     Linearly interpolate the collisional rate coefficients
!     for each collision partner at the required temperature
      IF(PARTNER_ID.EQ.1) THEN ! Collisions with H2
         DO I=1,COOLANT%NLEVEL
            DO J=1,COOLANT%NLEVEL
               C_COEFF=COOLANT%C_COEFF(PARTNER_ID,I,J,KLO)+(COOLANT%C_COEFF(PARTNER_ID,I,J,KHI)-COOLANT%C_COEFF(PARTNER_ID,I,J,KLO))*STEP
               COLLISIONAL_RATE(I,J)=COLLISIONAL_RATE(I,J)+C_COEFF*DENSITY*ABUNDANCE(nH2)
            END DO
         END DO
      END IF
      IF(PARTNER_ID.EQ.2) THEN ! Collisions with para-H2
         DO I=1,COOLANT%NLEVEL
            DO J=1,COOLANT%NLEVEL
               C_COEFF=COOLANT%C_COEFF(PARTNER_ID,I,J,KLO)+(COOLANT%C_COEFF(PARTNER_ID,I,J,KHI)-COOLANT%C_COEFF(PARTNER_ID,I,J,KLO))*STEP
               COLLISIONAL_RATE(I,J)=COLLISIONAL_RATE(I,J)+C_COEFF*DENSITY*ABUNDANCE(nH2)*PARA_FRACTION
            END DO
         END DO
      END IF
      IF(PARTNER_ID.EQ.3) THEN ! Collisions with ortho-H2
         DO I=1,COOLANT%NLEVEL
            DO J=1,COOLANT%NLEVEL
               C_COEFF=COOLANT%C_COEFF(PARTNER_ID,I,J,KLO)+(COOLANT%C_COEFF(PARTNER_ID,I,J,KHI)-COOLANT%C_COEFF(PARTNER_ID,I,J,KLO))*STEP
               COLLISIONAL_RATE(I,J)=COLLISIONAL_RATE(I,J)+C_COEFF*DENSITY*ABUNDANCE(nH2)*ORTHO_FRACTION
            END DO
         END DO
      END IF
      IF(PARTNER_ID.EQ.4) THEN ! Collisions with electrons
         DO I=1,COOLANT%NLEVEL
            DO J=1,COOLANT%NLEVEL
               C_COEFF=COOLANT%C_COEFF(PARTNER_ID,I,J,KLO)+(COOLANT%C_COEFF(PARTNER_ID,I,J,KHI)-COOLANT%C_COEFF(PARTNER_ID,I,J,KLO))*STEP
               COLLISIONAL_RATE(I,J)=COLLISIONAL_RATE(I,J)+C_COEFF*DENSITY*ABUNDANCE(nelect)
            END DO
         END DO
      END IF
      IF(PARTNER_ID.EQ.5) THEN ! Collisions with H
         DO I=1,COOLANT%NLEVEL
            DO J=1,COOLANT%NLEVEL
               C_COEFF=COOLANT%C_COEFF(PARTNER_ID,I,J,KLO)+(COOLANT%C_COEFF(PARTNER_ID,I,J,KHI)-COOLANT%C_COEFF(PARTNER_ID,I,J,KLO))*STEP
               COLLISIONAL_RATE(I,J)=COLLISIONAL_RATE(I,J)+C_COEFF*DENSITY*ABUNDANCE(nH)
            END DO
         END DO
      END IF
      IF(PARTNER_ID.EQ.6) THEN ! Collisions with He
         DO I=1,COOLANT%NLEVEL
            DO J=1,COOLANT%NLEVEL
               C_COEFF=COOLANT%C_COEFF(PARTNER_ID,I,J,KLO)+(COOLANT%C_COEFF(PARTNER_ID,I,J,KHI)-COOLANT%C_COEFF(PARTNER_ID,I,J,KLO))*STEP
               COLLISIONAL_RATE(I,J)=COLLISIONAL_RATE(I,J)+C_COEFF*DENSITY*ABUNDANCE(nHe)
            END DO
         END DO
      END IF
      IF(PARTNER_ID.EQ.7) THEN ! Collisions with protons
         DO I=1,COOLANT%NLEVEL
            DO J=1,COOLANT%NLEVEL
               C_COEFF=COOLANT%C_COEFF(PARTNER_ID,I,J,KLO)+(COOLANT%C_COEFF(PARTNER_ID,I,J,KHI)-COOLANT%C_COEFF(PARTNER_ID,I,J,KLO))*STEP
               COLLISIONAL_RATE(I,J)=COLLISIONAL_RATE(I,J)+C_COEFF*DENSITY*ABUNDANCE(nHx)
            END DO
         END DO
      END IF

   END DO ! End of loop over collision partners

   RETURN
END SUBROUTINE CALCULATE_COLLISIONAL_RATES
!=======================================================================
