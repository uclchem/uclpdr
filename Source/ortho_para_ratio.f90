!=======================================================================
!
!  Calculate the ortho-to-para ratio of H2 at thermal equilibrium
!  for the specified temperature, making use of energy level data
!  if available, or an approximation if not.
!
!-----------------------------------------------------------------------
FUNCTION ORTHO_PARA_RATIO(TEMPERATURE)

   USE HEALPIX_TYPES
   USE MAIN_MODULE

   IMPLICIT NONE

   REAL(KIND=DP) :: ORTHO_PARA_RATIO
   REAL(KIND=DP), INTENT(IN) :: TEMPERATURE

   INTEGER(KIND=I4B) :: I,J,N,ORTHO_INDEX,PARA_INDEX
   REAL(KIND=DP) :: I_ORTHO,I_PARA,ORTHO_FRACTION,PARA_FRACTION

!  Check if coolant data is available for the ortho and para forms
   ORTHO_INDEX=0; PARA_INDEX=0
   DO N=1,NCOOL
      IF(COOLANT(N)%NAME.EQ."o-H2") ORTHO_INDEX=N
      IF(COOLANT(N)%NAME.EQ."p-H2") PARA_INDEX=N
   END DO

!  Calculate the exact ortho/para ratio if molecular data is available
   IF(ORTHO_INDEX.NE.0 .AND. PARA_INDEX.NE.0) THEN

      I_ORTHO=1.0D0; I_PARA=0.0D0 ! Total nuclear spins of the two forms

!     Calculate the ortho/para ratio of H2 using the expression
!     from Poelman & Spaans (2005, A&A, 440, 559, equation 11)
      ORTHO_FRACTION=0.0D0
      DO I=1,COOLANT(ORTHO_INDEX)%NLEVEL
         ORTHO_FRACTION=ORTHO_FRACTION+COOLANT(ORTHO_INDEX)%WEIGHT(I) &
                     & *EXP(-COOLANT(ORTHO_INDEX)%ENERGY(I)/(KB*TEMPERATURE))
      END DO

      PARA_FRACTION=0.0D0
      DO I=1,COOLANT(PARA_INDEX)%NLEVEL
         PARA_FRACTION=PARA_FRACTION+COOLANT(PARA_INDEX)%WEIGHT(I) &
                    & *EXP(-COOLANT(PARA_INDEX)%ENERGY(I)/(KB*TEMPERATURE))
      END DO

      ORTHO_PARA_RATIO=(2*I_ORTHO+1)/(2*I_PARA+1)*(ORTHO_FRACTION/PARA_FRACTION)

   ELSE

!     Approximate the ortho/para ratio of H2 using the expression
!     given by Flower & Watt (1985, MNRAS, 213, 991, equation 2)
      ORTHO_PARA_RATIO=9.0D0*EXP(-170.5D0/TEMPERATURE)

!     Limit the ortho/para ratio to its statistical limit
      IF(ORTHO_PARA_RATIO.GT.3.0D0) ORTHO_PARA_RATIO=3.0D0

   END IF

END FUNCTION ORTHO_PARA_RATIO
!=======================================================================
