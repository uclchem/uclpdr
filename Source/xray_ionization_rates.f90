!=======================================================================
!
!  Calculate the rate coefficients for all X-ray primary ionization
!  reactions for the specified incident X-ray flux.
!
!=======================================================================
SUBROUTINE CALCULATE_XRAY_IONIZATION_RATES(ENERGY,FLUX,RATE,REACTANT,NREAC)

   USE HEALPIX_TYPES
   USE XRAY_FUNCTIONS
   USE XRAY_IONIZATION_MODULE

   IMPLICIT NONE

   INTEGER(KIND=I4B), INTENT(IN)    :: NREAC
   CHARACTER(LEN=10), INTENT(IN)    :: REACTANT(:,:)
   REAL(KIND=DP),     INTENT(IN)    :: ENERGY(:),FLUX(:)
   REAL(KIND=DP),     INTENT(INOUT) :: RATE(:)

   INTEGER(KIND=I4B) :: I,NPOINT

   NPOINT = SIZE(ENERGY)

!  Loop over all reactions
   DO I=1,NREAC
      IF(REACTANT(I,2).EQ."XRAY ") THEN
         RATE(I) = XRAY_IONIZATION_RATE(ENERGY,FLUX,REACTANT(I,1),NPOINT)
      END IF
   END DO

END SUBROUTINE CALCULATE_XRAY_IONIZATION_RATES
!=======================================================================

!=======================================================================
!
!  Calculate the photoionization rate for the required species given
!  the incident X-ray spectrum.
!
!-----------------------------------------------------------------------
FUNCTION XRAY_IONIZATION_RATE(ENERGY,FLUX,SPECIES,NPOINT) RESULT(RATE)

   USE HEALPIX_TYPES
   USE XRAY_FUNCTIONS

   IMPLICIT NONE

   INTEGER(KIND=I4B), INTENT(IN) :: NPOINT
   CHARACTER(LEN=10), INTENT(IN) :: SPECIES
   REAL(KIND=DP),     INTENT(IN) :: FLUX(1:NPOINT)
   REAL(KIND=DP),     INTENT(IN) :: ENERGY(1:NPOINT)

   REAL(KIND=DP) :: RATE
   REAL(KIND=DP) :: SIGMA(1:NPOINT)

   IF(SPECIES.EQ."H  ") THEN
      SIGMA = HYDROGEN_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)
   ELSE IF(SPECIES.EQ."H+ ") THEN
      SIGMA = HYDROGEN_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)
   ELSE IF(SPECIES.EQ."He ") THEN
      SIGMA = HELIUM_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)
   ELSE IF(SPECIES.EQ."He+") THEN
      SIGMA = HELIUM_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)
   ELSE IF(SPECIES.EQ."C  ") THEN
      SIGMA = CARBON_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)
   ELSE IF(SPECIES.EQ."C+ ") THEN
      SIGMA = CARBON_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)
   ELSE IF(SPECIES.EQ."N  ") THEN
      SIGMA = NITROGEN_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)
   ELSE IF(SPECIES.EQ."N+ ") THEN
      SIGMA = NITROGEN_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)
   ELSE IF(SPECIES.EQ."O  ") THEN
      SIGMA = OXYGEN_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)
   ELSE IF(SPECIES.EQ."O+ ") THEN
      SIGMA = OXYGEN_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)
   ELSE IF(SPECIES.EQ."Ne ") THEN
      SIGMA = NEON_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)
   ELSE IF(SPECIES.EQ."Ne+") THEN
      SIGMA = NEON_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)
   ELSE IF(SPECIES.EQ."Na ") THEN
      SIGMA = SODIUM_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)
   ELSE IF(SPECIES.EQ."Na+") THEN
      SIGMA = SODIUM_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)
   ELSE IF(SPECIES.EQ."Mg ") THEN
      SIGMA = MAGNESIUM_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)
   ELSE IF(SPECIES.EQ."Mg+") THEN
      SIGMA = MAGNESIUM_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)
   ELSE IF(SPECIES.EQ."Al ") THEN
      SIGMA = ALUMINIUM_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)
   ELSE IF(SPECIES.EQ."Al+") THEN
      SIGMA = ALUMINIUM_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)
   ELSE IF(SPECIES.EQ."Si ") THEN
      SIGMA = SILICON_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)
   ELSE IF(SPECIES.EQ."Si+") THEN
      SIGMA = SILICON_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)
   ELSE IF(SPECIES.EQ."P  ") THEN
      SIGMA = PHOSPHORUS_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)
   ELSE IF(SPECIES.EQ."P+ ") THEN
      SIGMA = PHOSPHORUS_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)
   ELSE IF(SPECIES.EQ."S  ") THEN
      SIGMA = SULPHUR_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)
   ELSE IF(SPECIES.EQ."S+ ") THEN
      SIGMA = SULPHUR_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)
   ELSE IF(SPECIES.EQ."Cl ") THEN
      SIGMA = CHLORINE_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)
   ELSE IF(SPECIES.EQ."Cl+") THEN
      SIGMA = CHLORINE_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)
   ELSE IF(SPECIES.EQ."Ar ") THEN
      SIGMA = ARGON_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)
   ELSE IF(SPECIES.EQ."Ar+") THEN
      SIGMA = ARGON_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)
   ELSE IF(SPECIES.EQ."Ca ") THEN
      SIGMA = CALCIUM_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)
   ELSE IF(SPECIES.EQ."Ca+") THEN
      SIGMA = CALCIUM_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)
   ELSE IF(SPECIES.EQ."Cr ") THEN
      SIGMA = CHROMIUM_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)
   ELSE IF(SPECIES.EQ."Cr+") THEN
      SIGMA = CHROMIUM_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)
   ELSE IF(SPECIES.EQ."Fe ") THEN
      SIGMA = IRON_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)
   ELSE IF(SPECIES.EQ."Fe+") THEN
      SIGMA = IRON_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)
   ELSE IF(SPECIES.EQ."Ni ") THEN
      SIGMA = NICKEL_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)
   ELSE IF(SPECIES.EQ."Ni+") THEN
      SIGMA = NICKEL_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)
   ELSE
      WRITE(6,*) 'ERROR! No X-ray photoionization cross section data available for species ',TRIM(SPECIES)
      STOP
   END IF

   RATE = INTEGRATE(SIGMA*FLUX/ENERGY,ENERGY,NPOINT)

   RETURN
END FUNCTION XRAY_IONIZATION_RATE
!=======================================================================
