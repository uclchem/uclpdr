!=======================================================================
!
!  Update the abundances for all coolant species, accounting for the
!  separate para- and ortho- forms of H2, if present, by determining
!  the ortho/para ratio at the required temperature.
!
!=======================================================================
SUBROUTINE UPDATE_COOLANT_ABUNDANCES(NPART,NCOOL,COOLANT,PARTICLE)

   USE HEALPIX_TYPES
   USE COOLANT_MODULE
   USE PARTICLE_MODULE
   USE FUNCTIONS_MODULE

   IMPLICIT NONE

   INTEGER(KIND=I4B),   INTENT(IN)    :: NPART
   INTEGER(KIND=I4B),   INTENT(IN)    :: NCOOL
   TYPE(COOLANT_TYPE),  INTENT(IN)    :: COOLANT(*)
   TYPE(PARTICLE_TYPE), INTENT(INOUT) :: PARTICLE(*)

   INTEGER(KIND=I4B) :: N,P
   REAL(KIND=DP) :: PARA_FRACTION,ORTHO_FRACTION

   DO N=1,NCOOL ! Loop over coolants
      DO P=1,NPART ! Loop over particles
         IF(COOLANT(N)%NAME.EQ."p-H2") THEN
!           Calculate the equilibrium H2 ortho/para ratio at the particle
!           gas temperature and the resulting fraction of H2 in para form
            PARA_FRACTION=1.0D0/(1.0D0+ORTHO_PARA_RATIO(PARTICLE(P)%GAS_TEMPERATURE))
            PARTICLE(P)%COOLANT(N)%DENSITY=PARTICLE(P)%ABUNDANCE(COOLANT(N)%INDEX) &
                                        & *PARA_FRACTION*PARTICLE(P)%GAS_DENSITY
         ELSE IF(COOLANT(N)%NAME.EQ."o-H2") THEN
!           Calculate the equilibrium H2 ortho/para ratio at the particle
!           gas temperature and the resulting fraction of H2 in ortho form
            ORTHO_FRACTION=1.0D0/(1.0D0+1.0D0/ORTHO_PARA_RATIO(PARTICLE(P)%GAS_TEMPERATURE))
            PARTICLE(P)%COOLANT(N)%DENSITY=PARTICLE(P)%ABUNDANCE(COOLANT(N)%INDEX) &
                                        & *ORTHO_FRACTION*PARTICLE(P)%GAS_DENSITY
         ELSE
            PARTICLE(P)%COOLANT(N)%DENSITY=PARTICLE(P)%ABUNDANCE(COOLANT(N)%INDEX) &
                                        & *PARTICLE(P)%GAS_DENSITY
         END IF
      END DO ! End of loop over particles
   END DO ! End of loop over coolants

   RETURN
END SUBROUTINE UPDATE_COOLANT_ABUNDANCES
!=======================================================================
