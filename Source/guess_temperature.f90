!=======================================================================
!
!  Calculate initial guesses for the gas temperature of each particle
!  by using an empirical approximation derived from fits to a grid of
!  1-dimensional model results. If specified, the dust temperature is
!  also set to a fixed user-supplied value.
!
!-----------------------------------------------------------------------
SUBROUTINE CALCULATE_GUESS_TEMPERATURES(NPART,PARTICLE)

   USE HEALPIX_TYPES
   USE PARTICLE_MODULE
   USE MAIN_MODULE, ONLY : GAS_TEMPERATURE_GUESS,DUST_TEMPERATURE_GUESS,TMIN,TMAX

   IMPLICIT NONE

   INTEGER(KIND=I4B),   INTENT(IN)    :: NPART
   TYPE(PARTICLE_TYPE), INTENT(INOUT) :: PARTICLE(*)

   INTEGER(KIND=I4B) :: P

   DO P=1,NPART ! Loop over particles

#ifdef GUESS_TEMP
      PARTICLE(P)%GAS_TEMPERATURE=10.0D0*(1+(100*PARTICLE(P)%FUV_FLUX)**(1.0D0/3.0D0)+(10000*PARTICLE(P)%XRAY_FLUX)**(1.0D0/2.0D0))
#else
      IF(GAS_TEMPERATURE_GUESS.GE.TMIN .AND. GAS_TEMPERATURE_GUESS.LE.TMAX) THEN
         PARTICLE(P)%GAS_TEMPERATURE=GAS_TEMPERATURE_GUESS
      ELSE
         WRITE(6,*) 'ERROR! User-supplied initial guess for the gas temperature is outside the allowed range specified by TMIN and TMAX'
         STOP
      END IF
#endif

#ifndef CALC_TDUST
      IF(DUST_TEMPERATURE_GUESS.GE.TMIN .AND. DUST_TEMPERATURE_GUESS.LE.TMAX) THEN
         PARTICLE(P)%DUST_TEMPERATURE=DUST_TEMPERATURE_GUESS
      ELSE
         WRITE(6,*) 'ERROR! User-supplied fixed dust temperature is outside the allowed range specified by TMIN and TMAX'
         STOP
      END IF
#endif

   END DO ! End of loop over particles

   RETURN
END SUBROUTINE CALCULATE_GUESS_TEMPERATURES
!=======================================================================
