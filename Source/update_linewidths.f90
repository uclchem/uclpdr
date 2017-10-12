!=======================================================================
!
!  Calculate the Doppler line width of each coolant species for all
!  particles. Contributions from both thermal and turbulent motions
!  are included.
!
!-----------------------------------------------------------------------
SUBROUTINE UPDATE_COOLANT_LINEWIDTHS(NPART,NCOOL,COOLANT,PARTICLE)

   USE HEALPIX_TYPES
   USE COOLANT_MODULE
   USE PARTICLE_MODULE
   USE GLOBAL_MODULE, ONLY : V_TURB

   IMPLICIT NONE

   INTEGER(KIND=I4B),   INTENT(IN)    :: NPART,NCOOL
   TYPE(COOLANT_TYPE),  INTENT(IN)    :: COOLANT(*)
   TYPE(PARTICLE_TYPE), INTENT(INOUT) :: PARTICLE(*)

   INTEGER(KIND=I4B) :: N,P
   REAL(KIND=DP) :: V_THERM

!  Calculate the mean thermal velocity of each coolant species at the relevant
!  gas temperature for each particle. Update the Doppler line widths (cm s^-1)
!  by adding the thermal and turbulent velocities in quadrature.
   DO N=1,NCOOL ! Loop over coolants
      DO P=1,NPART ! Loop over particles
         V_THERM = SQRT(2*KB*PARTICLE(P)%GAS_TEMPERATURE/(COOLANT(N)%MOLECULAR_MASS*MP)) ! v_thermal = (2kT/m)^1/2
         PARTICLE(P)%COOLANT(N)%LINEWIDTH = SQRT(V_THERM**2 + V_TURB**2) ! δv_D = (v_thermal^2 + v_turbulent^2)^1/2
      END DO ! End of loop over particles
   END DO ! End of loop over coolants

   RETURN
END SUBROUTINE UPDATE_COOLANT_LINEWIDTHS
!=======================================================================
