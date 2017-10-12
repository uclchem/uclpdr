!=======================================================================
!
!  Update the gas temperature based on the heating-cooling imbalance.
!
!-----------------------------------------------------------------------
SUBROUTINE UPDATE_TEMPERATURE(PARTICLE)

   USE HEALPIX_TYPES
   USE PARTICLE_MODULE
   USE MAIN_MODULE, ONLY : FCRIT,TDIFF,TMIN,TMAX

   IMPLICIT NONE

   TYPE(PARTICLE_TYPE), INTENT(INOUT) :: PARTICLE

   REAL(KIND=DP) :: DUMMY_TEMPERATURE,TEMPERATURE_DIFFERENCE
   REAL(KIND=DP) :: DIFFERENCE,RELATIVE_DIFFERENCE
   REAL(KIND=DP) :: TOTAL_COOLING,TOTAL_HEATING

!  Reset the temperature convergence flag
   PARTICLE%TEMPERATURE_CONVERGED=.FALSE.

!  Store the current temperature in a dummy variable
   DUMMY_TEMPERATURE=PARTICLE%GAS_TEMPERATURE

!  Calculate the absolute difference between the current and previous temperatures
   TEMPERATURE_DIFFERENCE=ABS(PARTICLE%GAS_TEMPERATURE-PARTICLE%PREVIOUS_TEMPERATURE)

!  Calculate the total cooling and heating rates
   TOTAL_COOLING=SUM(PARTICLE%COOLING_RATE)
   TOTAL_HEATING=SUM(PARTICLE%HEATING_RATE)

!  Calculate the difference between the total heating and total cooling rates
!  and the absolute value of the relative difference between the two rates
   DIFFERENCE=TOTAL_HEATING-TOTAL_COOLING
   RELATIVE_DIFFERENCE=2.0D0*ABS(DIFFERENCE)/ABS(TOTAL_HEATING+TOTAL_COOLING)

!  Check for convergence in the heating-cooling imbalance
   IF(RELATIVE_DIFFERENCE.LE.FCRIT) THEN
      PARTICLE%PREVIOUS_TEMPERATURE=DUMMY_TEMPERATURE
      PARTICLE%GAS_TEMPERATURE=DUMMY_TEMPERATURE
      PARTICLE%TEMPERATURE_CONVERGED=.TRUE.
      RETURN
   END IF

!  Determine the temperature bracket to begin searching within by first increasing
!  or decreasing the temperature by 30% according to the heating-cooling imbalance
   IF(.NOT.PARTICLE%BINARY_CHOP_SEARCH) THEN

!     If the heating continues to outweigh the cooling, increase the temperature by 30%
      IF(DIFFERENCE.GT.0 .AND. PARTICLE%PREVIOUS_DIFFERENCE.GE.0) THEN
         PARTICLE%GAS_TEMPERATURE=1.3D0*PARTICLE%GAS_TEMPERATURE
         PARTICLE%PREVIOUS_DIFFERENCE=DIFFERENCE
         PARTICLE%TLOW=DUMMY_TEMPERATURE ! Update the value of T_low
         PARTICLE%THIGH=PARTICLE%GAS_TEMPERATURE ! Update the value of T_high

!     If the cooling continues to outweigh the heating, decrease the temperature by 30%
      ELSE IF(DIFFERENCE.LT.0 .AND. PARTICLE%PREVIOUS_DIFFERENCE.LE.0) THEN
         PARTICLE%GAS_TEMPERATURE=0.7D0*PARTICLE%GAS_TEMPERATURE
         PARTICLE%PREVIOUS_DIFFERENCE=DIFFERENCE
         PARTICLE%TLOW=PARTICLE%GAS_TEMPERATURE ! Update the value of T_low
         PARTICLE%THIGH=DUMMY_TEMPERATURE ! Update the value of T_high

!     If the heating-cooling balance has reversed (either from net heating to net cooling or
!     vice-versa) then switch to the binary chop search method to determine the temperature
      ELSE
         PARTICLE%GAS_TEMPERATURE=(PARTICLE%THIGH+PARTICLE%TLOW)/2.0D0
         PARTICLE%PREVIOUS_DIFFERENCE=DIFFERENCE
         PARTICLE%BINARY_CHOP_SEARCH=.TRUE. ! From now on
      END IF

!  Perform a binary chop search (the min-max range was found by the 30% increase/decrease method)
   ELSE

      IF(DIFFERENCE.GT.0) THEN
         PARTICLE%GAS_TEMPERATURE=(PARTICLE%GAS_TEMPERATURE+PARTICLE%THIGH)/2.0D0
         PARTICLE%PREVIOUS_DIFFERENCE=DIFFERENCE
         PARTICLE%TLOW=DUMMY_TEMPERATURE ! Update the value of T_low
      END IF
      IF(DIFFERENCE.LT.0) THEN
         PARTICLE%GAS_TEMPERATURE=(PARTICLE%GAS_TEMPERATURE+PARTICLE%TLOW)/2.0D0
         PARTICLE%PREVIOUS_DIFFERENCE=DIFFERENCE
         PARTICLE%THIGH=DUMMY_TEMPERATURE ! Update the value of T_high
      END IF

   END IF

!  If the search routine is unable to converge on a temperature that satisfies the thermal balance
!  criterion, expand the min-max search bracket asymmetrically and begin to narrow the search again
!  If the repeated search fails to converge once more, force convergence at the current temperature
   IF(TEMPERATURE_DIFFERENCE.LE.TDIFF) THEN
      IF(.NOT.PARTICLE%BRACKET_EXPANDED) THEN
         PARTICLE%THIGH=PARTICLE%THIGH+SQRT(PI)
         PARTICLE%TLOW=PARTICLE%TLOW-SQRT(2.0)
         PARTICLE%BRACKET_EXPANDED=.TRUE.
      ELSE
         PARTICLE%PREVIOUS_TEMPERATURE=DUMMY_TEMPERATURE
         PARTICLE%GAS_TEMPERATURE=DUMMY_TEMPERATURE
         PARTICLE%TEMPERATURE_CONVERGED=.TRUE.
         RETURN
      END IF
   END IF

!  Check if the temperature falls outside of the allowed limits and force convergence if so
   IF(PARTICLE%GAS_TEMPERATURE.LE.TMIN .AND. DIFFERENCE.LT.0) THEN
      PARTICLE%GAS_TEMPERATURE=TMIN
      PARTICLE%TEMPERATURE_CONVERGED=.TRUE.
   END IF
   IF(PARTICLE%GAS_TEMPERATURE.GE.TMAX .AND. DIFFERENCE.GT.0) THEN
      PARTICLE%GAS_TEMPERATURE=TMAX
      PARTICLE%TEMPERATURE_CONVERGED=.TRUE.
   END IF

!  Replace the previous temperature with the current value
   PARTICLE%PREVIOUS_TEMPERATURE=DUMMY_TEMPERATURE

   RETURN
END SUBROUTINE UPDATE_TEMPERATURE
!=======================================================================
