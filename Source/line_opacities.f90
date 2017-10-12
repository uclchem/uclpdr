!=======================================================================
!
!  Calculate the line opacity of each coolant transition along each
!  HEALPix ray to the PDR surface (or simulation boundary) for all
!  particles. Calculations can be performed using either the Euler
!  or Trapezium integration scheme (specified by a compiler flag).
!
!-----------------------------------------------------------------------
SUBROUTINE CALCULATE_LINE_OPACITIES(NPART,NCOOL,NRAYS,COOLANT,PARTICLE)

   USE HEALPIX_TYPES
   USE COOLANT_MODULE
   USE PARTICLE_MODULE

   IMPLICIT NONE

   INTEGER(KIND=I4B),   INTENT(IN)    :: NPART,NCOOL,NRAYS
   TYPE(COOLANT_TYPE),  INTENT(IN)    :: COOLANT(*)
   TYPE(PARTICLE_TYPE), INTENT(INOUT) :: PARTICLE(*)

   INTEGER(KIND=I4B) :: I,J,L,M,N,P,ILEVEL,JLEVEL
   REAL(KIND=DP) :: STEP_SIZE,FACTOR1,FACTOR2,FACTOR3

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(DYNAMIC) &
!$OMP    SHARED(NPART,NCOOL,NRAYS,COOLANT,PARTICLE) &
!$OMP    PRIVATE(I,J,L,M,N,ILEVEL,JLEVEL,STEP_SIZE) &
!$OMP    PRIVATE(FACTOR1,FACTOR2,FACTOR3)
   DO P=1,NPART ! Loop over particles
      DO N=1,NCOOL ! Loop over coolants
         IF(COOLANT(N)%CONVERGED) CYCLE
         PARTICLE(P)%COOLANT(N)%OPACITY = 0.0D0
         DO ILEVEL=1,COOLANT(N)%NLEVEL ! Loop over levels (i)
            DO JLEVEL=1,COOLANT(N)%NLEVEL ! Loop over levels (j)
               IF(COOLANT(N)%A_COEFF(ILEVEL,JLEVEL).EQ.0) CYCLE
               FACTOR1 = (COOLANT(N)%A_COEFF(ILEVEL,JLEVEL)*C**3)/(8*PI*COOLANT(N)%FREQUENCY(ILEVEL,JLEVEL)**3) ! = A_ij.c^3/8π.nu_ij^3
               DO J=0,NRAYS-1 ! Loop over rays
                  DO I=2,PARTICLE(P)%NPOINT(J) ! Loop over evaluation points
                     L = PARTICLE(P)%EVALUATION_POINT(J,I-1)
                     M = PARTICLE(P)%EVALUATION_POINT(J,I)
                     STEP_SIZE = SQRT((PARTICLE(M)%COORDINATES(1)-PARTICLE(L)%COORDINATES(1))**2 &
                                  & + (PARTICLE(M)%COORDINATES(2)-PARTICLE(L)%COORDINATES(2))**2 &
                                  & + (PARTICLE(M)%COORDINATES(3)-PARTICLE(L)%COORDINATES(3))**2)
#ifdef EULER
!                    Line opacity of the coolant transition (i,j) using the Euler integration method
                     FACTOR2 = 1.0D0/(PARTICLE(M)%COOLANT(N)%LINEWIDTH) ! = 1/δv_D
                     FACTOR3 = (PARTICLE(M)%COOLANT(N)%POPULATION(JLEVEL)*COOLANT(N)%WEIGHT(ILEVEL)  & ! = (n_j.g_i - n_i.g_j)/g_j = n_i.(n_j.g_i/n_i.g_j - 1)
                            & - PARTICLE(M)%COOLANT(N)%POPULATION(ILEVEL)*COOLANT(N)%WEIGHT(JLEVEL)) &
                            &  /COOLANT(N)%WEIGHT(JLEVEL)
                     PARTICLE(P)%COOLANT(N)%OPACITY(J,ILEVEL,JLEVEL) = PARTICLE(P)%COOLANT(N)%OPACITY(J,ILEVEL,JLEVEL) &
                                                                   & + FACTOR1*FACTOR2*FACTOR3*STEP_SIZE ! dtau_ij = A_ij.c^3/8π.nu_ij^3 * 1/δv_D * n_i.(n_j.g_i/n_i.g_j - 1) * dr
#else
!                    Line opacity of the coolant transition (i,j) using the Trapezium integration method (default)
                     FACTOR2 = 2.0D0/(PARTICLE(L)%COOLANT(N)%LINEWIDTH+PARTICLE(M)%COOLANT(N)%LINEWIDTH) ! = 1/δv_D
                     FACTOR3 = ((PARTICLE(L)%COOLANT(N)%POPULATION(JLEVEL)*COOLANT(N)%WEIGHT(ILEVEL) & ! = (n_j.g_i - n_i.g_j)/g_j = n_i.(n_j.g_i/n_i.g_j - 1)
                            &   +PARTICLE(M)%COOLANT(N)%POPULATION(JLEVEL)*COOLANT(N)%WEIGHT(ILEVEL))/2.0  &
                            & - (PARTICLE(L)%COOLANT(N)%POPULATION(ILEVEL)*COOLANT(N)%WEIGHT(JLEVEL) &
                            &   +PARTICLE(M)%COOLANT(N)%POPULATION(ILEVEL)*COOLANT(N)%WEIGHT(JLEVEL))/2.0) &
                            &  /COOLANT(N)%WEIGHT(JLEVEL)
                     PARTICLE(P)%COOLANT(N)%OPACITY(J,ILEVEL,JLEVEL) = PARTICLE(P)%COOLANT(N)%OPACITY(J,ILEVEL,JLEVEL) &
                                                                   & + FACTOR1*FACTOR2*FACTOR3*STEP_SIZE ! dtau_ij = A_ij.c^3/8π.nu_ij^3 * 1/δv_D * n_i.(n_j.g_i/n_i.g_j - 1) * dr
#endif
                  END DO ! End of loop over evaluation points
               END DO ! End of loop over rays
            END DO ! End of loop over levels (j)
         END DO ! End of loop over levels (i)
      END DO ! End of loop over coolants
   END DO ! End of loop over particles
!$OMP END PARALLEL DO

   RETURN
END SUBROUTINE CALCULATE_LINE_OPACITIES
!=======================================================================
