!=======================================================================
!
!  Calculate the column density (cm^-2) of each species along each
!  HEALPix ray to the PDR surface (or simulation boundary) for all
!  particles. Calculations can be performed using either the Euler
!  or Trapezium integration scheme (specified by a compiler flag).
!
!-----------------------------------------------------------------------
SUBROUTINE CALCULATE_COLUMN_DENSITIES(NPART,NRAYS,NSPEC,PARTICLE)

   USE HEALPIX_TYPES
   USE PARTICLE_MODULE

   IMPLICIT NONE

   INTEGER(KIND=I4B),   INTENT(IN)    :: NPART,NRAYS,NSPEC
   TYPE(PARTICLE_TYPE), INTENT(INOUT) :: PARTICLE(*)

   INTEGER(KIND=I4B) :: I,J,K,L,M,P

   REAL(KIND=DP) :: STEP_SIZE

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(DYNAMIC) &
!$OMP    SHARED(NPART,NRAYS,NSPEC,PARTICLE) &
!$OMP    PRIVATE(I,J,K,L,M,STEP_SIZE)
   DO P=1,NPART ! Loop over particles

!     Calculate the total hydrogen column density (N_H) and
!     column density of each species along each HEALPix ray
      PARTICLE(P)%TOTAL_COLUMN = 0.0D0
      PARTICLE(P)%COLUMN_DENSITY = 0.0D0
      DO J=0,NRAYS-1 ! Loop over rays
         DO I=2,PARTICLE(P)%NPOINT(J) ! Loop over evaluation points
            L = PARTICLE(P)%EVALUATION_POINT(J,I-1)
            M = PARTICLE(P)%EVALUATION_POINT(J,I)
            STEP_SIZE = SQRT((PARTICLE(M)%COORDINATES(1)-PARTICLE(L)%COORDINATES(1))**2 &
                         & + (PARTICLE(M)%COORDINATES(2)-PARTICLE(L)%COORDINATES(2))**2 &
                         & + (PARTICLE(M)%COORDINATES(3)-PARTICLE(L)%COORDINATES(3))**2)
#ifdef EULER
!           Total hydrogen column density using the Euler integration method
            PARTICLE(P)%TOTAL_COLUMN(J) = PARTICLE(P)%TOTAL_COLUMN(J) &
                                      & + PARTICLE(M)%GAS_DENSITY &
                                      &  *STEP_SIZE
#else
!           Total hydrogen column density using the Trapezium integration method (default)
            PARTICLE(P)%TOTAL_COLUMN(J) = PARTICLE(P)%TOTAL_COLUMN(J) &
                                       & + (PARTICLE(L)%GAS_DENSITY+PARTICLE(M)%GAS_DENSITY)/2.0 &
                                       &  *STEP_SIZE
#endif
            DO K=1,NSPEC ! Loop over species
#ifdef EULER
!              Column density of each species using the Euler integration method
               PARTICLE(P)%COLUMN_DENSITY(J,K) = PARTICLE(P)%COLUMN_DENSITY(J,K) &
                                             & + PARTICLE(M)%ABUNDANCE(K) &
                                             &  *PARTICLE(M)%GAS_DENSITY &
                                             &  *STEP_SIZE
#else
!              Column density of each species using the Trapezium integration method (default)
               PARTICLE(P)%COLUMN_DENSITY(J,K) = PARTICLE(P)%COLUMN_DENSITY(J,K) &
                                             & + (PARTICLE(L)%ABUNDANCE(K)*PARTICLE(L)%GAS_DENSITY &
                                             &   +PARTICLE(M)%ABUNDANCE(K)*PARTICLE(M)%GAS_DENSITY)/2.0 &
                                             &  *STEP_SIZE
#endif
            END DO ! End of loop over species
         END DO ! End of loop over evaluation points
      END DO ! End of loop over rays
   END DO ! End of loop over particles
!$OMP END PARALLEL DO

   RETURN
END SUBROUTINE CALCULATE_COLUMN_DENSITIES
!=======================================================================
