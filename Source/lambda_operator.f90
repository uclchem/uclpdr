!=======================================================================
!
!  Calculate the lambda operator for each allowed coolant transition,
!  needed to solve the level populations using the Accelerated Lambda
!  Iteration (ALI) method.
!
!-----------------------------------------------------------------------
SUBROUTINE CALCULATE_LAMBDA_OPERATOR(NPART,NCOOL,NRAYS,COOLANT,PARTICLE)

   USE HEALPIX_TYPES
   USE COOLANT_MODULE
   USE PARTICLE_MODULE

   IMPLICIT NONE

   INTEGER(KIND=I4B),   INTENT(IN)    :: NPART,NCOOL,NRAYS
   TYPE(COOLANT_TYPE),  INTENT(IN)    :: COOLANT(*)
   TYPE(PARTICLE_TYPE), INTENT(INOUT) :: PARTICLE(*)

   INTEGER(KIND=I4B) :: I,J,K,L,M,N,P
   REAL(KIND=DP) :: DTAU_1,DTAU_2,FACTOR_1,FACTOR_2
   REAL(KIND=DP) :: ALI_ij,ALI_ij_RAY(0:NRAYS-1)

   DO P=1,NPART ! Loop over particles
      DO N=1,NCOOL ! Loop over coolants
         IF(COOLANT(N)%CONVERGED) CYCLE
         PARTICLE(P)%COOLANT(N)%LAMBDA = 0.0D0
         DO I=1,COOLANT(N)%NLEVEL ! Loop over levels (i)
            DO J=1,COOLANT(N)%NLEVEL ! Loop over levels (j)
               IF(COOLANT(N)%A_COEFF(I,J).EQ.0) CYCLE
               DO K=0,NRAYS-1 ! Loop over rays
                  IF(PARTICLE(P)%NPOINT(K).GT.2) THEN ! More than one neighbouring particle exists along this ray

                     L = PARTICLE(P)%EVALUATION_POINT(K,1)
                     M = PARTICLE(P)%EVALUATION_POINT(K,2)

!                    Calculate the optical depth difference between the particle of interest and its neighbour
                     IF(PARTICLE(M)%AV(K).LT.PARTICLE(L)%AV(K)) THEN ! Ray direction is toward the PDR surface
                         DTAU_1 = PARTICLE(L)%COOLANT(N)%OPACITY(K,I,J)-PARTICLE(M)%COOLANT(N)%OPACITY(K,I,J)
                     ELSE ! Ray direction is away from the PDR surface
                         DTAU_1 = PARTICLE(M)%COOLANT(N)%OPACITY(K,I,J)-PARTICLE(L)%COOLANT(N)%OPACITY(K,I,J)
                     END IF

                     L = PARTICLE(P)%EVALUATION_POINT(K,2)
                     M = PARTICLE(P)%EVALUATION_POINT(K,3)

!                    Calculate the optical depth difference between the neighbouring particle and the next one
                     IF(PARTICLE(M)%AV(K).LT.PARTICLE(L)%AV(K)) THEN ! Ray direction is toward the PDR surface
                         DTAU_2 = PARTICLE(L)%COOLANT(N)%OPACITY(K,I,J)-PARTICLE(M)%COOLANT(N)%OPACITY(K,I,J)
                     ELSE ! Ray direction is away from the PDR surface
                         DTAU_2 = PARTICLE(M)%COOLANT(N)%OPACITY(K,I,J)-PARTICLE(L)%COOLANT(N)%OPACITY(K,I,J)
                     END IF

                  ELSE IF(PARTICLE(P)%NPOINT(K).EQ.2) THEN ! One neighbouring particle exists along this ray

                     L = PARTICLE(P)%EVALUATION_POINT(K,1)
                     M = PARTICLE(P)%EVALUATION_POINT(K,2)
                     DTAU_1 = PARTICLE(L)%COOLANT(N)%OPACITY(K,I,J)-PARTICLE(M)%COOLANT(N)%OPACITY(K,I,J)
                     DTAU_2 = 0.0D0

                  ELSE ! There are no neighbouring particles along this ray (the particle is at the PDR surface)

                     DTAU_1 = 0.0D0
                     DTAU_2 = 0.0D0

                  END IF

!                 For large opacities, invoke the diffusion approximation
                  IF(DTAU_1.GT.100.0D0 .AND. DTAU_2.GT.100.0D0) THEN
                      ALI_ij_RAY(K) = 0.0D0
                  ELSE IF((DTAU_1+DTAU_2).NE.0) THEN
                      IF(DTAU_1.NE.0) THEN
                         FACTOR_1 = 2.0D0*(1.0D0-EXP(-DTAU_1))/DTAU_1
                      ELSE
                         FACTOR_1 = 2.0D0
                      END IF
                      IF(DTAU_2.NE.0) THEN
                         FACTOR_2 = 2.0D0*(1.0D0-EXP(-DTAU_2))/DTAU_2
                      ELSE
                         FACTOR_2 = 2.0D0
                      END IF
                      ALI_ij_RAY(K) = 1.0D0/(1.0D0+(FACTOR_1+FACTOR_2)/(DTAU_1+DTAU_2)) - 1.0D0
                  ELSE
                      ALI_ij_RAY(K) = -1.0D0
                  END IF

!!$                  IF(ISNAN(ALI_ij_RAY(K))) THEN
!!$                     WRITE(6,*) 'ERROR! NaN value calculated for ALI:'
!!$                     WRITE(6,*) 'P =',P,'N =',N,'I =',I,'J =',J,'K =',K
!!$                     WRITE(6,*) 'dt_1 =',DTAU_1,'dt_2 =',DTAU_2
!!$                     WRITE(6,*) 'A_1  =',FACTOR_1,'A_2  =',FACTOR_2
!!$                     WRITE(6,*) 'ALI  =',ALI_ij_RAY(K)
!!$                     STOP
!!$                  END IF

               END DO ! End of loop over rays

!              The total lambda operator value is the average of the values along each ray
               ALI_ij=SUM(ALI_ij_RAY)/REAL(NRAYS,KIND=DP)
               PARTICLE(P)%COOLANT(N)%LAMBDA(I,J) = ALI_ij + 1.0D0

            END DO ! End of loop over levels (j)
         END DO ! End of loop over levels (i)
      END DO ! End of loop over coolants
   END DO ! End of loop over particles

   RETURN
END SUBROUTINE CALCULATE_LAMBDA_OPERATOR
!=======================================================================
