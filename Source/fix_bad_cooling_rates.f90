!=======================================================================
!
!  Check for particles with cooling rates that deviate significantly
!  from those of their neighbours (i.e. those having rates more than
!  100 times higher/lower than their neighbours), and assign to them
!  new values based on the nearest-neighbour average.
!
!-----------------------------------------------------------------------
SUBROUTINE ADJUST_BAD_COOLING_RATES(NPART,NCOOL,NRAYS,COOLANT,PARTICLE)

   USE HEALPIX_TYPES
   USE COOLANT_MODULE
   USE PARTICLE_MODULE

   IMPLICIT NONE

   INTEGER(KIND=I4B),   INTENT(IN)    :: NPART,NCOOL,NRAYS
   TYPE(COOLANT_TYPE),  INTENT(IN)    :: COOLANT(*)
   TYPE(PARTICLE_TYPE), INTENT(INOUT) :: PARTICLE(*)

   INTEGER(KIND=I4B) :: K,L,M,N,P,NEIGHBOURS,NCHANGED(1:NCOOL)
   REAL(KIND=DP) :: THRESHOLD,MIN_RATE,MAX_RATE,AVG_RATE

   ! Specify the threshold value
   THRESHOLD = 100.0D0

   ! Initialize the counter
   NCHANGED = 0

   DO N=1,NCOOL ! Loop over coolants

      ! Ignore the Lyman-alpha cooling
      IF(COOLANT(N)%NAME.EQ."H") CYCLE

      DO P=1,NPART ! Loop over particles

         ! Find the minimum, maximum and log-average cooling rate from the nearest neighbour particles
         NEIGHBOURS = 0 ; MIN_RATE = 1.0D99 ; MAX_RATE = 0.0D0 ; AVG_RATE = 0.0D0
         DO K=0,NRAYS-1 ! Loop over rays
            IF(PARTICLE(P)%NPOINT(K).GE.2) THEN ! One or more neighbouring particles exist along this ray
               L = PARTICLE(P)%EVALUATION_POINT(K,2)

               ! If more particles exist along the ray, check that the neighbouring particle (L) does not have a "bad" rate
               IF(PARTICLE(P)%NPOINT(K).GE.3) THEN
                  M = PARTICLE(P)%EVALUATION_POINT(K,3)

                  ! If the neighouring particle (L) has a cooling rate similar to that of its neighbour (M) then include it
                  IF(PARTICLE(L)%COOLING_RATE(N).LE.PARTICLE(M)%COOLING_RATE(N)*THRESHOLD .OR. &
                     PARTICLE(L)%COOLING_RATE(N).GE.PARTICLE(M)%COOLING_RATE(N)/THRESHOLD) THEN
                     MIN_RATE = MIN(MIN_RATE,PARTICLE(L)%COOLING_RATE(N))
                     MAX_RATE = MAX(MAX_RATE,PARTICLE(L)%COOLING_RATE(N))
                     AVG_RATE = AVG_RATE + LOG10(ABS(PARTICLE(L)%COOLING_RATE(N)))
                     NEIGHBOURS = NEIGHBOURS + 1
                  END IF

               ! If only one neighbouring particle exists along the ray, include it regardless of whether its rate is "good" or "bad"
               ELSE
                  MIN_RATE = MIN(MIN_RATE,PARTICLE(L)%COOLING_RATE(N))
                  MAX_RATE = MAX(MAX_RATE,PARTICLE(L)%COOLING_RATE(N))
                  AVG_RATE = AVG_RATE + LOG10(ABS(PARTICLE(L)%COOLING_RATE(N)))
                  NEIGHBOURS = NEIGHBOURS + 1
               END IF

            END IF
         END DO ! End of loop over rays

         IF(NEIGHBOURS.GT.0) THEN
            AVG_RATE = 10.0D0**(AVG_RATE/REAL(NEIGHBOURS,KIND=DP))
         ELSE ! In the (unlikely) event of no neighbouring particles, use the cooling rate of the particle itself
            MIN_RATE = PARTICLE(P)%COOLING_RATE(N)
            MAX_RATE = PARTICLE(P)%COOLING_RATE(N)
            AVG_RATE = PARTICLE(P)%COOLING_RATE(N)
         END IF

         ! Check if the cooling rate for the current particle is more than THRESHOLD times
         ! higher/lower than the maximum/minimum cooling rate of its neighbours. If so, set
         ! its cooling rate to the nearest-neighbour log-average and increment the counter.
         IF(PARTICLE(P)%COOLING_RATE(N).GT.MAX_RATE*THRESHOLD .OR. &
          & PARTICLE(P)%COOLING_RATE(N).LT.MIN_RATE/THRESHOLD) THEN
            PARTICLE(P)%COOLING_RATE(N) = AVG_RATE
            NCHANGED(N) = NCHANGED(N) + 1
         END IF

      END DO ! End of loop over particles
   END DO ! End of loop over coolants

   WRITE(6,"(' Particles affected:',10(' [',A,': ',I4,']',:))") (TRIM(ADJUSTL(COOLANT(N)%NAME)),NCHANGED(N),N=1,NCOOL)

   RETURN
END SUBROUTINE ADJUST_BAD_COOLING_RATES
!=======================================================================
