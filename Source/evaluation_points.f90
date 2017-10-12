!=======================================================================
!
!  Find the evaluation points along each HEALPix ray for each particle.
!
!-----------------------------------------------------------------------
SUBROUTINE FIND_EVALUATION_POINTS(NPART,NRAYS,PARTICLE)

#ifdef OPENMP
   USE OMP_LIB
#endif
   USE HEALPIX_TYPES
   USE PARTICLE_MODULE

   IMPLICIT NONE

   INTEGER(KIND=I4B),   INTENT(IN)    :: NPART,NRAYS
   TYPE(PARTICLE_TYPE), INTENT(INOUT) :: PARTICLE(*)

   INTEGER(KIND=I4B) :: I,J,L,P
   INTEGER(KIND=I4B) :: EVALUATION_POINT(0:NRAYS-1,1:NPART)
   REAL(KIND=DP) :: RAY_VECTOR(0:13,1:3),UNIT_VECTOR(1:3)
   REAL(KIND=DP) :: DISTANCE(0:NRAYS-1,1:NPART)
   REAL(KIND=DP), PARAMETER :: R = 1.0D0/SQRT(3.0D0)

!  Sorting variables (Shell's method; see Numerical Recipes, p.323)
   INTEGER(KIND=I4B) :: INC,INDEX
   REAL(KIND=DP) :: VALUE

!  Define the unit vector direction of each ray in Cartesian coordinates
   DATA RAY_VECTOR(0,:)/-1,0,0/
   DATA RAY_VECTOR(1,:)/+1,0,0/
   DATA RAY_VECTOR(2,:)/0,-1,0/
   DATA RAY_VECTOR(3,:)/0,+1,0/
   DATA RAY_VECTOR(4,:)/0,0,-1/
   DATA RAY_VECTOR(5,:)/0,0,+1/

   DATA RAY_VECTOR(6,:) /-R,-R,-R/
   DATA RAY_VECTOR(7,:) /-R,+R,-R/
   DATA RAY_VECTOR(8,:) /+R,+R,-R/
   DATA RAY_VECTOR(9,:) /+R,-R,-R/
   DATA RAY_VECTOR(10,:)/-R,-R,+R/
   DATA RAY_VECTOR(11,:)/-R,+R,+R/
   DATA RAY_VECTOR(12,:)/+R,+R,+R/
   DATA RAY_VECTOR(13,:)/+R,-R,+R/

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(DYNAMIC) &
!$OMP SHARED(NPART,NRAYS,PARTICLE,RAY_VECTOR) &
!$OMP PRIVATE(EVALUATION_POINT,UNIT_VECTOR,DISTANCE) &
!$OMP PRIVATE(I,J,L,INC,INDEX,VALUE)
   DO P=1,NPART

!     Allocate the array of total number of evaluation points along each ray
      ALLOCATE(PARTICLE(P)%NPOINT(0:NRAYS-1))

!     Find the evaluation points along each ray
      EVALUATION_POINT = 0
      EVALUATION_POINT(:,1) = P
      PARTICLE(P)%NPOINT = 1
      DISTANCE = 0.0D0
      DO J=0,NRAYS-1 ! Loop over rays
         I = 2
         DO L=1,NPART ! Search loop over all particles
            IF(L.EQ.P) CYCLE ! The particle of interest (P) has already been included in the list
            IF(PARTICLE(L)%TYPE.LT.2) CYCLE ! Do not include ionized (type 1) or undefined (type 0) particles

!           Calculate the distance of the current particle (L) from the particle of interest (P)
            DISTANCE(J,I) = SQRT((PARTICLE(L)%COORDINATES(1)-PARTICLE(P)%COORDINATES(1))**2 &
                             & + (PARTICLE(L)%COORDINATES(2)-PARTICLE(P)%COORDINATES(2))**2 &
                             & + (PARTICLE(L)%COORDINATES(3)-PARTICLE(P)%COORDINATES(3))**2)
            IF(DISTANCE(J,I).EQ.0) CYCLE

!           Calculate the unit vector direction of the particle (L) from the particle of interest (P)
            UNIT_VECTOR(1) = (PARTICLE(L)%COORDINATES(1)-PARTICLE(P)%COORDINATES(1))
            UNIT_VECTOR(2) = (PARTICLE(L)%COORDINATES(2)-PARTICLE(P)%COORDINATES(2))
            UNIT_VECTOR(3) = (PARTICLE(L)%COORDINATES(3)-PARTICLE(P)%COORDINATES(3))
            UNIT_VECTOR = UNIT_VECTOR/DISTANCE(J,I)

!           If the particle falls in the direction of the current ray, add it to the list of evaluation points
            IF(UNIT_VECTOR(1).EQ.RAY_VECTOR(J,1) .AND. &
               UNIT_VECTOR(2).EQ.RAY_VECTOR(J,2) .AND. &
               UNIT_VECTOR(3).EQ.RAY_VECTOR(J,3)) THEN
               EVALUATION_POINT(J,I) = L
               PARTICLE(P)%NPOINT(J) = PARTICLE(P)%NPOINT(J)+1
               I = I+1
            END IF
         END DO ! End of search loop over all particles
      END DO ! End of loop over rays

!     Sort the evaluation points by increasing distance from the particle of interest using Shell's method
!     (diminishing increment sort). The algorithm is taken directly from Numerical Recipes (p.323)
      DO J=0,NRAYS-1 ! Loop over rays
!        Determine the starting increment
         INC = 1
 1       INC = 3*INC+1
         IF(INC.LE.PARTICLE(P)%NPOINT(J)) GOTO 1
!        Loop over the partial sorts
 2       CONTINUE
            INC = INC/3
!           Outer loop of straight insertion
            DO I=INC+1,PARTICLE(P)%NPOINT(J)
               VALUE = DISTANCE(J,I) ; INDEX = EVALUATION_POINT(J,I)
               L = I
!              Inner loop of straight insertion
 3             IF(DISTANCE(J,L-INC).GT.VALUE) THEN
                  DISTANCE(J,L) = DISTANCE(J,L-INC) ; EVALUATION_POINT(J,L) = EVALUATION_POINT(J,L-INC)
                  L = L-INC
                  IF(L.LE.INC) GOTO 4
               GOTO 3
               END IF
 4             DISTANCE(J,L) = VALUE ; EVALUATION_POINT(J,L) = INDEX
            END DO
         IF(INC.GT.1) GOTO 2
      END DO ! End of loop over rays

!     Allocate memory for the required number of evaluation points and store their indices
      ALLOCATE(PARTICLE(P)%EVALUATION_POINT(0:NRAYS-1,1:MAXVAL(PARTICLE(P)%NPOINT)))
      PARTICLE(P)%EVALUATION_POINT = 0
      DO J=0,NRAYS-1
         PARTICLE(P)%EVALUATION_POINT(J,1:PARTICLE(P)%NPOINT(J)) = EVALUATION_POINT(J,1:PARTICLE(P)%NPOINT(J))
      END DO

   END DO ! End of loop over particles
!$OMP END PARALLEL DO

   RETURN
END SUBROUTINE FIND_EVALUATION_POINTS
!=======================================================================
