!=======================================================================
!
!  Calculate the level populations of a given coolant by constructing
!  the matrix of transition rates and solving to find the populations
!  assuming statistical equilibrium. The resulting set of statistical
!  equilibrium equations take the form:
!
!     n_i.∑_j R_ij = ∑_j n_j.R_ji
!  or:
!     n_i.∑_j R_ij - ∑_j n_j.R_ji = 0
!
!  where n_i is the population density (cm^-3) of level i and R_ij is
!  the transition rate (s^-1) from level i to level j. By rearranging
!  these equations, they can then be put into the matrix form: A.n=0,
!  where A is a coefficient matrix and n is a vector containing the N
!  population densities of all the levels. The right-hand side of the
!  matrix equation is then a vector of length N, composed of zeroes.
!
!  The elements of the coefficient matrix A are specified as follows:
!
!     A_ij = -R_ji    (j≠i)
!     A_ii = ∑_j R_ij (j≠i)
!
!  Since this set of equilibrium equations is not independent, one of
!  the equations has to be replaced by the conservation equation:
!
!     ∑_j n_j = n_tot
!
!  where n_tot is the density of the coolant species in all levels.
!
!  Therefore, the last row of the coefficient matrix is replaced with
!  this summation over all levels, and the last right-hand-side value
!  is set to the total density of the coolant species (cm^-3).
!
!  This system of linear equations is then solved using Gauss-Jordan
!  elimination to determine the population densities of all N levels.
!
!-----------------------------------------------------------------------
SUBROUTINE CALCULATE_LEVEL_POPULATIONS(NRAYS,ID,COOLANT,PARTICLE)

   USE HEALPIX_TYPES
   USE COOLANT_MODULE
   USE PARTICLE_MODULE
   USE SUBROUTINES_MODULE

   IMPLICIT NONE

   INTEGER(KIND=I4B),   INTENT(IN)    :: NRAYS
   INTEGER(KIND=I4B),   INTENT(IN)    :: ID
   TYPE(COOLANT_TYPE),  INTENT(IN)    :: COOLANT
   TYPE(PARTICLE_TYPE), INTENT(INOUT) :: PARTICLE

   INTEGER(KIND=I4B) :: I,J
   REAL(KIND=DP)     :: SUM
   REAL(KIND=DP), ALLOCATABLE :: R(:,:),A(:,:),B(:)

!  Allocate the arrays
   ALLOCATE(R(1:COOLANT%NLEVEL,1:COOLANT%NLEVEL))
   ALLOCATE(A(1:COOLANT%NLEVEL,1:COOLANT%NLEVEL))
   ALLOCATE(B(1:COOLANT%NLEVEL))

!  Construct the matrix of transition rates R_ij (s^-1)
   CALL CONSTRUCT_TRANSITION_MATRIX(NRAYS,ID,COOLANT,PARTICLE,R)

!  Fill the coefficient matrix A and the right-hand-side vector b
   DO I=1,COOLANT%NLEVEL
      SUM=0.0D0
      DO J=1,COOLANT%NLEVEL
         IF(J.EQ.I) CYCLE
         A(I,J) = -R(J,I)
         SUM=SUM+R(I,J)
      END DO
      A(I,I)=SUM
   END DO
   B=0.0D0

!  Replace the last equilibrium equation in the transition matrix with
!  the conservation equation (i.e. the sum of the population densities
!  over all levels), and replace the last entry in the right-hand-side
!  vector with the total density of the coolant species.
   A(COOLANT%NLEVEL,:)=1 ! Sum over all levels, ∑_j n_j
   B(COOLANT%NLEVEL)=PARTICLE%COOLANT(ID)%DENSITY

!  Call the Gauss-Jordan solver (the solution is returned in vector b)
   CALL GAUSS_JORDAN(COOLANT%NLEVEL,A,B)

!  Replace negative or NaN values caused by numerical noise around zero
   DO I=1,COOLANT%NLEVEL
      IF(.NOT.B(I).GE.0) B(I)=0.0D0
   END DO

!  Store the previously-calculated population densities for comparison
   PARTICLE%COOLANT(ID)%PREVIOUS_POPULATION=PARTICLE%COOLANT(ID)%POPULATION

!  Store the new population densities
   PARTICLE%COOLANT(ID)%POPULATION=B

   DEALLOCATE(R,A,B)
         
   RETURN
END SUBROUTINE CALCULATE_LEVEL_POPULATIONS
!=======================================================================

!=======================================================================
!
!  Standard linear equation solver using Gauss-Jordon elimination taken
!  directly from Numerical Recipes (Chapter B2).
!
!  A is an NxN input coefficient matrix, B is an input vector of size N
!  containing the right-hand-side values. On output, A is replaced by its
!  matrix inverse and B is replaced by the corresponding set of solution
!  values.
!
!-----------------------------------------------------------------------
SUBROUTINE GAUSS_JORDAN(N,A,B)

   USE HEALPIX_TYPES
   USE SWAP_FUNCTION

   IMPLICIT NONE

   INTEGER(KIND=I4B), INTENT(IN)    :: N
   REAL(KIND=DP),     INTENT(INOUT) :: A(:,:)
   REAL(KIND=DP),     INTENT(INOUT) :: B(:)

   INTEGER(KIND=I4B) :: I,J,K,L,IROW,ICOL
   INTEGER(KIND=I4B), ALLOCATABLE :: IPIV(:),INDEX_ROW(:),INDEX_COL(:)
   REAL(KIND=DP) :: MAX,DUMMY,PIVINV

   ALLOCATE(IPIV(1:N),INDEX_ROW(1:N),INDEX_COL(1:N))

   ICOL=0
   IROW=0
   IPIV=0

   DO I=1,N ! Main loop over columns to be reduced
      MAX=0.0D0
      DO J=1,N
         IF(IPIV(J).NE.1) THEN
            DO K=1,N
               IF(IPIV(K).EQ.0) THEN
                  IF(ISNAN(A(J,K))) THEN
                     WRITE(6,*)
                     WRITE(6,*) 'ERROR! NaN found in coefficient matrix A of Gauss-Jordan routine'
                     WRITE(6,"(' A(',I3,',',I3,') =',F4.1)") J,K,A(J,K)
                     WRITE(6,*) 'A ='
                     DO L=1,N
                        WRITE(6,"(100(PD9.1))") A(L,:)
                     END DO
                     WRITE(6,*)
                     STOP
                  END IF
                  IF(ABS(A(J,K)).GE.MAX) THEN
                     MAX=ABS(A(J,K))
                     IROW=J
                     ICOL=K
                  END IF
               ELSE IF(IPIV(K).GT.1) THEN
                  WRITE(6,*)
                  WRITE(6,*) 'ERROR! Singular matrix found in Gauss-Jordan routine (#1)'
                  WRITE(6,*)
                  STOP
               END IF
            END DO
         END IF
      END DO
      IPIV(ICOL)=IPIV(ICOL)+1
      IF(IROW.NE.ICOL) THEN
         CALL SWAP(A(IROW,:),A(ICOL,:))
         CALL SWAP(B(IROW),B(ICOL))
      END IF
      INDEX_ROW(I)=IROW
      INDEX_COL(I)=ICOL
      IF(A(ICOL,ICOL).EQ.0.0D0) THEN
         WRITE(6,*)
         WRITE(6,*) 'ERROR! Singular matrix found in Gauss-Jordan routine (#2)'
         WRITE(6,*)
         STOP
      END IF
      PIVINV=1.0D0/A(ICOL,ICOL)
      A(ICOL,ICOL)=1.0D0
      A(ICOL,:)=A(ICOL,:)*PIVINV
      B(ICOL)=B(ICOL)*PIVINV
      DO L=1,N
         IF(L.NE.ICOL) THEN
            DUMMY=A(L,ICOL)
            A(L,ICOL)=0.0D0
            A(L,:)=A(L,:)-A(ICOL,:)*DUMMY
            B(L)=B(L)-B(ICOL)*DUMMY
         END IF
      END DO
   END DO ! End of main loop over columns

!  Unscramble the solution by interchanging pairs of columns
!  in the reverse order to which the permutation was built up
   DO L=N,1,-1
      CALL SWAP(A(:,INDEX_ROW(L)),A(:,INDEX_COL(L)))
   END DO

   DEALLOCATE(IPIV,INDEX_ROW,INDEX_COL)

   RETURN
END SUBROUTINE GAUSS_JORDAN
!=======================================================================
