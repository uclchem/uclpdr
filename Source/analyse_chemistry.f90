!=======================================================================
!
!  Analyse the formation/destruction routes of the specified species
!  that are dominant in the given particle by examining the reaction
!  network and all rates involved. Save a breakdown of these results
!  to the specified output file.
!
!-----------------------------------------------------------------------
SUBROUTINE ANALYSE_CHEMISTRY(PREFIX,SUFFIX,NPART,NSPEC,NREAC,TIME, &
                           & SPECIES,REACTANT,PRODUCT,PARTICLE, &
                           & AV_LIST,SPECIES_LIST)

   USE HEALPIX_TYPES
   USE PARTICLE_MODULE
   USE NUM2STR_FUNCTION
   USE GLOBAL_MODULE, ONLY : nH,nH2,nHe,nelect

   IMPLICIT NONE

   CHARACTER(LEN=*),    INTENT(IN) :: PREFIX,SUFFIX
   INTEGER(KIND=I4B),   INTENT(IN) :: NPART,NSPEC,NREAC
   REAL(KIND=DP),       INTENT(IN) :: TIME,AV_LIST(:)
   CHARACTER(LEN=*),    INTENT(IN) :: SPECIES(:),SPECIES_LIST(:)
   CHARACTER(LEN=*),    INTENT(IN) :: REACTANT(:,:),PRODUCT(:,:)
   TYPE(PARTICLE_TYPE), INTENT(IN) :: PARTICLE(:)

   INTEGER(KIND=I4B)  :: I,J,K,L,M,N,P,INDEX,IP,ID
   INTEGER(KIND=I4B)  :: NPR(1:NREAC),NDR(1:NREAC)
   INTEGER(KIND=I4B)  :: LENGTH(4),NMAX(1),MULTIPLE,TOTAL
   REAL(KIND=DP)      :: ABUNDANCE1,ABUNDANCE2,PERCENTAGE
   REAL(KIND=DP)      :: PRODUCTION_RATE(1:NREAC),TOTAL_PRODUCTION
   REAL(KIND=DP)      :: DESTRUCTION_RATE(1:NREAC),TOTAL_DESTRUCTION
   CHARACTER(LEN=20)  :: LHS,RHS
   CHARACTER(LEN=256) :: FILENAME

   INTEGER(KIND=I4B), SAVE :: UNIT=1

!  Specify the output file name
   FILENAME=TRIM(ADJUSTL(PREFIX))//TRIM(ADJUSTL('.chem'))//TRIM(ADJUSTL(SUFFIX))
   WRITE(6,*) 'Writing file [',TRIM(FILENAME),']'
   OPEN(UNIT=UNIT,FILE=FILENAME,STATUS='REPLACE')

   DO INDEX=LBOUND(AV_LIST,1),UBOUND(AV_LIST,1) ! Loop over list of visual extinctions to analyse

!  Find index P of the particle with AV closest to
!  the current entry in the visual extinction list
   P=1
   DO I=1,NPART
      IF(ABS(PARTICLE(I)%AV(0)-AV_LIST(INDEX)).LT.ABS(PARTICLE(P)%AV(0)-AV_LIST(INDEX))) P=I
   END DO

   DO L=LBOUND(SPECIES_LIST,1),UBOUND(SPECIES_LIST,1) ! Loop over list of species to analyse

!     Find index I of the species that corresponds
!     to the current entry L in the species list
      DO I=1,NSPEC
         IF(SPECIES(I).EQ.SPECIES_LIST(L)) EXIT
      END DO

!     Check that species I corresponds to the desired species
!     Continue to the next entry in the species list if not
      IF(SPECIES(I).NE.SPECIES_LIST(L)) CYCLE

      WRITE(UNIT,1) TRIM(NUM2STR(P)),TIME,PARTICLE(P)%AV(0),PARTICLE(P)%GAS_DENSITY, &
                  & TRIM(NUM2STR(PARTICLE(P)%GAS_TEMPERATURE,1)),TRIM(SPECIES(I))

!     Reset the formation/destruction rate counters
      IP=0
      ID=0

!     Reset the total formation/destruction rates
      TOTAL_PRODUCTION=0.0D0
      TOTAL_DESTRUCTION=0.0D0

!     Check all reactions to find relevant ones
      DO J=1,NREAC

!        Reset the reactant abundances
         ABUNDANCE1=0.0D0
         ABUNDANCE2=0.0D0

!-----------------------------------------------------------------------
!        Formation routes for species I
!-----------------------------------------------------------------------

         IF(ANY(PRODUCT(J,:).EQ.SPECIES(I))) THEN

            IP=IP+1

!           Store the reactant abundances
            DO N=1,NSPEC
               IF(SPECIES(N).EQ.REACTANT(J,1)) ABUNDANCE1=PARTICLE(P)%ABUNDANCE(N)
               IF(SPECIES(N).EQ.REACTANT(J,2)) ABUNDANCE2=PARTICLE(P)%ABUNDANCE(N)
            END DO

            MULTIPLE=0
            IF(PRODUCT(J,1).EQ.SPECIES(I)) MULTIPLE=MULTIPLE+1
            IF(PRODUCT(J,2).EQ.SPECIES(I)) MULTIPLE=MULTIPLE+1
            IF(PRODUCT(J,3).EQ.SPECIES(I)) MULTIPLE=MULTIPLE+1
            IF(PRODUCT(J,4).EQ.SPECIES(I)) MULTIPLE=MULTIPLE+1

!           Calculate the reaction rate
            IF(REACTANT(J,2).EQ."PHOTON" .OR. &
             & REACTANT(J,2).EQ."CRP   " .OR. &
             & REACTANT(J,2).EQ."CRPHOT" .OR. &
             & REACTANT(J,2).EQ."XRAY  " .OR. &
             & REACTANT(J,2).EQ."FREEZE" .OR. &
             & REACTANT(J,2).EQ."CRH   " .OR. &
             & REACTANT(J,2).EQ."PHOTD " .OR. &
             & REACTANT(J,2).EQ."THERM ") THEN
               PRODUCTION_RATE(IP)=PARTICLE(P)%RATE(J)*ABUNDANCE1*MULTIPLE
            ELSE IF(REACTANT(J,2).EQ."XRSEC ") THEN
               IF(REACTANT(J,1).EQ."H2 ") THEN
                  PRODUCTION_RATE(IP)=PARTICLE(P)%RATE(J)*ABUNDANCE1/(41.9D0*(1.0D0+6.72D0*(1.83D0*PARTICLE(P)%ABUNDANCE(nelect)/(1.0D0+0.83D0*PARTICLE(P)%ABUNDANCE(nelect)))**0.824D0)*(PARTICLE(P)%ABUNDANCE(nH2)+0.53D0*PARTICLE(P)%ABUNDANCE(nH)))
               ELSE IF(REACTANT(J,1).EQ."He ") THEN
                  PRODUCTION_RATE(IP)=PARTICLE(P)%RATE(J)*ABUNDANCE1/(487.D0*(1.0D0+12.5D0*PARTICLE(P)%ABUNDANCE(nelect)**0.994D0)*(PARTICLE(P)%ABUNDANCE(nHe)))
               ELSE
                  PRODUCTION_RATE(IP)=PARTICLE(P)%RATE(J)*ABUNDANCE1/(39.8D0*(1.0D0+12.2D0*PARTICLE(P)%ABUNDANCE(nelect)**0.866D0)*(PARTICLE(P)%ABUNDANCE(nH)+1.89D0*PARTICLE(P)%ABUNDANCE(nH2)))
               END IF
            ELSE IF(REACTANT(J,2).EQ."XRLYA ") THEN
               PRODUCTION_RATE(IP)=PARTICLE(P)%RATE(J)*ABUNDANCE1/(39.8D0*(1.0D0+12.2D0*PARTICLE(P)%ABUNDANCE(nelect)**0.866D0)*(PARTICLE(P)%ABUNDANCE(nH)+1.89D0*PARTICLE(P)%ABUNDANCE(nH2)))
            ELSE IF(REACTANT(J,2).EQ."XRPHOT") THEN
               PRODUCTION_RATE(IP)=PARTICLE(P)%RATE(J)*ABUNDANCE1/(41.9D0*(1.0D0+6.72D0*(1.83D0*PARTICLE(P)%ABUNDANCE(nelect)/(1.0D0+0.83D0*PARTICLE(P)%ABUNDANCE(nelect)))**0.824D0)*(PARTICLE(P)%ABUNDANCE(nH2)+0.53D0*PARTICLE(P)%ABUNDANCE(nH)))
            ELSE
               PRODUCTION_RATE(IP)=PARTICLE(P)%RATE(J)*ABUNDANCE1*ABUNDANCE2*PARTICLE(P)%GAS_DENSITY*MULTIPLE
            END IF

            TOTAL_PRODUCTION=TOTAL_PRODUCTION+PRODUCTION_RATE(IP)
            NPR(IP)=J
         END IF

!-----------------------------------------------------------------------
!        Destruction routes for species I
!-----------------------------------------------------------------------

         IF(ANY(REACTANT(J,:).EQ.SPECIES(I))) THEN

            ID=ID+1

!           Store the reactant abundances
            DO N=1,NSPEC
               IF(SPECIES(N).EQ.REACTANT(J,1)) ABUNDANCE1=PARTICLE(P)%ABUNDANCE(N)
               IF(SPECIES(N).EQ.REACTANT(J,2)) ABUNDANCE2=PARTICLE(P)%ABUNDANCE(N)
            END DO

            MULTIPLE=0
            IF(REACTANT(J,1).EQ.SPECIES(I)) MULTIPLE=MULTIPLE+1
            IF(REACTANT(J,2).EQ.SPECIES(I)) MULTIPLE=MULTIPLE+1
            IF(REACTANT(J,3).EQ.SPECIES(I)) MULTIPLE=MULTIPLE+1

!           Calculate the reaction rate
            IF(REACTANT(J,2).EQ."PHOTON" .OR. &
             & REACTANT(J,2).EQ."CRP   " .OR. &
             & REACTANT(J,2).EQ."CRPHOT" .OR. &
             & REACTANT(J,2).EQ."XRAY  " .OR. &
             & REACTANT(J,2).EQ."FREEZE" .OR. &
             & REACTANT(J,2).EQ."CRH   " .OR. &
             & REACTANT(J,2).EQ."PHOTD " .OR. &
             & REACTANT(J,2).EQ."THERM ") THEN
               DESTRUCTION_RATE(ID)=PARTICLE(P)%RATE(J)*ABUNDANCE1*MULTIPLE
            ELSE IF(REACTANT(J,2).EQ."XRSEC ") THEN
               IF(REACTANT(J,1).EQ."H2 ") THEN
                  DESTRUCTION_RATE(ID)=PARTICLE(P)%RATE(J)*ABUNDANCE1/(41.9D0*(1.0D0+6.72D0*(1.83D0*PARTICLE(P)%ABUNDANCE(nelect)/(1.0D0+0.83D0*PARTICLE(P)%ABUNDANCE(nelect)))**0.824D0)*(PARTICLE(P)%ABUNDANCE(nH2)+0.53D0*PARTICLE(P)%ABUNDANCE(nH)))
               ELSE IF(REACTANT(J,1).EQ."He ") THEN
                  DESTRUCTION_RATE(ID)=PARTICLE(P)%RATE(J)*ABUNDANCE1/(487.D0*(1.0D0+12.5D0*PARTICLE(P)%ABUNDANCE(nelect)**0.994D0)*(PARTICLE(P)%ABUNDANCE(nHe)))
               ELSE
                  DESTRUCTION_RATE(ID)=PARTICLE(P)%RATE(J)*ABUNDANCE1/(39.8D0*(1.0D0+12.2D0*PARTICLE(P)%ABUNDANCE(nelect)**0.866D0)*(PARTICLE(P)%ABUNDANCE(nH)+1.89D0*PARTICLE(P)%ABUNDANCE(nH2)))
               END IF
            ELSE IF(REACTANT(J,2).EQ."XRLYA ") THEN
               DESTRUCTION_RATE(ID)=PARTICLE(P)%RATE(J)*ABUNDANCE1/(39.8D0*(1.0D0+12.2D0*PARTICLE(P)%ABUNDANCE(nelect)**0.866D0)*(PARTICLE(P)%ABUNDANCE(nH)+1.89D0*PARTICLE(P)%ABUNDANCE(nH2)))
            ELSE IF(REACTANT(J,2).EQ."XRPHOT") THEN
               DESTRUCTION_RATE(ID)=PARTICLE(P)%RATE(J)*ABUNDANCE1/(41.9D0*(1.0D0+6.72D0*(1.83D0*PARTICLE(P)%ABUNDANCE(nelect)/(1.0D0+0.83D0*PARTICLE(P)%ABUNDANCE(nelect)))**0.824D0)*(PARTICLE(P)%ABUNDANCE(nH2)+0.53D0*PARTICLE(P)%ABUNDANCE(nH)))
            ELSE
               DESTRUCTION_RATE(ID)=PARTICLE(P)%RATE(J)*ABUNDANCE1*ABUNDANCE2*PARTICLE(P)%GAS_DENSITY*MULTIPLE
            END IF

            TOTAL_DESTRUCTION=TOTAL_DESTRUCTION+DESTRUCTION_RATE(ID)
            NDR(ID)=J
         END IF

      END DO ! End of loop over reactions

!     Prevent divide-by-zero errors by setting minimum finite
!     values for the total formation and destruction rates
      IF(TOTAL_PRODUCTION.LT.1.0D-99) TOTAL_PRODUCTION=1.0D-99
      IF(TOTAL_DESTRUCTION.LT.1.0D-99) TOTAL_DESTRUCTION=1.0D-99

!-----------------------------------------------------------------------
!     Find the lengths of the longest reactant and product names
!-----------------------------------------------------------------------

      LENGTH=0
      DO N=1,IP
!        Only consider reactions that represent > 0.5% of the total rate
         IF(1.0D2*(PRODUCTION_RATE(N)/TOTAL_PRODUCTION).GT.0.5D0) THEN
            LENGTH(1)=MAX(LENGTH(1),LEN_TRIM(REACTANT(NPR(N),1)),LEN_TRIM(SPECIES(I)))
            LENGTH(2)=MAX(LENGTH(2),LEN_TRIM(TRIM(ADJUSTL(REACTANT(NPR(N),2)))//'   '// &
                                           & TRIM(ADJUSTL(REACTANT(NPR(N),3)))))
         END IF
      END DO
      DO N=1,ID
!        Only consider reactions that represent > 0.5% of the total rate
         IF(1.0D2*(DESTRUCTION_RATE(N)/TOTAL_DESTRUCTION).GT.0.5D0) THEN
            LENGTH(2)=MAX(LENGTH(2),LEN_TRIM(TRIM(ADJUSTL(REACTANT(NDR(N),2)))//'   '// &
                                           & TRIM(ADJUSTL(REACTANT(NDR(N),3)))))
            LENGTH(3)=MAX(LENGTH(3),LEN_TRIM(PRODUCT(NDR(N),1)),LEN_TRIM(SPECIES(I)))
         END IF
      END DO

!-----------------------------------------------------------------------
!     Output the formation reactions and their rates
!-----------------------------------------------------------------------

      NMAX=0
      TOTAL=0

!     List the formation reactions in order of decreasing importance
      DO M=1,IP

!        Find the location of the maximum value
         NMAX=MAXLOC(PRODUCTION_RATE(1:IP))
         N=NMAX(1)

!        Exit the loop once the reaction rates reach zero
         IF(PRODUCTION_RATE(N).EQ.0.0D0) EXIT

!        Calculate the percentage of the total formation
!        rate that is contributed by the current reaction
         PERCENTAGE=1.0D2*(PRODUCTION_RATE(N)/TOTAL_PRODUCTION)
         IF(PERCENTAGE.LT.0.5D0) PERCENTAGE=1.0D0
         TOTAL=TOTAL+NINT(PERCENTAGE)

!        Find the first instance of the species of interest in the list of products
         NMAX=MAXLOC(PRODUCT(NPR(N),:),(PRODUCT(NPR(N),:).EQ.SPECIES(I)))
         K=NMAX(1)

!        Create the left- and right-hand sides of the reaction string
         LHS=REACTANT(NPR(N),1)
         DO J=2,3 ! Loop over reactants
            IF(REACTANT(NPR(N),J).NE." ") THEN
               IF(LEN_TRIM(LHS).EQ.LEN_TRIM(REACTANT(NPR(N),1))) THEN
                  LHS=LHS(1:LENGTH(1))//' + '//REACTANT(NPR(N),J)
               ELSE
                  LHS=TRIM(ADJUSTL(LHS))//' + '//REACTANT(NPR(N),J)
               END IF
            END IF
         END DO ! End of loop over reactants

         RHS=PRODUCT(NPR(N),K)
         DO J=1,4 ! Loop over products
            IF(J.EQ.K) CYCLE ! Skip this one since it has already been added to the string
            IF(PRODUCT(NPR(N),J).NE." ") THEN
               IF(LEN_TRIM(RHS).EQ.LEN_TRIM(PRODUCT(NPR(N),K))) THEN
                  RHS=RHS(1:LENGTH(3))//' + '//PRODUCT(NPR(N),J)
               ELSE
                  RHS=TRIM(ADJUSTL(RHS))//' + '//PRODUCT(NPR(N),J)
               END IF
            END IF
         END DO ! End of loop over products

!        Right-align the left-hand side of the reaction string
         LHS=REPEAT(' ',LEN(LHS)-LENGTH(1)-LENGTH(2)-3)//ADJUSTL(LHS)

!        Left-align the right-hand side of the reaction string
         RHS=ADJUSTL(RHS)

!        Write the reaction and its percentage of the total formation rate
         WRITE(UNIT,2) LHS,RHS,NINT(PERCENTAGE)

!        Exit the loop once the sum of the reaction rates reaches 100%
         IF(TOTAL.GE.100) EXIT

!        Set the formation rate of this reaction to zero
!        to prevent it from being included more than once
         PRODUCTION_RATE(N)=0.0D0

      END DO ! End of loop over formation reactions
      WRITE(UNIT,*)

!-----------------------------------------------------------------------
!     Output the destruction reactions and their rates
!-----------------------------------------------------------------------

      NMAX=0
      TOTAL=0

!     List the destruction reactions in order of decreasing importance
      DO M=1,ID

!        Find the location of the maximum value
         NMAX=MAXLOC(DESTRUCTION_RATE(1:ID))
         N=NMAX(1)

!        Exit the loop once the reaction rates reach zero
         IF(DESTRUCTION_RATE(N).EQ.0.0D0) EXIT

!        Calculate the percentage of the total destruction
!        rate that is contributed by the current reaction
         PERCENTAGE=1.0D2*(DESTRUCTION_RATE(N)/TOTAL_DESTRUCTION)
         IF(PERCENTAGE.LT.0.5D0) PERCENTAGE=1.0D0
         TOTAL=TOTAL+NINT(PERCENTAGE)

!        Find the first instance of the species of interest in the list of reactants
         NMAX=MAXLOC(REACTANT(NDR(N),:),(REACTANT(NDR(N),:).EQ.SPECIES(I)))
         K=NMAX(1)

!        Create the left- and right-hand sides of the reaction string
         LHS=REACTANT(NDR(N),K)
         DO J=1,3 ! Loop over reactants
            IF(J.EQ.K) CYCLE ! Skip this one since it has already been added to the string
            IF(REACTANT(NDR(N),J).NE." ") THEN
               IF(LEN_TRIM(LHS).EQ.LEN_TRIM(REACTANT(NDR(N),K))) THEN
                  LHS=LHS(1:LENGTH(1))//' + '//REACTANT(NDR(N),J)
               ELSE
                  LHS=TRIM(ADJUSTL(LHS))//' + '//REACTANT(NDR(N),J)
               END IF
            END IF
         END DO ! End of loop over reactants

         RHS=PRODUCT(NDR(N),1)
         DO J=2,4 ! Loop over products
            IF(PRODUCT(NDR(N),J).NE." ") THEN
               IF(LEN_TRIM(RHS).EQ.LEN_TRIM(PRODUCT(NDR(N),1))) THEN
                  RHS=RHS(1:LENGTH(3))//' + '//PRODUCT(NDR(N),J)
               ELSE
                  RHS=TRIM(ADJUSTL(RHS))//' + '//PRODUCT(NDR(N),J)
               END IF
            END IF
         END DO ! End of loop over products

!        Right-align the left-hand side of the reaction string
         LHS=REPEAT(' ',LEN(LHS)-LENGTH(1)-LENGTH(2)-3)//ADJUSTL(LHS)

!        Left-align the right-hand side of the reaction string
         RHS=ADJUSTL(RHS)

!        Write the reaction and its percentage of the total destruction rate
         WRITE(UNIT,2) LHS,RHS,-NINT(PERCENTAGE)

!        Exit the loop once the sum of the reaction rates reaches 100%
         IF(TOTAL.GE.100) EXIT

!        Set the destruction rate of this reaction to zero
!        to prevent it from being included more than once
         DESTRUCTION_RATE(N)=0.0D0

      END DO ! End of loop over destruction reactions
      WRITE(UNIT,*)

!-----------------------------------------------------------------------
!     Output the abundance and total formation/destruction rates
!-----------------------------------------------------------------------

      WRITE(UNIT,3) PARTICLE(P)%ABUNDANCE(I),TOTAL_PRODUCTION,TOTAL_DESTRUCTION
      IF(L.LT.UBOUND(SPECIES_LIST,1)) WRITE(UNIT,4)

   END DO ! End of loop over species list
   WRITE(UNIT,5)

   END DO ! End of loop over visual extinction list
   CLOSE(UNIT)

 1 FORMAT('Particle = ',A,', t =',ES7.1E1,' yr, A_V =',ES7.1E1,' mag',/, &
        & 'n_H =',ES7.1E1,' cm^-3, T_gas = ',A,' K',/, &
        & 'Species = ',A,/)
 2 FORMAT(A20,' --> ',A20,'Rate:',SP,I4,'%')
 3 FORMAT('Total Abundance  =',ES10.3,/, &
        & 'Formation Rate   =',ES10.3,' s^-1',/, &
        & 'Destruction Rate =',ES10.3,' s^-1',/)
 4 FORMAT(80('-'),/)
 5 FORMAT(80('='),/)

   RETURN
END SUBROUTINE ANALYSE_CHEMISTRY
!=======================================================================
