!=======================================================================
!
!  Read in the chemical reactions, the reactants and products involved,
!  the Arrhenius equation parameters that describe them and the minimum
!  and maximum temperatures within which they are valid. The specified
!  file is assumed to contain comma separated values (CSV) format data.
!  This is in line with the new Rate06 standard, removing the need for
!  file-dependent FORMAT statements.
!
!-----------------------------------------------------------------------
SUBROUTINE READ_REACTIONS(FILENAME,NREAC,REACTANT,PRODUCT,ALPHA,BETA,GAMMA,RTMIN,RTMAX,DUPLICATE)

   USE HEALPIX_TYPES
   USE NUM2STR_FUNCTION
   USE GLOBAL_MODULE, ONLY : nRGR,nRH2,nRHD,nRCO,nRCI,nRSI,nR_H2x_1,nR_H2x_2, &
                           & nR_H3x_1,nR_H3x_2,nR_H3Ox_1,nR_H3Ox_2,nR_H3Ox_3, &
                           & nR_Hex_1,nR_Hex_2,nR_Hex_3,nR_Hex_4,nR_HCOx_1

   IMPLICIT NONE

   CHARACTER(LEN=*),  INTENT(IN)  :: FILENAME
   INTEGER(KIND=I4B), INTENT(IN)  :: NREAC
   CHARACTER(LEN=10), INTENT(OUT) :: REACTANT(1:NREAC,1:3),PRODUCT(1:NREAC,1:4)
   REAL(KIND=DP),     INTENT(OUT) :: ALPHA(1:NREAC),BETA(1:NREAC),GAMMA(1:NREAC)
   REAL(KIND=DP),     INTENT(OUT) :: RTMIN(1:NREAC),RTMAX(1:NREAC)
   INTEGER(KIND=I4B), INTENT(OUT) :: DUPLICATE(1:NREAC)

   INTEGER(KIND=I4B) :: I,J,INDEX,IER
   CHARACTER(LEN=1)  :: CLEM

!  Open the input file
   OPEN(UNIT=1,FILE=FILENAME,IOSTAT=IER,ACTION='READ',STATUS='OLD')

!  Produce an error message if the file does not exist (or cannot be opened for whatever reason)
   IF(IER.NE.0) THEN
      WRITE(6,*) 'ERROR! Cannot open reaction rate file ',TRIM(FILENAME),' for input'
      WRITE(6,*)
      CLOSE(1)
      STOP
   END IF

!  Initialize the arrays
   REACTANT=""
   PRODUCT=""
   ALPHA=0.0D0
   BETA=0.0D0
   GAMMA=0.0D0
   RTMIN=0.0D0
   RTMAX=0.0D0
   DUPLICATE=0

!  Initialize the reaction index numbers. If they are not assigned
!  subsequently an attempt to access those reactions will generate
!  an error and the code will crash. This is a useful bug catch.
   nRGR=0 ; nRH2=0 ; nRHD=0 ; nRCO=0 ; nRCI=0 ; nRSI=0
   nR_H2x_1=0 ; nR_H2x_2=0 ; nR_H3x_1=0 ; nR_H3x_2=0
   nR_Hex_1=0 ; nR_Hex_2=0 ; nR_Hex_3=0 ; nR_Hex_4=0
   nR_H3Ox_1=0; nR_H3Ox_2=0; nR_H3Ox_3=0; nR_HCOx_1=0

   DO I=1,NREAC

!     Read the reactant and product names, Arrhenius equation parameters and minimum/maximum temperatures within which they are valid
      READ(1,*,IOSTAT=IER) INDEX,(REACTANT(I,J),J=1,3),(PRODUCT(I,J),J=1,4), &
                           & ALPHA(I),BETA(I),GAMMA(I),CLEM,RTMIN(I),RTMAX(I)

!     Produce an error message if an unexpected data format is encountered
      IF(IER.NE.0) THEN
         WRITE(6,*) 'ERROR! Unexpected data format on line ',TRIM(NUM2STR(I)),' of reaction rate file'
         WRITE(6,*)
         CLOSE(1)
         STOP
      END IF

!     Assign the various index numbers to their correct reactions
      IF((REACTANT(I,1).EQ."H  " .AND. REACTANT(I,2).EQ."H     " .AND. REACTANT(I,3).EQ." ") .AND. &
       & (PRODUCT(I,1).EQ."H2  " .AND. PRODUCT(I,2).EQ."   "))  nRGR=I
      IF((REACTANT(I,1).EQ."H  " .AND. REACTANT(I,2).EQ."H     " .AND. REACTANT(I,3).EQ."#") .AND. &
       & (PRODUCT(I,1).EQ."H2  " .AND. PRODUCT(I,2).EQ."#  "))  nRGR=I
      IF((REACTANT(I,1).EQ."H2 " .AND. REACTANT(I,2).EQ."PHOTON" .AND. REACTANT(I,3).EQ." ") .AND. &
       & (PRODUCT(I,1).EQ."H   " .AND. PRODUCT(I,2).EQ."H  "))  nRH2=I
      IF((REACTANT(I,1).EQ."HD " .AND. REACTANT(I,2).EQ."PHOTON" .AND. REACTANT(I,3).EQ." ") .AND. &
      & ((PRODUCT(I,1).EQ."H   " .AND. PRODUCT(I,2).EQ."D  ") .OR. &
      &  (PRODUCT(I,1).EQ."D   " .AND. PRODUCT(I,2).EQ."H  "))) nRHD=I
      IF((REACTANT(I,1).EQ."CO " .AND. REACTANT(I,2).EQ."PHOTON" .AND. REACTANT(I,3).EQ." ") .AND. &
      & ((PRODUCT(I,1).EQ."C   " .AND. PRODUCT(I,2).EQ."O  ") .OR. &
      &  (PRODUCT(I,1).EQ."O   " .AND. PRODUCT(I,2).EQ."C  "))) nRCO=I
      IF((REACTANT(I,1).EQ."C  " .AND. REACTANT(I,2).EQ."PHOTON" .AND. REACTANT(I,3).EQ." ") .AND. &
      & ((PRODUCT(I,1).EQ."C+  " .AND. PRODUCT(I,2).EQ."e- ") .OR. &
      &  (PRODUCT(I,1).EQ."e-  " .AND. PRODUCT(I,2).EQ."C+ "))) nRCI=I
      IF((REACTANT(I,1).EQ."S  " .AND. REACTANT(I,2).EQ."PHOTON" .AND. REACTANT(I,3).EQ." ") .AND. &
      & ((PRODUCT(I,1).EQ."S+  " .AND. PRODUCT(I,2).EQ."e- ") .OR. &
      &  (PRODUCT(I,1).EQ."e-  " .AND. PRODUCT(I,2).EQ."S+ "))) nRSI=I

      IF((ANY(REACTANT(I,:).EQ."H2+ ") .AND. ANY(REACTANT(I,:).EQ."e- ") .AND. REACTANT(I,3).EQ." ") .AND. &
      &  (ANY(PRODUCT(I,:).EQ."H    ") .AND. ANY(PRODUCT(I,:).EQ."H   "))) nR_H2x_1=I
      IF((ANY(REACTANT(I,:).EQ."H3+ ") .AND. ANY(REACTANT(I,:).EQ."e- ") .AND. REACTANT(I,3).EQ." ") .AND. &
      &  (ALL(PRODUCT(I,1:3).EQ."H  ") .AND. ANY(PRODUCT(I,:).EQ."H   "))) nR_H3x_1=I
      IF((ANY(REACTANT(I,:).EQ."H3+ ") .AND. ANY(REACTANT(I,:).EQ."e- ") .AND. REACTANT(I,3).EQ." ") .AND. &
      &  (ANY(PRODUCT(I,:).EQ."H2   ") .AND. ANY(PRODUCT(I,:).EQ."H   "))) nR_H3x_2=I
      IF((ANY(REACTANT(I,:).EQ."H3O+") .AND. ANY(REACTANT(I,:).EQ."e- ") .AND. REACTANT(I,3).EQ." ") .AND. &
      &  (ANY(PRODUCT(I,:).EQ."OH   ") .AND. ANY(PRODUCT(I,:).EQ."H   "))) nR_H3Ox_1=I
      IF((ANY(REACTANT(I,:).EQ."H3O+") .AND. ANY(REACTANT(I,:).EQ."e- ") .AND. REACTANT(I,3).EQ." ") .AND. &
      &  (ANY(PRODUCT(I,:).EQ."OH   ") .AND. ANY(PRODUCT(I,:).EQ."H2  "))) nR_H3Ox_2=I
      IF((ANY(REACTANT(I,:).EQ."H3O+") .AND. ANY(REACTANT(I,:).EQ."e- ") .AND. REACTANT(I,3).EQ." ") .AND. &
      &  (ANY(PRODUCT(I,:).EQ."H2O  ") .AND. ANY(PRODUCT(I,:).EQ."H   "))) nR_H3Ox_3=I
      IF((ANY(REACTANT(I,:).EQ."HCO+") .AND. ANY(REACTANT(I,:).EQ."e- ") .AND. REACTANT(I,3).EQ." ") .AND. &
      &  (ANY(PRODUCT(I,:).EQ."CO   ") .AND. ANY(PRODUCT(I,:).EQ."H   "))) nR_HCOx_1=I
      IF((ANY(REACTANT(I,:).EQ."H2+ ") .AND. ANY(REACTANT(I,:).EQ."H  ") .AND. REACTANT(I,3).EQ." ") .AND. &
      &  (ANY(PRODUCT(I,:).EQ."H2   ") .AND. ANY(PRODUCT(I,:).EQ."H+  "))) nR_H2x_2=I
      IF((ANY(REACTANT(I,:).EQ."He+ ") .AND. ANY(REACTANT(I,:).EQ."H2 ") .AND. REACTANT(I,3).EQ." ") .AND. &
      &  (ANY(PRODUCT(I,:).EQ."He   ") .AND. ANY(PRODUCT(I,:).EQ."H+  "))) nR_Hex_1=I
      IF((ANY(REACTANT(I,:).EQ."He+ ") .AND. ANY(REACTANT(I,:).EQ."H2 ") .AND. REACTANT(I,3).EQ." ") .AND. &
      &  (ANY(PRODUCT(I,:).EQ."He   ") .AND. ANY(PRODUCT(I,:).EQ."H2+ "))) nR_Hex_2=I
      IF((ANY(REACTANT(I,:).EQ."He+ ") .AND. ANY(REACTANT(I,:).EQ."CO ") .AND. REACTANT(I,3).EQ." ") .AND. &
      &  (ANY(PRODUCT(I,:).EQ."He   ") .AND. ANY(PRODUCT(I,:).EQ."C+  "))) nR_Hex_3=I
      IF((ANY(REACTANT(I,:).EQ."He+ ") .AND. ANY(REACTANT(I,:).EQ."CO ") .AND. REACTANT(I,3).EQ." ") .AND. &
      &  (ANY(PRODUCT(I,:).EQ."He   ") .AND. ANY(PRODUCT(I,:).EQ."O   "))) nR_Hex_4=I

!     Check for duplicate reactions and set the DUPLICATE counter to the appropriate value.
!     Adjust the minimum temperature of each reaction so that temperature ranges are adjacent.
      IF(I.GT.1) THEN
         IF(REACTANT(I,1).EQ.REACTANT(I-1,1) .AND. REACTANT(I,2).EQ.REACTANT(I-1,2) .AND. &
          & REACTANT(I,3).EQ.REACTANT(I-1,3) .AND. &
          & PRODUCT(I,1).EQ.PRODUCT(I-1,1) .AND. PRODUCT(I,2).EQ.PRODUCT(I-1,2) .AND. &
          & PRODUCT(I,3).EQ.PRODUCT(I-1,3) .AND. PRODUCT(I,4).EQ.PRODUCT(I-1,4)) THEN
            IF(DUPLICATE(I-1).EQ.0) DUPLICATE(I-1)=1
            DUPLICATE(I)=DUPLICATE(I-1)+1
            RTMIN(I)=RTMAX(I-1)
         ELSE
            DUPLICATE(I)=0
         END IF
      ELSE
         DUPLICATE(I)=0
      END IF

!     Check for large negative gamma values as they can cause problems when
!     calculating abundances. Produce a warning message for any that occur.
      IF(GAMMA(I).LT.-1.0D2) THEN
         WRITE(10,"('NOTE: Large negative gamma factor for rate',I5,' (',F8.1,')')") INDEX,GAMMA(I)
      END IF

   END DO
   CLOSE(1)

   RETURN
END SUBROUTINE READ_REACTIONS
!=======================================================================
