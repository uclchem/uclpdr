!=======================================================================
!
!  Read in the chemical species, their initial fractional abundances
!  and their molecular masses. The specified file is assumed to contain
!  entries in comma separated values (CSV) format. This is in line with
!  the new Rate06 standard, removing the need for file-dependent FORMAT
!  statements.
!
!-----------------------------------------------------------------------
SUBROUTINE READ_SPECIES(FILENAME,NSPEC,SPECIES,INITIAL_ABUNDANCE,MOLECULAR_MASS)

   USE HEALPIX_TYPES
   USE GLOBAL_MODULE
   USE NUM2STR_FUNCTION

   IMPLICIT NONE

   CHARACTER(LEN=*),  INTENT(IN)  :: FILENAME
   INTEGER(KIND=I4B), INTENT(IN)  :: NSPEC
   CHARACTER(LEN=10), INTENT(OUT) :: SPECIES(1:NSPEC)
   REAL(KIND=DP),     INTENT(OUT) :: INITIAL_ABUNDANCE(1:NSPEC)
   REAL(KIND=DP),     INTENT(OUT) :: MOLECULAR_MASS(1:NSPEC)

   INTEGER(KIND=I4B) :: I,INDEX,IER

!  Open the input file
   OPEN(UNIT=1,FILE=FILENAME,IOSTAT=IER,ACTION='READ',STATUS='OLD')

!  Produce an error message if the file does not exist (or cannot be opened for whatever reason)
   IF(IER.NE.0) THEN
      WRITE(6,*) 'ERROR! Cannot open chemical species file ',TRIM(FILENAME),' for input'
      WRITE(6,*)
      CLOSE(1)
      STOP
   ENDIF

!  Initialize the arrays
   SPECIES=""
   INITIAL_ABUNDANCE=0.0D0
   MOLECULAR_MASS=0.0D0

!  Initialize the species index numbers. If they are not assigned
!  subsequently any attempt to access those species will generate
!  an error and the code will crash. This is a useful bug catch.
   nH=0 ; nHx=0 ; nD=0 ; nDx=0 ; nH2=0 ; nH2x=0 ; nH3x=0 ; nHD=0
   nH2Dx=0 ; nHe=0 ; nHex=0 ; nC=0 ; nCx=0 ; nN=0 ; nNx=0 ; nO=0
   nOx=0 ; nF=0 ; nFx=0 ; nNa=0 ; nNax=0 ; nMg=0 ; nMgx=0 ; nSi=0
   nSix=0 ; nS=0 ; nSx=0 ; nCl=0 ; nClx=0 ; nCa=0 ; nCax=0 ; nCaxx=0
   nFe=0 ; nFex=0 ; nCH=0 ; nCHx=0 ; nCH2=0 ; nCH2x=0 ; nCH3x=0
   nOH=0 ; nH2O=0 ; nH2Ox=0 ; nH3Ox=0 ; nNH=0 ; nNH2=0 ; nNH3=0
   nCO=0 ; nHCOx=0 ; nCS=0 ; nH2v=0 ; nPAH=0 ; nPAHx=0 ; nPAHm=0
   nelect=0

   DO I=1,NSPEC

!     Read the species name, initial abundance and molecular mass
      READ(1,*,IOSTAT=IER) INDEX,SPECIES(I),INITIAL_ABUNDANCE(I),MOLECULAR_MASS(I)

!     Produce an error message if an unexpected data format is encountered
      IF(IER.NE.0) THEN
         WRITE(6,*) 'ERROR! Unexpected data format on line ',TRIM(NUM2STR(I)),' of chemical species file'
         WRITE(6,*)
         CLOSE(1)
         STOP
      ENDIF

!     Assign the various index numbers to their correct species
      IF(SPECIES(I).EQ."H         ") nH      = I
      IF(SPECIES(I).EQ."H+        ") nHx     = I
      IF(SPECIES(I).EQ."D         ") nD      = I
      IF(SPECIES(I).EQ."D+        ") nDx     = I
      IF(SPECIES(I).EQ."H2        ") nH2     = I
      IF(SPECIES(I).EQ."H2*       ") nH2v    = I
      IF(SPECIES(I).EQ."H2+       ") nH2x    = I
      IF(SPECIES(I).EQ."H3+       ") nH3x    = I
      IF(SPECIES(I).EQ."HD        ") nHD     = I
      IF(SPECIES(I).EQ."H2D+      ") nH2Dx   = I
      IF(SPECIES(I).EQ."He        ") nHe     = I
      IF(SPECIES(I).EQ."HE        ") nHe     = I
      IF(SPECIES(I).EQ."He+       ") nHex    = I
      IF(SPECIES(I).EQ."HE+       ") nHex    = I
      IF(SPECIES(I).EQ."C         ") nC      = I
      IF(SPECIES(I).EQ."C+        ") nCx     = I
      IF(SPECIES(I).EQ."N         ") nN      = I
      IF(SPECIES(I).EQ."N+        ") nNx     = I
      IF(SPECIES(I).EQ."O         ") nO      = I
      IF(SPECIES(I).EQ."O+        ") nOx     = I
      IF(SPECIES(I).EQ."F         ") nF      = I
      IF(SPECIES(I).EQ."F+        ") nFx     = I
      IF(SPECIES(I).EQ."Na        ") nNa     = I
      IF(SPECIES(I).EQ."NA        ") nNa     = I
      IF(SPECIES(I).EQ."Na+       ") nNax    = I
      IF(SPECIES(I).EQ."NA+       ") nNax    = I
      IF(SPECIES(I).EQ."Mg        ") nMg     = I
      IF(SPECIES(I).EQ."MG        ") nMg     = I
      IF(SPECIES(I).EQ."Mg+       ") nMgx    = I
      IF(SPECIES(I).EQ."MG+       ") nMgx    = I
      IF(SPECIES(I).EQ."Si        ") nSi     = I
      IF(SPECIES(I).EQ."SI        ") nSi     = I
      IF(SPECIES(I).EQ."Si+       ") nSix    = I
      IF(SPECIES(I).EQ."SI+       ") nSix    = I
      IF(SPECIES(I).EQ."S         ") nS      = I
      IF(SPECIES(I).EQ."S+        ") nSx     = I
      IF(SPECIES(I).EQ."Cl        ") nCl     = I
      IF(SPECIES(I).EQ."CL        ") nCl     = I
      IF(SPECIES(I).EQ."Cl+       ") nClx    = I
      IF(SPECIES(I).EQ."CL+       ") nClx    = I
      IF(SPECIES(I).EQ."Ca        ") nCa     = I
      IF(SPECIES(I).EQ."CA        ") nCa     = I
      IF(SPECIES(I).EQ."Ca+       ") nCax    = I
      IF(SPECIES(I).EQ."CA+       ") nCax    = I
      IF(SPECIES(I).EQ."Ca++      ") nCaxx   = I
      IF(SPECIES(I).EQ."CA++      ") nCaxx   = I
      IF(SPECIES(I).EQ."Fe        ") nFe     = I
      IF(SPECIES(I).EQ."FE        ") nFe     = I
      IF(SPECIES(I).EQ."Fe+       ") nFex    = I
      IF(SPECIES(I).EQ."FE+       ") nFex    = I
      IF(SPECIES(I).EQ."CH        ") nCH     = I
      IF(SPECIES(I).EQ."CH+       ") nCHx    = I
      IF(SPECIES(I).EQ."CH2       ") nCH2    = I
      IF(SPECIES(I).EQ."CH2+      ") nCH2x   = I
      IF(SPECIES(I).EQ."CH3+      ") nCH3x   = I
      IF(SPECIES(I).EQ."OH        ") nOH     = I
      IF(SPECIES(I).EQ."H2O       ") nH2O    = I
      IF(SPECIES(I).EQ."H2O+      ") nH2Ox   = I
      IF(SPECIES(I).EQ."H3O+      ") nH3Ox   = I
      IF(SPECIES(I).EQ."NH        ") nNH     = I
      IF(SPECIES(I).EQ."NH2       ") nNH2    = I
      IF(SPECIES(I).EQ."NH3       ") nNH3    = I
      IF(SPECIES(I).EQ."CO        ") nCO     = I
      IF(SPECIES(I).EQ."HCO+      ") nHCOx   = I
      IF(SPECIES(I).EQ."CS        ") nCS     = I
      IF(SPECIES(I).EQ."PAH       ") nPAH    = I
      IF(SPECIES(I).EQ."PAH0      ") nPAH    = I
      IF(SPECIES(I).EQ."PAH+      ") nPAHx   = I
      IF(SPECIES(I).EQ."PAH-      ") nPAHm   = I
      IF(SPECIES(I).EQ."e-        ") nelect  = I
      IF(SPECIES(I).EQ."ELECTR    ") nelect  = I

   ENDDO
   CLOSE(1)

!  Check that the final species in the file is electron
!  Print a warning message to screen and logfile if not
   IF(SPECIES(NSPEC).NE."e-") THEN
      WRITE(6,*) 'WARNING! Last entry in species file is not e-'
      WRITE(10,*)'WARNING! Last entry in species file is not e-'
   ENDIF

!  Check that the total hydrogen nuclei abundance adds up to 1
!  If not, modify the abundance of H2 (only consider H, H+ & H3+)
   IF((INITIAL_ABUNDANCE(nH)+INITIAL_ABUNDANCE(nHx)+2.0D0*INITIAL_ABUNDANCE(nH2)).NE.1.0D0) THEN
      INITIAL_ABUNDANCE(nH2)=0.5D0*(1.0D0-INITIAL_ABUNDANCE(nH)-INITIAL_ABUNDANCE(nHx)-INITIAL_ABUNDANCE(nH3x))
   ENDIF

!  Calculate the intial electron abundance, if not specified, as the sum of the metal ion abundances
   IF(INITIAL_ABUNDANCE(nelect).LE.0.0D0) THEN
      INITIAL_ABUNDANCE(nelect)=0.0D0
      IF(nHx.NE.0)  INITIAL_ABUNDANCE(nelect)=INITIAL_ABUNDANCE(nelect)+INITIAL_ABUNDANCE(nHx)
      IF(nCx.NE.0)  INITIAL_ABUNDANCE(nelect)=INITIAL_ABUNDANCE(nelect)+INITIAL_ABUNDANCE(nCx)
      IF(nSx.NE.0)  INITIAL_ABUNDANCE(nelect)=INITIAL_ABUNDANCE(nelect)+INITIAL_ABUNDANCE(nSx)
      IF(nNax.NE.0) INITIAL_ABUNDANCE(nelect)=INITIAL_ABUNDANCE(nelect)+INITIAL_ABUNDANCE(nNax)
      IF(nMgx.NE.0) INITIAL_ABUNDANCE(nelect)=INITIAL_ABUNDANCE(nelect)+INITIAL_ABUNDANCE(nMgx)
      IF(nSix.NE.0) INITIAL_ABUNDANCE(nelect)=INITIAL_ABUNDANCE(nelect)+INITIAL_ABUNDANCE(nSix)
      IF(nClx.NE.0) INITIAL_ABUNDANCE(nelect)=INITIAL_ABUNDANCE(nelect)+INITIAL_ABUNDANCE(nClx)
      IF(nCax.NE.0) INITIAL_ABUNDANCE(nelect)=INITIAL_ABUNDANCE(nelect)+INITIAL_ABUNDANCE(nCax)
      IF(nFex.NE.0) INITIAL_ABUNDANCE(nelect)=INITIAL_ABUNDANCE(nelect)+INITIAL_ABUNDANCE(nFex)
      IF(nPAHx.NE.0)INITIAL_ABUNDANCE(nelect)=INITIAL_ABUNDANCE(nelect)+INITIAL_ABUNDANCE(nPAHx)
      IF(nPAHm.NE.0)INITIAL_ABUNDANCE(nelect)=INITIAL_ABUNDANCE(nelect)-INITIAL_ABUNDANCE(nPAHm)
   ENDIF

   RETURN
END SUBROUTINE READ_SPECIES
!=======================================================================
