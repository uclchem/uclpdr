!=======================================================================
MODULE FUV_FUNCTIONS

   INTERFACE

      SUBROUTINE READ_AVAILABLE_CROSS_SECTIONS(FILENAME,NXSEC,CROSS_SECTION)
         USE HEALPIX_TYPES
         USE CROSS_SECTION_MODULE
         IMPLICIT NONE
         CHARACTER(LEN=*),         INTENT(IN)  :: FILENAME
         INTEGER(KIND=I4B),        INTENT(OUT) :: NXSEC
         TYPE(CROSS_SECTION_TYPE), ALLOCATABLE, INTENT(OUT) :: CROSS_SECTION(:)
      END SUBROUTINE READ_AVAILABLE_CROSS_SECTIONS

      FUNCTION PHOTOREACTION_RATE(CHI,AV,CROSS_SECTION) RESULT(TOTAL_RATE)
         USE HEALPIX_TYPES
         USE CROSS_SECTION_MODULE
         IMPLICIT NONE
         REAL(KIND=DP),            INTENT(IN) :: CHI,AV
         TYPE(CROSS_SECTION_TYPE), INTENT(IN) :: CROSS_SECTION
         REAL(KIND=DP) :: TOTAL_RATE
      END FUNCTION PHOTOREACTION_RATE

      FUNCTION DRAINE_FIELD_FLUX(CHI,AV,WAVELENGTH,NPOINT) RESULT(FLUX)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
         REAL(KIND=DP),     INTENT(IN)  :: CHI,AV
         REAL(KIND=DP),     INTENT(IN)  :: WAVELENGTH(1:NPOINT)
         REAL(KIND=DP) :: FLUX(1:NPOINT)
      END FUNCTION DRAINE_FIELD_FLUX

      FUNCTION SCATTER(AV,LAMBDA)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         REAL(KIND=DP), INTENT(IN) :: AV,LAMBDA
         REAL(KIND=DP) :: SCATTER
      END FUNCTION SCATTER

      FUNCTION INTEGRATE(F,X,N) RESULT(F_DX)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         REAL(KIND=DP) :: F_DX
         INTEGER(KIND=I4B), INTENT(IN) :: N
         REAL(KIND=DP),     INTENT(IN) :: F(1:N),X(1:N)
      END FUNCTION INTEGRATE

   END INTERFACE

END MODULE FUV_FUNCTIONS
!=======================================================================

!=======================================================================
SUBROUTINE READ_AVAILABLE_CROSS_SECTIONS(FILENAME,REACTANT,PRODUCT, &
                                       & NREAC,NXSEC,CROSS_SECTION)

   USE HEALPIX_TYPES
   USE NUM2STR_FUNCTION
   USE CROSS_SECTION_MODULE

   IMPLICIT NONE

   CHARACTER(LEN=*),         INTENT(IN)  :: FILENAME
   INTEGER(KIND=I4B),        INTENT(IN)  :: NREAC,NXSEC
   CHARACTER(LEN=10),        INTENT(IN)  :: REACTANT(:,:),PRODUCT(:,:)
   TYPE(CROSS_SECTION_TYPE), INTENT(OUT) :: CROSS_SECTION(:)

   INTEGER(KIND=I4B) :: I,N,IER
   CHARACTER(LEN=10) :: CURRENT_REACTANT(1:3),CURRENT_PRODUCT(1:4)

!  Open the input file
   OPEN(UNIT=1,FILE=FILENAME,IOSTAT=IER,ACTION='READ',STATUS='OLD')

!  Produce an error message if the file does not exist (or cannot be opened for whatever reason)
   IF(IER.NE.0) THEN
      WRITE(6,*) 'ERROR! Cannot open photoreaction data file ',TRIM(FILENAME),' for input'
      WRITE(6,*)
      CLOSE(1)
      STOP
   END IF

!  Initialize the photoreaction cross section properties
   DO N=1,NXSEC
      CROSS_SECTION(N)%SPECIES=""
      CROSS_SECTION(N)%FILENAME=""
      CROSS_SECTION(N)%INDEX=0
      CROSS_SECTION(N)%NLINE=0
      CROSS_SECTION(N)%NCONT=0
      CROSS_SECTION(N)%SCALING_FACTOR=0.0D0
   END DO

   DO N=1,NXSEC

!     Initialize the current reactant and product names
      CURRENT_REACTANT = "" ; CURRENT_PRODUCT = ""

!     Read the reactant and product names, rate scaling factor and the associated filename for the cross section data
      READ(1,*,IOSTAT=IER) CURRENT_REACTANT,CURRENT_PRODUCT,CROSS_SECTION(N)%SCALING_FACTOR,CROSS_SECTION(N)%FILENAME

!     Produce an error message if an unexpected data format is encountered
      IF(IER.NE.0) THEN
         WRITE(6,*) 'ERROR! Unexpected data format on line ',TRIM(NUM2STR(N)),' of photoreaction data file ',TRIM(FILENAME)
         WRITE(6,*)
         CLOSE(1)
         STOP
      END IF

!     Determine the reaction index of the photoreaction, if it is present in the chemical network
      DO I=1,NREAC
         IF(ANY(REACTANT(I,:).EQ.CURRENT_REACTANT(1))  .AND. ANY(REACTANT(I,:).EQ.CURRENT_REACTANT(2))) THEN
            IF(ANY(PRODUCT(I,:).EQ.CURRENT_PRODUCT(1)) .AND. ANY(PRODUCT(I,:).EQ.CURRENT_PRODUCT(2)) .AND. &
             & ANY(PRODUCT(I,:).EQ.CURRENT_PRODUCT(3)) .AND. ANY(PRODUCT(I,:).EQ.CURRENT_PRODUCT(4))) THEN
               CROSS_SECTION(N)%SPECIES = REACTANT(I,1)
               CROSS_SECTION(N)%INDEX = I
            END IF
         END IF
      END DO

   END DO
   CLOSE(1)

   CALL READ_CROSS_SECTIONS(NXSEC,CROSS_SECTION)

   RETURN
END SUBROUTINE READ_AVAILABLE_CROSS_SECTIONS
!=======================================================================

!=======================================================================
FUNCTION PHOTOREACTION_RATE(CHI,AV,CROSS_SECTION) RESULT(TOTAL_RATE)

   USE HEALPIX_TYPES
   USE FUV_FUNCTIONS
   USE CROSS_SECTION_MODULE

   IMPLICIT NONE

   REAL(KIND=DP),            INTENT(IN) :: CHI,AV
   TYPE(CROSS_SECTION_TYPE), INTENT(IN) :: CROSS_SECTION

   REAL(KIND=DP) :: TOTAL_RATE

   REAL(KIND=DP) :: LINE_RATE
   REAL(KIND=DP) :: CONTINUUM_RATE
   REAL(KIND=DP), ALLOCATABLE :: FLUX(:)

!  Initialize the line and continuum absorption photoreaction rates
   LINE_RATE = 0.0D0
   CONTINUUM_RATE = 0.0D0

   IF(CROSS_SECTION%NLINE.GT.0) THEN
      ALLOCATE(FLUX(1:CROSS_SECTION%NLINE))
      FLUX = DRAINE_FIELD_FLUX(CHI,AV,CROSS_SECTION%LINE_WAVELENGTH,CROSS_SECTION%NLINE)
      LINE_RATE = SUM(FLUX*CROSS_SECTION%LINE_CROSS_SECTION)
      DEALLOCATE(FLUX)
   END IF

   IF(CROSS_SECTION%NCONT.GT.0) THEN
      ALLOCATE(FLUX(1:CROSS_SECTION%NCONT))
      FLUX = DRAINE_FIELD_FLUX(CHI,AV,CROSS_SECTION%CONTINUUM_WAVELENGTH,CROSS_SECTION%NCONT)
      CONTINUUM_RATE = INTEGRATE(FLUX*CROSS_SECTION%CONTINUUM_CROSS_SECTION, &
                               & CROSS_SECTION%CONTINUUM_WAVELENGTH,CROSS_SECTION%NCONT)
      DEALLOCATE(FLUX)
   END IF

   TOTAL_RATE = LINE_RATE + CONTINUUM_RATE

   RETURN
END FUNCTION PHOTOREACTION_RATE
!=======================================================================

!=======================================================================
FUNCTION DRAINE_FIELD_FLUX(CHI,AV,WAVELENGTH,NPOINT) RESULT(FLUX)

   USE HEALPIX_TYPES
   USE FUNCTIONS_MODULE

   IMPLICIT NONE

   INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
   REAL(KIND=DP),     INTENT(IN)  :: CHI,AV
   REAL(KIND=DP),     INTENT(IN)  :: WAVELENGTH(1:NPOINT)

   REAL(KIND=DP) :: FLUX(1:NPOINT)

   INTEGER(KIND=I4B) :: I

!  Initialize the fluxes
   FLUX = 0.0D0

   DO I=1,NPOINT
      IF(WAVELENGTH(I).LT.911.75D0) THEN
         FLUX(I) = 0.0D0
      ELSE IF(WAVELENGTH(I).LT.2000.0D0) THEN
         FLUX(I) = 3.202833D15/(WAVELENGTH(I)**3) - 5.154206D18/(WAVELENGTH(I)**4) + 2.054625D21/(WAVELENGTH(I)**5)
      ELSE
         FLUX(I) = 7.315D2*WAVELENGTH(I)**0.7D0
      END IF
      FLUX(I) = CHI*FLUX(I)*SCATTER(AV,WAVELENGTH(I))
   END DO

   RETURN
END FUNCTION DRAINE_FIELD_FLUX
!=======================================================================
