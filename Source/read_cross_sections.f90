!=======================================================================
!
!  Read in the photoreaction cross sections for each available species,
!  including values for both line and continuum absorption, if present.
!  The specified files are all assumed to contain entries in the format
!  adopted by Ewine van Dishoeck on her photorates website:
!
!     http://home.strw.leidenuniv.nl/~ewine/photo/
!
!  therefore allowing new files to be downloaded directly from the site.
!
!-----------------------------------------------------------------------
SUBROUTINE READ_CROSS_SECTIONS(NXSEC,CROSS_SECTION)

   USE HEALPIX_TYPES
   USE CROSS_SECTION_MODULE

   IMPLICIT NONE

   INTEGER(KIND=I4B), INTENT(IN) :: NXSEC
   TYPE(CROSS_SECTION_TYPE), INTENT(INOUT) :: CROSS_SECTION(*)

   INTEGER(KIND=I4B) :: I,L,N,INDEX,IER
   INTEGER(KIND=I4B) :: NLINE,NCONT

   DO N=1,NXSEC ! Loop over cross section datafiles

!     Open the input file
      OPEN(UNIT=1,FILE='Datafiles/Photoreaction-Rates/'//CROSS_SECTION(N)%FILENAME,IOSTAT=IER,ACTION='READ',STATUS='OLD')
      READ(1,*,IOSTAT=IER)

!     Produce an error message if the file does not exist (or cannot be opened for whatever reason)
      IF(IER.NE.0) THEN
         WRITE(6,*) 'ERROR! Cannot open cross section data file ',TRIM(CROSS_SECTION(N)%FILENAME),' for input'
         WRITE(6,*)
         CLOSE(1)
         STOP
      END IF

!     Read the number of line absorption values and allocate
!     the wavelength and cross section arrays accordingly
      READ(1,*,IOSTAT=IER) NLINE
      IF(NLINE.LT.0) THEN
         WRITE(6,*) 'ERROR! Incorrect number of line absorption entries in cross section data file ',TRIM(CROSS_SECTION(N)%FILENAME),' (NLINE=',NLINE,')'
         CLOSE(1)
         STOP
      END IF
      CROSS_SECTION(N)%NLINE=NLINE
      IF(NLINE.GT.0) THEN
         ALLOCATE(CROSS_SECTION(N)%LINE_WAVELENGTH(1:NLINE))
         ALLOCATE(CROSS_SECTION(N)%LINE_CROSS_SECTION(1:NLINE))
      ELSE
         ALLOCATE(CROSS_SECTION(N)%LINE_WAVELENGTH(1:1))
         ALLOCATE(CROSS_SECTION(N)%LINE_CROSS_SECTION(1:1))
      ENDIF

!     Initialize the wavelengths and cross section values
      CROSS_SECTION(N)%LINE_WAVELENGTH=0.0D0
      CROSS_SECTION(N)%LINE_CROSS_SECTION=0.0D0

!     Read the wavelength (Å) and cross section (cm^2 Å) for each line absorption
      DO L=1,NLINE ! Loop over lines
         READ(1,*,IOSTAT=IER) I,CROSS_SECTION(N)%LINE_WAVELENGTH(I),CROSS_SECTION(N)%LINE_CROSS_SECTION(I)
      END DO ! End of loop over lines

!     Read the number of continuum absorption values and allocate
!     the wavelength and cross section arrays accordingly
      READ(1,*,IOSTAT=IER) NCONT
      READ(1,*,IOSTAT=IER)
      IF(NCONT.LT.0) THEN
         WRITE(6,*) 'ERROR! Incorrect number of continuum absorption entries in cross section data file ',TRIM(CROSS_SECTION(N)%FILENAME),' (NCONT=',NCONT,')'
         CLOSE(1)
         STOP
      END IF
      CROSS_SECTION(N)%NCONT=NCONT
      IF(NCONT.GT.0) THEN
         ALLOCATE(CROSS_SECTION(N)%CONTINUUM_WAVELENGTH(1:NCONT))
         ALLOCATE(CROSS_SECTION(N)%CONTINUUM_CROSS_SECTION(1:NCONT))
      ELSE
         ALLOCATE(CROSS_SECTION(N)%CONTINUUM_WAVELENGTH(1:1))
         ALLOCATE(CROSS_SECTION(N)%CONTINUUM_CROSS_SECTION(1:1))
      ENDIF

!     Initialize the wavelengths and cross section values
      CROSS_SECTION(N)%CONTINUUM_WAVELENGTH=0.0D0
      CROSS_SECTION(N)%CONTINUUM_CROSS_SECTION=0.0D0

!     Read the wavelength (Å) and cross section (cm^2) values for continuous absorption
      DO L=1,NCONT ! Loop over entries
         READ(1,*,IOSTAT=IER) I,CROSS_SECTION(N)%CONTINUUM_WAVELENGTH(I),CROSS_SECTION(N)%CONTINUUM_CROSS_SECTION(I)
      END DO ! End of loop over entries

      CLOSE(1)

   END DO ! End of loop over cross section datafiles
   
   RETURN
END SUBROUTINE READ_CROSS_SECTIONS
!=======================================================================
