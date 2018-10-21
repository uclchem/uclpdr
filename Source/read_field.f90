!=======================================================================
!
!  Read the external FUV and X-ray radiation field SEDs incident along
!  each ray. The specified file is assumed to contain entries in comma
!  separated values (CSV) format, removing the need for file-dependent
!  FORMAT statements. Each line in the file describes the unattenuated
!  FUV flux as a function of wavelength or X-ray flux as a function of
!  keV that is incident along the specified ray ID (the first value on
!  each line).
!
!-----------------------------------------------------------------------
SUBROUTINE READ_FIELD(FILENAME,NPART,NRAYS,PARTICLE)

   USE HEALPIX_TYPES
   USE PARTICLE_MODULE
   USE NUM2STR_FUNCTION

   IMPLICIT NONE

   CHARACTER(LEN=*),    INTENT(IN)    :: FILENAME
   INTEGER(KIND=I4B),   INTENT(IN)    :: NPART,NRAYS
   TYPE(PARTICLE_TYPE), INTENT(INOUT) :: PARTICLE(*)

   INTEGER(KIND=I4B) :: I,P,RAY_ID,IER
   REAL(KIND=DP) :: DRAINE_FIELD,XRAY_FIELD

!  Open the input file
   OPEN(UNIT=1,FILE='Input/'//FILENAME,IOSTAT=IER,ACTION='READ',STATUS='OLD')
   READ(1,*,IOSTAT=IER) ! Skip the first line of comments (column headers, etc.)

!  Produce an error message if the file does not exist (or cannot be opened for whatever reason)
   IF(IER.NE.0) THEN
      WRITE(6,*) 'ERROR! Cannot open external radiation field file ',TRIM(FILENAME),' for input'
      WRITE(6,*)
      CLOSE(1)
      STOP
   END IF

!  Read the FUV radiation field incident along the specified ray(s)
   DO I=1,NRAYS
      READ(1,*,IOSTAT=IER,END=1) RAY_ID,DRAINE_FIELD

!     Exit the loop if the next comment line (describing the X-ray field) has been reached
      IF(IER.NE.0) EXIT

!     Check that the ray ID is valid and produce an error message if not
      IF(RAY_ID.LT.0 .OR. RAY_ID.GE.NRAYS) THEN
         WRITE(6,*) 'ERROR! Invalid ray ID (',TRIM(NUM2STR(RAY_ID)),') in external radiation file'
         WRITE(6,*) 'Ray ID values must be in the range 0 <= ID < Nrays'
         WRITE(6,*)
         CLOSE(1)
         STOP
      END IF

!     Assign the value to each particle
      DO P=1,NPART
         PARTICLE(P)%FUV_SURFACE(RAY_ID) = DRAINE_FIELD
      END DO
   END DO

!  Skip the next line of comments describing the X-ray field
   READ(1,*,IOSTAT=IER,END=1)

!  Read the X-ray radiation field incident along the specified ray(s)
   DO I=1,NRAYS
      READ(1,*,IOSTAT=IER,END=1) RAY_ID,XRAY_FIELD

!     Produce an error message if an unexpected data format is encountered
      IF(IER.NE.0) THEN
         WRITE(6,*) 'ERROR! Unexpected data format in external radiation file'
         WRITE(6,*)
         CLOSE(1)
         STOP
      END IF

!     Check that the ray ID is valid and produce an error message if not
      IF(RAY_ID.LT.0 .OR. RAY_ID.GE.NRAYS) THEN
         WRITE(6,*) 'ERROR! Invalid ray ID (',TRIM(NUM2STR(RAY_ID)),') in external radiation file'
         WRITE(6,*) 'Ray ID values must be in the range 0 <= ID < Nrays'
         WRITE(6,*)
         CLOSE(1)
         STOP
      END IF

!     Assign the value to each particle
      DO P=1,NPART
         PARTICLE(P)%XRAY_SURFACE(RAY_ID) = XRAY_FIELD
      END DO
   END DO

!  Clode the file once the end is reached
1  CLOSE(1)

   RETURN
END SUBROUTINE READ_FIELD
!=======================================================================
