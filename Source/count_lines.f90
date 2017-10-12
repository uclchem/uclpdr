!=======================================================================
!
!  Return the number of lines contained in the specified file.
!
!-----------------------------------------------------------------------
FUNCTION COUNT_LINES(FILENAME) RESULT(LINES)

   USE HEALPIX_TYPES
   IMPLICIT NONE

   INTEGER(KIND=I4B) :: LINES
   CHARACTER(LEN=*), INTENT(IN) :: FILENAME
   INTEGER(KIND=I4B) :: IER

!  Open the file
   OPEN(UNIT=1,FILE=FILENAME,IOSTAT=IER,ACTION='READ',STATUS='OLD')

!  Produce an error message if the file does not exist (or cannot be opened for whatever reason)
   IF(IER.NE.0) THEN
      WRITE(6,*) 'ERROR! Cannot open file ',TRIM(FILENAME),' in function COUNT_LINES'
      WRITE(6,*)
      CLOSE(1)
      STOP
   ENDIF

!  Count the number of lines in the file
   LINES=0
   DO
      READ(1,*,IOSTAT=IER)
      IF(IER.NE.0) EXIT
      LINES=LINES+1
   ENDDO

!  Close the file
   CLOSE(1)

   RETURN
END FUNCTION COUNT_LINES
!=======================================================================
