!=======================================================================
!
!  Return the number of non-overlapping occurrences
!  of a given substring within the specified string.
!
!-----------------------------------------------------------------------
FUNCTION COUNT_SUBSTRING(STRING,SUBSTRING) RESULT(COUNT)

   IMPLICIT NONE

   CHARACTER(LEN=*), INTENT(IN) :: STRING,SUBSTRING
   INTEGER :: I,COUNT,POSITION

   COUNT=0
   IF(LEN(SUBSTRING).EQ.0) RETURN

   I=1
   DO 
      POSITION=INDEX(STRING(I:),SUBSTRING)
      IF(POSITION.EQ.0) RETURN
      COUNT=COUNT+1
      I=I+POSITION+LEN(SUBSTRING)
   END DO

END FUNCTION COUNT_SUBSTRING
!=======================================================================
