PROGRAM CONVERT_PARTICLES

   IMPLICIT NONE

   INTEGER :: I,J,K,IER
   INTEGER :: NPARTICLE
   REAL(KIND=8) :: T_GAS,T_DUST
   REAL(KIND=8), ALLOCATABLE :: POSITION(:,:),DENSITY(:)
   CHARACTER(LEN=512) :: FILENAME
   CHARACTER :: TYPE

   T_GAS  = 50.0D0
   T_DUST = 20.0D0

!  Read the input file name from the command line arguments
   CALL GETARG(1, FILENAME)

!  Open the file
   OPEN(UNIT=1,FILE=FILENAME,IOSTAT=IER,ACTION='READ',STATUS='OLD')

!  Produce an error message if the file does not exist (or cannot be opened for whatever reason)
   IF(IER.NE.0) THEN
      WRITE(6,*) 'ERROR! Cannot open file ',TRIM(FILENAME)
      WRITE(6,*)
      CLOSE(1)
      STOP
   ENDIF

!  Count the number of lines in the file
   NPARTICLE=0
   DO
      READ(1,*,IOSTAT=IER)
      IF(IER.NE.0) EXIT
      NPARTICLE=NPARTICLE+1
   ENDDO

!  Close the file
   CLOSE(1)

   ALLOCATE(POSITION(1:NPARTICLE,1:3),DENSITY(1:NPARTICLE))

!  Reopen the file, read in the particle positions and calculate their densities
   OPEN(UNIT=1,FILE=FILENAME,IOSTAT=IER,ACTION='READ',STATUS='OLD')
   DO I=1,NPARTICLE
      READ(1,*) POSITION(I,:)
      POSITION(I,:) = POSITION(I,:)*3.08567758E18 ! Convert from pc to cm
      DENSITY(I) = DENSITY_PROFILE(POSITION(I,:))
   END DO

!  Close the file
   CLOSE(1)

!  Open and write the output file
   OPEN(UNIT=1,FILE='cloud.dat',STATUS='REPLACE')
   WRITE(1,"(A)") 'Particle # ,      x (cm)      ,      y (cm)      ,      z (cm)      ,  n_H (c^m-3)  ,   T_gas (K)   ,  T_dust (K)   , chi (Draine)  , Type (I|P|D)'
   DO I=1,NPARTICLE
      TYPE = 'P'
      WRITE(1,"(I10,SP,3(' ,',ES17.9),SS,4(' ,',ES14.7),' ,',A2)") I,POSITION(I,:),DENSITY(I),T_GAS,T_DUST,0,TYPE
   END DO
   CLOSE(1)

   STOP

CONTAINS

FUNCTION DENSITY_PROFILE(POSITION) RESULT(DENSITY)

   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: POSITION(:)
   REAL(KIND=8) :: DENSITY

!!$   DENSITY = 1.0D3
   DENSITY = 1.0D4
!!$   DENSITY = 10.0D0**5.5

END FUNCTION DENSITY_PROFILE

END PROGRAM CONVERT_PARTICLES
