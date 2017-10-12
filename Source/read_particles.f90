!=======================================================================
!
!  Read in the cloud particles, their coordinates, gas number densities,
!  gas and dust temperatures, local FUV fluxes and particle types (i.e.,
!  ionized, PDR or dark cloud). The specified file is assumed to contain
!  entries in comma separated values (CSV) format, removing the need for
!  file-dependent FORMAT statements.
!
!-----------------------------------------------------------------------
SUBROUTINE READ_PARTICLES(FILENAME,NPART,PARTICLE)

   USE HEALPIX_TYPES
   USE PARTICLE_MODULE

   IMPLICIT NONE

   CHARACTER(LEN=*),    INTENT(IN)    :: FILENAME
   INTEGER(KIND=I4B),   INTENT(IN)    :: NPART
   TYPE(PARTICLE_TYPE), INTENT(INOUT) :: PARTICLE(*)

   INTEGER(KIND=I4B) :: P,INDEX,IER
   CHARACTER(LEN=1)  :: TYPE

!  Open the input file
   OPEN(UNIT=1,FILE='Input/'//FILENAME,IOSTAT=IER,ACTION='READ',STATUS='OLD')
   READ(1,*,IOSTAT=IER) ! Skip the first line of comments (column headers, etc.)

!  Produce an error message if the file does not exist (or cannot be opened for whatever reason)
   IF(IER.NE.0) THEN
      WRITE(6,*) 'ERROR! Cannot open cloud particle file ',TRIM(FILENAME),' for input'
      WRITE(6,*)
      CLOSE(1)
      STOP
   ENDIF

!  Initialize the arrays
   DO P=1,NPART
      PARTICLE(P)%COORDINATES=0.0D0
      PARTICLE(P)%GAS_DENSITY=0.0D0
      PARTICLE(P)%DUST_DENSITY=0.0D0
      PARTICLE(P)%GAS_TEMPERATURE=0.0D0
      PARTICLE(P)%DUST_TEMPERATURE=0.0D0
      PARTICLE(P)%FUV_FLUX=0.0D0
      PARTICLE(P)%XRAY_FLUX=0.0D0
      PARTICLE(P)%TOTAL_COLUMN=0.0D0
      PARTICLE(P)%AV=0.0D0
      PARTICLE(P)%FUV_SURFACE=0.0D0
      PARTICLE(P)%XRAY_SURFACE=0.0D0
      PARTICLE(P)%TYPE=0
   ENDDO

   DO P=1,NPART

!     Read the particle coordinates, gas number density, gas and dust temperatures, FUV flux and particle type
      READ(1,*,IOSTAT=IER) INDEX,PARTICLE(P)%COORDINATES,PARTICLE(P)%GAS_DENSITY, &
      & PARTICLE(P)%GAS_TEMPERATURE,PARTICLE(P)%DUST_TEMPERATURE,PARTICLE(P)%FUV_FLUX,TYPE

!     Produce an error message if an unexpected data format is encountered
      IF(IER.NE.0) THEN
         WRITE(6,*) 'ERROR! Unexpected data format on line ',P+1,' of cloud particle file'
         WRITE(6,*)
         CLOSE(1)
         STOP
      ENDIF

!     Calculate the dust number density based on a fixed dust-to-gas mass ratio
      PARTICLE(P)%DUST_DENSITY=2.0D-12*PARTICLE(P)%GAS_DENSITY

!     Assign the particle type (Ionized=1, PDR=2, Dark=3) and check for invalid
!     entries. Produce a warning message if any occur and assign them as TYPE=0
      IF(TYPE.EQ."I" .OR. TYPE.EQ."i") THEN
         PARTICLE(P)%TYPE=1
      ELSE IF(TYPE.EQ."P" .OR. TYPE.EQ."p") THEN
         PARTICLE(P)%TYPE=2
      ELSE IF(TYPE.EQ."D" .OR. TYPE.EQ."d") THEN
         PARTICLE(P)%TYPE=3
      ELSE
         PARTICLE(P)%TYPE=0
         WRITE(6,"(' WARNING! Unrecognized type for particle',I10,' (type=',A,')')") INDEX,TYPE
         WRITE(10,"('WARNING! Unrecognized type for particle',I10,' (type=',A,')')") INDEX,TYPE
      ENDIF

   ENDDO
   CLOSE(1)

   RETURN
END SUBROUTINE READ_PARTICLES
!=======================================================================
