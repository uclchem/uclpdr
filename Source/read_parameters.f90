!=======================================================================
!
!  Read in the physical parameters that describe the model simulation.
!  Determine the number of particles in the model cloud, the number of
!  species and reactions in the chemical network, and the total number
!  of available photorate cross sections from the relevant input files.
!
!-----------------------------------------------------------------------
SUBROUTINE READ_PARAMETERS(FILENAME,PARTICLE_FILE,RADIATION_FILE,OUTPUT_PREFIX)

   USE HEALPIX_TYPES
   USE MAIN_MODULE
   USE HEALPIX_MODULE,   ONLY : HEALPIX_LEVEL,THETA_CRITICAL
   USE CHEMISTRY_MODULE, ONLY : RELATIVE_ABUNDANCE_TOLERANCE,ABSOLUTE_ABUNDANCE_TOLERANCE
   USE GLOBAL_MODULE,    ONLY : END_TIME,ZETA,T_CMB,V_TURB,METALLICITY,AV_FAC,UV_FAC,OMEGA,GRAIN_RADIUS
   USE FUNCTIONS_MODULE
   USE NUM2STR_FUNCTION

   IMPLICIT NONE

   CHARACTER(LEN=*), INTENT(IN)  :: FILENAME
   CHARACTER(LEN=*), INTENT(OUT) :: PARTICLE_FILE,RADIATION_FILE,OUTPUT_PREFIX

   INTEGER(KIND=I4B) :: I,N,IER

!  Open the input file
   OPEN(UNIT=1,FILE=FILENAME,IOSTAT=IER,ACTION='READ',STATUS='OLD')

!  Produce an error message if the file does not exist (or cannot be opened for whatever reason)
   IF(IER.NE.0) THEN
      WRITE(6,*) 'ERROR! Cannot open model parameters file ',TRIM(FILENAME),' for input'
      WRITE(6,*)
      CLOSE(1)
      STOP
   END IF

!  Set the line counter
   I=1

!  Read the parameters, converting units where necessary
   READ(1,*,IOSTAT=IER,ERR=1) ; I=I+1 ; READ(1,*,IOSTAT=IER,ERR=1) ; I=I+1 ; READ(1,*,IOSTAT=IER,ERR=1) ; I=I+1
   READ(1,*,IOSTAT=IER,ERR=1) PARTICLE_FILE  ; I=I+1
   READ(1,*,IOSTAT=IER,ERR=1) RADIATION_FILE ; I=I+1
   READ(1,*,IOSTAT=IER,ERR=1) OUTPUT_PREFIX  ; I=I+1
   READ(1,*,IOSTAT=IER,ERR=1) ; I=I+1 ; READ(1,*,IOSTAT=IER,ERR=1) ; I=I+1 ; READ(1,*,IOSTAT=IER,ERR=1) ; I=I+1
   READ(1,*,IOSTAT=IER,ERR=1) NRAYS ; I=I+1
   READ(1,*,IOSTAT=IER,ERR=1) HEALPIX_LEVEL  ; I=I+1 ; IF(NRAYS.EQ.0) NRAYS = 12*4**HEALPIX_LEVEL ! Use HEALPix
   READ(1,*,IOSTAT=IER,ERR=1) THETA_CRITICAL ; I=I+1
   READ(1,*,IOSTAT=IER,ERR=1) ; I=I+1 ; READ(1,*,IOSTAT=IER,ERR=1) ; I=I+1 ; READ(1,*,IOSTAT=IER,ERR=1) ; I=I+1
   READ(1,*,IOSTAT=IER,ERR=1) END_TIME ; I=I+1
   READ(1,*,IOSTAT=IER,ERR=1) ZETA ; I=I+1 ; ZETA = ZETA/1.3D-17 ! Scale by the standard rate
   READ(1,*,IOSTAT=IER,ERR=1) V_TURB ; I=I+1 ; V_TURB = V_TURB*1.0D5 ! Convert from km/s to cm/s
   READ(1,*,IOSTAT=IER,ERR=1) METALLICITY ; I=I+1
   READ(1,*,IOSTAT=IER,ERR=1) T_CMB ; I=I+1
   READ(1,*,IOSTAT=IER,ERR=1) ; I=I+1 ; READ(1,*,IOSTAT=IER,ERR=1) ; I=I+1 ; READ(1,*,IOSTAT=IER,ERR=1) ; I=I+1
   READ(1,*,IOSTAT=IER,ERR=1) AV_FAC ; I=I+1
   READ(1,*,IOSTAT=IER,ERR=1) UV_FAC ; I=I+1
   READ(1,*,IOSTAT=IER,ERR=1) OMEGA  ; I=I+1
   READ(1,*,IOSTAT=IER,ERR=1) GRAIN_RADIUS ; I=I+1 ; GRAIN_RADIUS = GRAIN_RADIUS*1.0D-4 ! Convert from µm to cm
   READ(1,*,IOSTAT=IER,ERR=1) ; I=I+1 ; READ(1,*,IOSTAT=IER,ERR=1) ; I=I+1 ; READ(1,*,IOSTAT=IER,ERR=1) ; I=I+1
   READ(1,*,IOSTAT=IER,ERR=1) GAS_TEMPERATURE_GUESS  ; I=I+1
   READ(1,*,IOSTAT=IER,ERR=1) DUST_TEMPERATURE_GUESS ; I=I+1
   READ(1,*,IOSTAT=IER,ERR=1) ; I=I+1 ; READ(1,*,IOSTAT=IER,ERR=1) ; I=I+1 ; READ(1,*,IOSTAT=IER,ERR=1) ; I=I+1
   READ(1,*,IOSTAT=IER,ERR=1) CHEMISTRY_ITERATIONS   ; I=I+1
   READ(1,*,IOSTAT=IER,ERR=1) TEMPERATURE_ITERATIONS ; I=I+1
   READ(1,*,IOSTAT=IER,ERR=1) POPULATION_ITERATIONS  ; I=I+1
   READ(1,*,IOSTAT=IER,ERR=1) ; I=I+1 ; READ(1,*,IOSTAT=IER,ERR=1) ; I=I+1 ; READ(1,*,IOSTAT=IER,ERR=1) ; I=I+1
   READ(1,*,IOSTAT=IER,ERR=1) RELATIVE_ABUNDANCE_TOLERANCE ; I=I+1
   READ(1,*,IOSTAT=IER,ERR=1) ABSOLUTE_ABUNDANCE_TOLERANCE ; I=I+1
   READ(1,*,IOSTAT=IER,ERR=1) ABUNDANCE_LIMIT ; I=I+1
   READ(1,*,IOSTAT=IER,ERR=1) CHEMISTRY_CONVERGENCE_CRITERION  ; I=I+1 ; CHEMISTRY_CONVERGENCE_CRITERION = CHEMISTRY_CONVERGENCE_CRITERION*1.0D-2 ! Convert from % to fraction
   READ(1,*,IOSTAT=IER,ERR=1) ; I=I+1 ; READ(1,*,IOSTAT=IER,ERR=1) ; I=I+1 ; READ(1,*,IOSTAT=IER,ERR=1) ; I=I+1
   READ(1,*,IOSTAT=IER,ERR=1) POPULATION_LIMIT ; I=I+1
   READ(1,*,IOSTAT=IER,ERR=1) POPULATION_CONVERGENCE_CRITERION ; I=I+1 ; POPULATION_CONVERGENCE_CRITERION = POPULATION_CONVERGENCE_CRITERION*1.0D-2 ! Convert from % to fraction
   READ(1,*,IOSTAT=IER,ERR=1) POPULATION_ITERATION_LIMIT  ; I=I+1
   READ(1,*,IOSTAT=IER,ERR=1) POPULATION_PERCENTAGE_LIMIT ; I=I+1
   READ(1,*,IOSTAT=IER,ERR=1) ; I=I+1 ; READ(1,*,IOSTAT=IER,ERR=1) ; I=I+1 ; READ(1,*,IOSTAT=IER,ERR=1) ; I=I+1
   READ(1,*,IOSTAT=IER,ERR=1) TMIN  ; I=I+1
   READ(1,*,IOSTAT=IER,ERR=1) TMAX  ; I=I+1
   READ(1,*,IOSTAT=IER,ERR=1) TDIFF ; I=I+1
   READ(1,*,IOSTAT=IER,ERR=1) FCRIT ; I=I+1 ; FCRIT = FCRIT*1.0D-2 ! Convert from % to fraction
   READ(1,*,IOSTAT=IER,ERR=1) ; I=I+1 ; READ(1,*,IOSTAT=IER,ERR=1) ; I=I+1 ; READ(1,*,IOSTAT=IER,ERR=1) ; I=I+1
   READ(1,*,IOSTAT=IER,ERR=1) NCOOL ; I=I+1

!  Allocate the coolant species array and read in the filename(s)
   ALLOCATE(COOLANT(1:NCOOL))
   DO N=1,NCOOL ! Loop over coolants
      READ(1,*,IOSTAT=IER,ERR=1) COOLANT(N)%FILENAME ; I=I+1
   END DO

!  Produce an error message if an unexpected data format is encountered
 1 IF(IER.NE.0) THEN
      WRITE(6,*) 'ERROR! Unexpected data format on line ',TRIM(NUM2STR(I)),' of model parameters file'
      WRITE(6,*)
      CLOSE(1)
      STOP
   END IF

   CLOSE(1)

!  Check the specified parameter values are valid
   IF(HEALPIX_LEVEL.LT.0 .OR. HEALPIX_LEVEL.GT.4) THEN
      WRITE(6,*) 'HEALPix level of refinement must be between 0 and 4'
      WRITE(6,*)
      STOP
   END IF
   IF(THETA_CRITICAL.LT.0 .OR. THETA_CRITICAL.GE.PI/2) THEN
      WRITE(6,*) 'HEALPix search angle must be in the range: 0 < theta < π/2'
      WRITE(6,*)
      STOP
   END IF
   IF(END_TIME.LT.0 .OR. END_TIME.GT.1.0D9) THEN
      WRITE(6,*) 'Simulation duration must be in the range: 0 < t < 1E9 yr'
      WRITE(6,*)
      STOP
   END IF
   IF(ZETA.LT.0) THEN
      WRITE(6,*) 'Cosmic-ray ionization rate must be a positive value'
      WRITE(6,*)
      STOP
   END IF

!  Determine the number of particles in the model cloud, the number of
!  species and reactions in the chemical network, and the total number
!  of available photorate cross sections from the relevant input files
   NPART = COUNT_LINES('Input/'//PARTICLE_FILE)-1
   NSPEC = COUNT_LINES('Datafiles/Chemical-Network/species.dat')
   NREAC = COUNT_LINES('Datafiles/Chemical-Network/rates.dat')
   NXSEC = COUNT_LINES('Datafiles/Photoreaction-Rates/photorates.dat')
   NHEAT = 10 ! Fixed number of heating mechanisms considered

   RETURN
END SUBROUTINE READ_PARAMETERS
!=======================================================================
