!=======================================================================
!
!   UCL-PDR: a one-dimensional astrochemistry and radiative transfer
!            code for modelling FUV- and X-ray-illuminated regions.
!
!-----------------------------------------------------------------------
PROGRAM UCL_PDR

#ifdef OPENMP
   USE OMP_LIB
#endif
   USE HEALPIX_TYPES
   USE MAIN_MODULE
   USE CHEMISTRY_MODULE
   USE GLOBAL_MODULE
   USE OUTPUT_MODULE
   USE FUNCTIONS_MODULE
   USE SUBROUTINES_MODULE
   USE NUM2STR_FUNCTION

   IMPLICIT NONE

   INTEGER(KIND=I4B)  :: I,J,P,N ! Loop variables
   INTEGER(KIND=I4B)  :: TEMPERATURE_ITERATION,CHEMISTRY_ITERATION,POPULATION_ITERATION ! Iteration counters
   INTEGER(KIND=I4B)  :: CHEMISTRY_STATUS,LEVEL_POPULATION_STATUS,THERMAL_BALANCE_STATUS ! Status flags
   REAL(KIND=DP)      :: CPU_START,CPU_CHEM,CPU_POP,CPU_END ! CPU times returned by the CPU_TIME subroutine

   CHARACTER(LEN=256) :: PARTICLE_FILE  ! Particle data file name
   CHARACTER(LEN=256) :: RADIATION_FILE ! Radiation field file name
   CHARACTER(LEN=128) :: FILE_PREFIX ! Output file name prefix
   CHARACTER(LEN=8)   :: FILE_SUFFIX ! Output file name suffix

   INTEGER(KIND=I4B)  :: NCHEM ! Number of particles whose abundances are to be calculated in the current iteration
   INTEGER(KIND=I4B), ALLOCATABLE :: DUMMY_INDEX(:) ! Index of each particle whose abundances are to be calculated
   REAL(KIND=DP),     ALLOCATABLE :: DUMMY_ABUNDANCE(:,:),DUMMY_RATE(:,:) ! Temporary arrays to store the relevant
   REAL(KIND=DP),     ALLOCATABLE :: DUMMY_DENSITY(:),DUMMY_TEMPERATURE(:) ! properties of particles to be updated

   REAL(KIND=DP),     ALLOCATABLE :: PERCENTAGE_POPULATION_CONVERGED(:)
   REAL(KIND=DP)                  :: PERCENTAGE_TEMPERATURE_CONVERGED
   REAL(KIND=DP)                  :: PERCENTAGE_CHEMISTRY_CONVERGED

   LOGICAL :: FILE_EXISTS

#ifdef OPENMP
!  Define the variables needed by OpenMP
   INTEGER(KIND=I4B) :: NTHREAD,THREAD ! Number of threads and current thread number
#endif

#ifdef OPENMP
   CPU_START = OMP_GET_WTIME()
#else
   CALL CPU_TIME(CPU_START)
#endif

   CALL PRINT_SPLASHSCREEN
   WRITE(6,*) 'Starting model...'

!  Open and set up the log files
   CALL OPEN_LOGFILES

   WRITE(6,*) 'Reading input files...'

!  Read the model parameters input file
   CALL READ_PARAMETERS('Input/model-parameters.dat',PARTICLE_FILE,RADIATION_FILE,FILE_PREFIX)

!  Allocate the chemical network arrays
   ALLOCATE(SPECIES(1:NSPEC),INITIAL_ABUNDANCE(1:NSPEC),MOLECULAR_MASS(1:NSPEC))
   ALLOCATE(REACTANT(1:NREAC,1:3),PRODUCT(1:NREAC,1:4))
   ALLOCATE(ALPHA(1:NREAC),BETA(1:NREAC),GAMMA(1:NREAC))
   ALLOCATE(RTMIN(1:NREAC),RTMAX(1:NREAC),DUPLICATE(1:NREAC))
   ALLOCATE(CROSS_SECTION(1:NXSEC))

!  Allocate the particle property arrays
   ALLOCATE(PARTICLE(1:NPART))
   DO P=1,NPART
      ALLOCATE(PARTICLE(P)%TOTAL_COLUMN(0:NRAYS-1),PARTICLE(P)%AV(0:NRAYS-1))
      ALLOCATE(PARTICLE(P)%FUV_SURFACE(0:NRAYS-1),PARTICLE(P)%XRAY_SURFACE(0:NRAYS-1))
      ALLOCATE(PARTICLE(P)%RATE(1:NREAC),PARTICLE(P)%ABUNDANCE(1:NSPEC))
      ALLOCATE(PARTICLE(P)%COLUMN_DENSITY(0:NRAYS-1,1:NSPEC))
      ALLOCATE(PARTICLE(P)%COOLING_RATE(1:NCOOL),PARTICLE(P)%HEATING_RATE(1:NHEAT))
      ALLOCATE(PARTICLE(P)%PREVIOUS_ABUNDANCE(1:NSPEC))
   END DO

!  Allocate the dummy arrays to pass to the chemistry routine
   ALLOCATE(DUMMY_INDEX(1:NPART),DUMMY_DENSITY(1:NPART),DUMMY_TEMPERATURE(1:NPART))
   ALLOCATE(DUMMY_ABUNDANCE(1:NSPEC,1:NPART),DUMMY_RATE(1:NREAC,1:NPART))

!  Read the cloud particle properties
   CALL READ_PARTICLES(PARTICLE_FILE,NPART,PARTICLE)

!  Read the external radiation field properties
   CALL READ_FIELD(RADIATION_FILE,NPART,NRAYS,PARTICLE)

!  Read the species of the chemical network, their initial abundances and molecular masses
   CALL READ_SPECIES('Datafiles/Chemical-Network/species.dat',NSPEC,SPECIES,INITIAL_ABUNDANCE,MOLECULAR_MASS)

!  Read the reactions of the chemical network, their reactants, products and Arrhenius equation parameters
   CALL READ_REACTIONS('Datafiles/Chemical-Network/rates.dat',NREAC,REACTANT,PRODUCT,ALPHA,BETA,GAMMA,RTMIN,RTMAX,DUPLICATE)

!  Read the coolant species, their energy level and radiative transition properties and their collisional rate coefficients
   CALL READ_COOLANTS(NCOOL,COOLANT)

!  Read the available photoreaction cross section data
   CALL READ_AVAILABLE_CROSS_SECTIONS('Datafiles/Photoreaction-Rates/photorates.dat',REACTANT,PRODUCT,NREAC,NXSEC,CROSS_SECTION)

!  Assign the correct index number for each coolant species, producing
!  error message if the species is not present in the chemical network
   DO N=1,NCOOL
      DO I=1,NSPEC
!        Ignore "p-" or "o-" prefixes in the coolant name
         IF(TRIM(ADJUSTL(SPECIES(I))).EQ.TRIM(ADJUSTL(COOLANT(N)%NAME(INDEX(COOLANT(N)%NAME,'-')+1:)))) THEN
            COOLANT(N)%INDEX = I
            EXIT
         END IF
      END DO
      IF(COOLANT(N)%INDEX.EQ.0) THEN
         WRITE(6,*) 'ERROR! Coolant species ',TRIM(ADJUSTL(COOLANT(N)%NAME)),' is not present in the chemical network'
         STOP
      END IF
   END DO

!  Allocate the level population arrays
   DO P=1,NPART
      ALLOCATE(PARTICLE(P)%COOLANT(1:NCOOL))
      DO N=1,NCOOL
         PARTICLE(P)%COOLANT(N)%NAME = COOLANT(N)%NAME
         PARTICLE(P)%COOLANT(N)%NLEVEL = COOLANT(N)%NLEVEL
         ALLOCATE(PARTICLE(P)%COOLANT(N)%POPULATION(1:COOLANT(N)%NLEVEL))
         ALLOCATE(PARTICLE(P)%COOLANT(N)%PREVIOUS_POPULATION(1:COOLANT(N)%NLEVEL))
         ALLOCATE(PARTICLE(P)%COOLANT(N)%EMISSIVITY(1:COOLANT(N)%NLEVEL,1:COOLANT(N)%NLEVEL))
         ALLOCATE(PARTICLE(P)%COOLANT(N)%OPACITY(0:NRAYS-1,1:COOLANT(N)%NLEVEL,1:COOLANT(N)%NLEVEL))
#ifdef USE_ALI
         ALLOCATE(PARTICLE(P)%COOLANT(N)%LAMBDA(1:COOLANT(N)%NLEVEL,1:COOLANT(N)%NLEVEL))
#endif
      END DO
   END DO
   ALLOCATE(PERCENTAGE_POPULATION_CONVERGED(1:NCOOL))

!  Find the evaluation points along each HEALPix ray for each particle
   WRITE(6,*) 'Finding evaluation points...'
   CALL FIND_EVALUATION_POINTS(NPART,NRAYS,PARTICLE)

   WRITE(6,*) 'Setting up initial parameters...'

!  Scale the initial abundances by the metallicity
   DO i=1,NSPEC
      IF (.NOT. ANY((/nh,nh2,nhx/) .eq. i)) INITIAL_ABUNDANCE(i) = INITIAL_ABUNDANCE(i)*METALLICITY
   END DO

!  Specify the initial abundances for each particle
   DO P=1,NPART
      PARTICLE(P)%ABUNDANCE = INITIAL_ABUNDANCE
   END DO

!  Calculate the column densities based on the initial abundances
   CALL CALCULATE_COLUMN_DENSITIES(NPART,NRAYS,NSPEC,PARTICLE)
   

!  Calculate the visual extinction for each particle based on the total column density (N_H)
   DO P=1,NPART
      PARTICLE(P)%AV = PARTICLE(P)%TOTAL_COLUMN*AV_FAC
   END DO

!  Calculate the local FUV flux for each particle based on the hard-coded incident fluxes along each ray
   DO P=1,NPART
      PARTICLE(P)%FUV_FLUX = 0.0D0
      DO J=0,NRAYS-1
         PARTICLE(P)%FUV_FLUX = PARTICLE(P)%FUV_FLUX + PARTICLE(P)%FUV_SURFACE(J)*EXP(-UV_FAC*PARTICLE(P)%AV(J))
      END DO
      IF(PARTICLE(P)%FUV_FLUX.LT.1.0D-99) PARTICLE(P)%FUV_FLUX = 0.0D0 ! Impose a lower cut-off of 1E-99
      IF(PARTICLE(P)%FUV_FLUX.GE.1.0D-10) PARTICLE(P)%TYPE = 2
      IF(PARTICLE(P)%FUV_FLUX.LT.1.0D-10) PARTICLE(P)%TYPE = 3
   END DO

!  Calculate the local X-ray flux for each particle
   WRITE(6,*) 'Calculating X-ray properties...'
   CALL CALCULATE_XRAY_PROPERTIES(NPART,NRAYS,PARTICLE)

!  Calculate the dust temperature for each particle
#ifdef CALC_TDUST
   WRITE(6,*) 'Calculating dust temperatures...'
   CALL CALCULATE_DUST_TEMPERATURES(NPART,NRAYS,PARTICLE)
#endif

!  Calculate the initial guess temperature for each particle
   CALL CALCULATE_GUESS_TEMPERATURES(NPART,PARTICLE)

   WRITE(6,*)
   FILE_SUFFIX = '.ini'

!  Write the initial cloud properties to the output file
   CALL WRITE_PROPERTIES(FILE_PREFIX,FILE_SUFFIX,NPART,PARTICLE)

!  Write the initial abundances to the output file
   CALL WRITE_ABUNDANCES(FILE_PREFIX,FILE_SUFFIX,NPART,NSPEC,SPECIES,PARTICLE)

!  Write the visual extinction along each ray to the output file
   CALL WRITE_EXTINCTION(FILE_PREFIX,'.out',NPART,NRAYS,PARTICLE)

!  Initialize the thermal balance iteration properties and set
!  the temperature convergence flag to false for each particle
   PARTICLE%BRACKET_EXPANDED = .FALSE.
   PARTICLE%BINARY_CHOP_SEARCH = .FALSE.
   PARTICLE%TEMPERATURE_CONVERGED = .FALSE.
   PARTICLE%PREVIOUS_DIFFERENCE = 0.0D0
   PARTICLE%PREVIOUS_TEMPERATURE = 0.0D0
   PERCENTAGE_TEMPERATURE_CONVERGED = 0

!  Delete any previous thermal balance iteration output files before starting
   INQUIRE(FILE=TRIM(ADJUSTL(FILE_PREFIX))//'.emis.0001',EXIST=FILE_EXISTS)
   IF(FILE_EXISTS) CALL SYSTEM('rm '//TRIM(ADJUSTL(FILE_PREFIX))//'.emis.[0-9][0-9][0-9][0-9]')
   INQUIRE(FILE=TRIM(ADJUSTL(FILE_PREFIX))//'.cool.0001',EXIST=FILE_EXISTS)
   IF(FILE_EXISTS) CALL SYSTEM('rm '//TRIM(ADJUSTL(FILE_PREFIX))//'.cool.[0-9][0-9][0-9][0-9]')
   INQUIRE(FILE=TRIM(ADJUSTL(FILE_PREFIX))//'.heat.0001',EXIST=FILE_EXISTS)
   IF(FILE_EXISTS) CALL SYSTEM('rm '//TRIM(ADJUSTL(FILE_PREFIX))//'.heat.[0-9][0-9][0-9][0-9]')
   INQUIRE(FILE=TRIM(ADJUSTL(FILE_PREFIX))//'.temp.0001',EXIST=FILE_EXISTS)
   IF(FILE_EXISTS) CALL SYSTEM('rm '//TRIM(ADJUSTL(FILE_PREFIX))//'.temp.[0-9][0-9][0-9][0-9]')

   DO TEMPERATURE_ITERATION=1,TEMPERATURE_ITERATIONS ! Start of thermal balance iteration loop

         WRITE(6,*)
         WRITE(6,*) 'Temperature iteration #',TRIM(NUM2STR(TEMPERATURE_ITERATION))

      !  Reset the chemistry convergence flag for each particle
         PARTICLE%CHEMISTRY_CONVERGED = .FALSE.
         PERCENTAGE_CHEMISTRY_CONVERGED = 0

      !  Delete any previous chemistry iteration output files before starting
         INQUIRE(FILE=TRIM(ADJUSTL(FILE_PREFIX))//'.abun.0001',EXIST=FILE_EXISTS)
         IF(FILE_EXISTS) CALL SYSTEM('rm '//TRIM(ADJUSTL(FILE_PREFIX))//'.abun.[0-9][0-9][0-9][0-9]')

#ifdef OPENMP
   CPU_CHEM = OMP_GET_WTIME()
#else
   CALL CPU_TIME(CPU_CHEM)
#endif

      DO CHEMISTRY_ITERATION=1,CHEMISTRY_ITERATIONS ! Start of chemistry iteration loop

         WRITE(6,*)
         WRITE(6,*) 'Chemistry iteration #',TRIM(NUM2STR(CHEMISTRY_ITERATION))
         WRITE(6,*) 'Calculating reaction rates...'

   !$OMP PARALLEL DEFAULT(NONE) &
   !$OMP    SHARED(NPART,NRAYS,NSPEC,NREAC,NXSEC,CROSS_SECTION,PARTICLE) &
   !$OMP    SHARED(REACTANT,PRODUCT,ALPHA,BETA,GAMMA,RTMIN,RTMAX,DUPLICATE)
   !$OMP DO
   !     Calculate the reaction rate coefficients for each particle
         DO P=1,NPART
            CALL CALCULATE_REACTION_RATES(NRAYS,NSPEC,NREAC,NXSEC,CROSS_SECTION,PARTICLE(P)%GAS_TEMPERATURE,PARTICLE(P)%DUST_TEMPERATURE, &
                                        & PARTICLE(P)%FUV_FLUX,PARTICLE(P)%XRAY_ENERGY_DEPOSITION_RATE,PARTICLE(P)%FUV_SURFACE,PARTICLE(P)%AV, &
                                        & PARTICLE(P)%COLUMN_DENSITY,REACTANT,PRODUCT,ALPHA,BETA,GAMMA,RTMIN,RTMAX,DUPLICATE, &
                                        & PARTICLE(P)%RATE)
            CALL CALCULATE_XRAY_IONIZATION_RATES(PARTICLE(P)%XRAY_ENERGIES,PARTICLE(P)%XRAY_FLUXES,PARTICLE(P)%RATE,REACTANT,NREAC)
         END DO
!$OMP END DO
!$OMP END PARALLEL

!     Determine the list of particles whose abundances need to be updated during this iteration
!     and load their relevant properties into dummy arrays to be passed to the chemistry routine
      CALL UPDATE_ABUNDANCES(TEMPERATURE_ITERATION,CHEMISTRY_ITERATION,NPART,PARTICLE,INITIAL_ABUNDANCE, &
                           & NCHEM,DUMMY_INDEX,DUMMY_ABUNDANCE,DUMMY_RATE,DUMMY_DENSITY,DUMMY_TEMPERATURE)

!     Calculate the species abundances for the selected particles
      WRITE(6,*) 'Calculating chemical abundances for ',TRIM(NUM2STR(NCHEM)),' particles...'
      CALL CALCULATE_ABUNDANCES(DUMMY_ABUNDANCE,DUMMY_RATE,DUMMY_DENSITY,DUMMY_TEMPERATURE,NCHEM,NSPEC,NREAC)

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP    SHARED(NCHEM,PARTICLE,DUMMY_INDEX,DUMMY_ABUNDANCE) &
!$OMP    PRIVATE(P)
!$OMP DO
!     Transfer the updated abundances back into the particle array
      DO I=1,NCHEM
         P = DUMMY_INDEX(I)
         PARTICLE(P)%ABUNDANCE = DUMMY_ABUNDANCE(:,I)
      END DO
!$OMP END DO
!$OMP END PARALLEL

!     Calculate the column densities based on the updated abundances
      WRITE(6,*) 'Calculating column densities...'
      CALL CALCULATE_COLUMN_DENSITIES(NPART,NRAYS,NSPEC,PARTICLE)

!     Write the abundances from the current chemistry iteration to the output file
      WRITE(FILE_SUFFIX,"(I6)") CHEMISTRY_ITERATION
      FILE_SUFFIX = '.'//REPEAT('0',4-LEN_TRIM(ADJUSTL(FILE_SUFFIX)))//TRIM(ADJUSTL(FILE_SUFFIX))
      CALL WRITE_ABUNDANCES(FILE_PREFIX,FILE_SUFFIX,NPART,NSPEC,SPECIES,PARTICLE)

!     Check for convergence of the abundances for each particle
      CALL CHECK_CHEMISTRY_CONVERGENCE(NPART,NSPEC,PARTICLE,PERCENTAGE_CHEMISTRY_CONVERGED)

      WRITE(6,"(' Abundances converged for ',I3,'% of the particles')") INT(PERCENTAGE_CHEMISTRY_CONVERGED)

!     Stop iterations once the abundances have converged for all particles
      IF(INT(PERCENTAGE_CHEMISTRY_CONVERGED).GE.100) EXIT

   END DO ! End of chemistry iteration loop

#ifdef OPENMP
   CPU_CHEM = OMP_GET_WTIME()-CPU_CHEM
#else
   CALL CPU_TIME(CPU_END)
   CPU_CHEM = CPU_END-CPU_CHEM
#endif

!  Update the coolant species abundances and line widths
   CALL UPDATE_COOLANT_ABUNDANCES(NPART,NCOOL,COOLANT,PARTICLE)
   CALL UPDATE_COOLANT_LINEWIDTHS(NPART,NCOOL,COOLANT,PARTICLE)

   WRITE(6,*)
   WRITE(6,*) 'Calculating population densities at LTE...'

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP    SHARED(NPART,NCOOL,COOLANT,PARTICLE) &
!$OMP    PRIVATE(N)
!$OMP DO
!  Calculate the level populations at LTE for all coolant species
   DO P=1,NPART
      DO N=1,NCOOL
         CALL CALCULATE_LTE_POPULATIONS(COOLANT(N)%NLEVEL,COOLANT(N)%ENERGY,COOLANT(N)%WEIGHT, &
                                      & PARTICLE(P)%COOLANT(N)%DENSITY,PARTICLE(P)%GAS_TEMPERATURE, &
                                      & PARTICLE(P)%COOLANT(N)%POPULATION)
      END DO
   END DO
!$OMP END DO
!$OMP END PARALLEL

!  Calculate the line opacities for each coolant transition
   WRITE(6,*) 'Calculating emission line opacities...'
   CALL CALCULATE_LINE_OPACITIES(NPART,NCOOL,NRAYS,COOLANT,PARTICLE)
#ifdef USE_ALI
   CALL CALCULATE_LAMBDA_OPERATOR(NPART,NCOOL,NRAYS,COOLANT,PARTICLE)
#endif

!  Write the initial population densities to the output file
   FILE_SUFFIX = '.ini'
   CALL WRITE_POPULATIONS(FILE_PREFIX,FILE_SUFFIX,NPART,NCOOL,COOLANT,PARTICLE)

!  Reset the population convergence flag for each coolant of each particle
   DO P=1,NPART
      PARTICLE(P)%COOLANT%CONVERGED = .FALSE.
   END DO
   COOLANT%CONVERGED = .FALSE.
   PERCENTAGE_POPULATION_CONVERGED = 0

!  Delete any previous population density iteration output files before starting
   INQUIRE(FILE=TRIM(ADJUSTL(FILE_PREFIX))//'.pop.0001',EXIST=FILE_EXISTS)
   IF(FILE_EXISTS) CALL SYSTEM('rm '//TRIM(ADJUSTL(FILE_PREFIX))//'.pop.[0-9][0-9][0-9][0-9]')

#ifdef OPENMP
   CPU_POP = OMP_GET_WTIME()
#else
   CALL CPU_TIME(CPU_POP)
#endif

   DO POPULATION_ITERATION=1,POPULATION_ITERATIONS ! Start of population density iteration loop

      WRITE(6,*)
      WRITE(6,*) 'Population iteration #',TRIM(NUM2STR(POPULATION_ITERATION))
      WRITE(6,*) 'Temperature iteration #',TRIM(NUM2STR(TEMPERATURE_ITERATION))
      WRITE(6,*) 'Calculating population densities...'

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP    SHARED(NPART,NCOOL,NRAYS,COOLANT,PARTICLE) &
!$OMP    SHARED(POPULATION_ITERATION,POPULATION_ITERATION_LIMIT) &
!$OMP    SHARED(PERCENTAGE_POPULATION_CONVERGED) &
!$OMP    SHARED(POPULATION_PERCENTAGE_LIMIT) &
!$OMP    PRIVATE(N)
!$OMP DO
!     Calculate the level populations for all coolant species using the LVG approximation
      DO P=1,NPART
         DO N=1,NCOOL

            IF(COOLANT(N)%CONVERGED) CYCLE

!           If convergence is proving difficult for a coolant species, attempt to speed it
!           up by fixing the population densities of particles that have already converged
            IF(POPULATION_ITERATION.GT.2*POPULATION_ITERATION_LIMIT .AND. &
             & PERCENTAGE_POPULATION_CONVERGED(N).GT.POPULATION_PERCENTAGE_LIMIT .AND. &
             & PARTICLE(P)%COOLANT(N)%CONVERGED) CYCLE

!!$!           If convergence is proving particularly difficult for CO, then attempt to speed it
!!$!           up by fixing the population densities of particles that have already converged
!!$            IF(POPULATION_ITERATION.GT.4*POPULATION_ITERATION_LIMIT .AND. &
!!$             & TRIM(ADJUSTL(COOLANT(N)%NAME)).EQ."CO" .AND. &
!!$             & PARTICLE(P)%COOLANT(N)%CONVERGED) CYCLE

            CALL CALCULATE_LEVEL_POPULATIONS(NRAYS,N,COOLANT(N),PARTICLE(P))

!           If convergence is proving difficult for a coolant species, attempt to reduce
!           oscillations in its population densities by averaging them between iterations
            IF(POPULATION_ITERATION.GT.POPULATION_ITERATION_LIMIT) THEN
               PARTICLE(P)%COOLANT(N)%POPULATION = (PARTICLE(P)%COOLANT(N)%PREVIOUS_POPULATION &
                                                & + PARTICLE(P)%COOLANT(N)%POPULATION)/2.0D0
            END IF

!!$!           If convergence is proving difficult for a coolant species, attempt to speed it
!!$!           up by recalculating the level populations without updating the line opacities
!!$            IF(POPULATION_ITERATION.GT.2*POPULATION_ITERATION_LIMIT .AND. &
!!$             & MODULO(POPULATION_ITERATION,2).EQ.0) THEN
!!$               CALL CALCULATE_LEVEL_POPULATIONS(NRAYS,N,COOLANT(N),PARTICLE(P))
!!$            END IF

!!$!           If convergence is proving particularly difficult for CO, then attempt to speed it
!!$!           up by by recalculating the level populations without updating the line opacities
!!$            IF(POPULATION_ITERATION.GT.2*POPULATION_ITERATION_LIMIT .AND. &
!!$             & TRIM(ADJUSTL(COOLANT(N)%NAME)).EQ."CO" .AND. &
!!$             & MODULO(POPULATION_ITERATION,2).EQ.0) THEN
!!$               CALL CALCULATE_LEVEL_POPULATIONS(NRAYS,N,COOLANT(N),PARTICLE(P))
!!$            END IF

         END DO
      END DO
!$OMP END DO
!$OMP END PARALLEL

!     Calculate the line opacities based on the updated population densities
      WRITE(6,*) 'Calculating emission line opacities...'
      CALL CALCULATE_LINE_OPACITIES(NPART,NCOOL,NRAYS,COOLANT,PARTICLE)
#ifdef USE_ALI
      CALL CALCULATE_LAMBDA_OPERATOR(NPART,NCOOL,NRAYS,COOLANT,PARTICLE)
#endif

!!$!     Write the population densities from the current iteration to the output file
!!$      WRITE(FILE_SUFFIX,"(I6)") POPULATION_ITERATION
!!$      FILE_SUFFIX = '.'//REPEAT('0',4-LEN_TRIM(ADJUSTL(FILE_SUFFIX)))//TRIM(ADJUSTL(FILE_SUFFIX))
!!$      CALL WRITE_POPULATIONS(FILE_PREFIX,FILE_SUFFIX,NPART,NCOOL,COOLANT,PARTICLE)

!     Check for convergence of the population densities for each coolant
      CALL CHECK_POPULATION_CONVERGENCE(NPART,NCOOL,COOLANT,PARTICLE,PERCENTAGE_POPULATION_CONVERGED)

      WRITE(6,"(' Convergence status:',20(' [',A,': ',I3,'%]',:))") 'T',INT(PERCENTAGE_TEMPERATURE_CONVERGED),(TRIM(ADJUSTL(COOLANT(N)%NAME)),INT(PERCENTAGE_POPULATION_CONVERGED(N)),N=1,NCOOL)

!     Stop iterations once all coolants are converged for all particles
      IF(ALL(INT(PERCENTAGE_POPULATION_CONVERGED).GE.100)) EXIT

   END DO ! End of population density iteration loop

#ifdef OPENMP
   CPU_POP = OMP_GET_WTIME()-CPU_POP
#else
   CALL CPU_TIME(CPU_END)
   CPU_POP = CPU_END-CPU_POP
#endif

!!$!  If some particles remain unconverged after the maximum number of iterations has been reached,
!!$!  calculate approximate values for their population densities using nearest neighbour weighted
!!$!  (inverse distance weighted) interpolation between the populations of neighbouring particles
!!$   IF(ANY(INT(PERCENTAGE_POPULATION_CONVERGED).LT.100)) THEN
!!$      WRITE(6,*)
!!$      WRITE(6,*) 'Using interpolated population densities for the unconverged particles...'
!!$      PERCENTAGE_POPULATION_CONVERGED = 0
!!$      DO P=1,NPART
!!$         DO N=1,NCOOL
!!$            IF(.NOT.PARTICLE(P)%COOLANT(N)%CONVERGED) THEN
!!$               PERCENTAGE_POPULATION_CONVERGED(N) = PERCENTAGE_POPULATION_CONVERGED(N) + 1
!!$               CALL USE_NEARBY_POPULATIONS(NPART,P,N,PARTICLE)
!!$            END IF
!!$         END DO
!!$      END DO
!!$      WRITE(6,"(' Particles affected:',10(' [',A,': ',I4,']',:))") (TRIM(ADJUSTL(COOLANT(N)%NAME)),INT(PERCENTAGE_POPULATION_CONVERGED(N)),N=1,NCOOL)
!!$      WRITE(FILE_SUFFIX,"(I6)") POPULATION_ITERATION
!!$      FILE_SUFFIX = '.'//REPEAT('0',4-LEN_TRIM(ADJUSTL(FILE_SUFFIX)))//TRIM(ADJUSTL(FILE_SUFFIX))
!!$      CALL WRITE_POPULATIONS(FILE_PREFIX,FILE_SUFFIX,NPART,NCOOL,COOLANT,PARTICLE)
!!$   END IF

   WRITE(6,*)
   WRITE(6,*) 'Calculating cooling rates...'

! $OMP PARALLEL DEFAULT(NONE) &
! $OMP SHARED(NPART,NCOOL,nH,nelect,PARTICLE) &
! $OMP PRIVATE(N)
! $OMP DO
!  Calculate the cooling rate due to the Lyman-alpha emission for each particle
!  using the analytical expression of Spitzer (1978) neglecting photon trapping
   DO P=1,NPART
      DO N=1,NCOOL
         IF(PARTICLE(P)%COOLANT(N)%NAME.EQ."H") THEN
            PARTICLE(P)%COOLANT(N)%EMISSIVITY(2,1) = 7.3D-19*(PARTICLE(P)%ABUNDANCE(nelect)*PARTICLE(P)%GAS_DENSITY) &
                                                          & *(PARTICLE(P)%ABUNDANCE(nH)*PARTICLE(P)%GAS_DENSITY) &
                                                          & *EXP(-118400.0D0/PARTICLE(P)%GAS_TEMPERATURE)
         END IF
      END DO
   END DO
! $OMP END DO
! $OMP END PARALLEL

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(NPART,NCOOL,PARTICLE) &
!$OMP PRIVATE(N)
!$OMP DO
!  Calculate the cooling rates for each particle
   DO P=1,NPART
      DO N=1,NCOOL
         PARTICLE(P)%COOLING_RATE(N) = SUM(PARTICLE(P)%COOLANT(N)%EMISSIVITY,MASK=.NOT.ISNAN(PARTICLE(P)%COOLANT(N)%EMISSIVITY))
      END DO
   END DO
   ! !  Calculate the cooling rates for each particle
   ! DO P=1,NPART
   !    PARTICLE(P)%COOLING_RATE(1)=1.0d-16
   !    DO N=2,NCOOL
   !       PARTICLE(P)%COOLING_RATE(N) = 0.0!SUM(PARTICLE(P)%COOLANT(N)%EMISSIVITY,MASK=.NOT.ISNAN(PARTICLE(P)%COOLANT(N)%EMISSIVITY))
   !    END DO
   ! END DO
!$OMP END DO
!$OMP END PARALLEL

!  Check for particles with cooling rates that deviate significantly from those of
!  their neighbours and assign them new values using the nearest-neighbour average.
   WRITE(6,*)
   WRITE(6,*) 'Using interpolated cooling rates for the unconverged particles...'
   CALL ADJUST_BAD_COOLING_RATES(NPART,NCOOL,NRAYS,COOLANT,PARTICLE)
   WRITE(6,*)

   WRITE(6,*) 'Calculating heating rates...'

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(NPART,NSPEC,NREAC,PARTICLE)
!$OMP DO
!  Calculate the heating rates for each particle
   DO P=1,NPART
      CALL CALCULATE_HEATING_RATES(NSPEC,NREAC,PARTICLE(P))
   END DO
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(NPART,PARTICLE)
!$OMP DO
!  Check if thermal balance has been reached and update the temperatures if not
   DO P=1,NPART
      CALL UPDATE_TEMPERATURE(PARTICLE(P))
   END DO
!$OMP END DO
!$OMP END PARALLEL

   WRITE(6,*)
   WRITE(FILE_SUFFIX,"(I6)") TEMPERATURE_ITERATION
   FILE_SUFFIX = '.'//REPEAT('0',4-LEN_TRIM(ADJUSTL(FILE_SUFFIX)))//TRIM(ADJUSTL(FILE_SUFFIX))

!  Write the line emissivities from the current iteration to the output file
   CALL WRITE_EMISSIVITIES(FILE_PREFIX,FILE_SUFFIX,NPART,NCOOL,COOLANT,PARTICLE)

!  Write the cooling rates from the current iteration to the output file
   CALL WRITE_COOLING_RATES(FILE_PREFIX,FILE_SUFFIX,NPART,NCOOL,COOLANT,PARTICLE)

!  Write the heating rates from the current iteration to the output file
   CALL WRITE_HEATING_RATES(FILE_PREFIX,FILE_SUFFIX,NPART,NHEAT,PARTICLE)

!  Write the thermal balance properties from the current iteration to the output file
   CALL WRITE_THERMAL_BALANCE(FILE_PREFIX,FILE_SUFFIX,NPART,PARTICLE)

!  Write the convergence status flags from the current iteration to the output file
   CALL WRITE_CONVERGENCE_STATUS(FILE_PREFIX,FILE_SUFFIX,NPART,NCOOL,COOLANT,PARTICLE)

!  Stop iterations once temperatures have converged for all particles
   IF(INT(PERCENTAGE_TEMPERATURE_CONVERGED).GE.100) EXIT

!  Update the percentage of particles whose temperatures have converged
   PERCENTAGE_TEMPERATURE_CONVERGED = COUNT(PARTICLE%TEMPERATURE_CONVERGED)/REAL(NPART,KIND=DP)*100

   WRITE(6,*)
   WRITE(6,"(' Temperature converged for ',I3,'% of the particles')") INT(PERCENTAGE_TEMPERATURE_CONVERGED)

   END DO ! End of thermal balance iteration loop

   WRITE(6,*)
   WRITE(6,"(' Temperature converged for ',I3,'% of the particles')") INT(PERCENTAGE_TEMPERATURE_CONVERGED)

!!$!  If some particles still fail to meet the thermal balance criterion after global convergence has been
!!$!  declared (or after the maximum number of iterations has been reached), assign them temperatures that
!!$!  are more consistent with those of surrounding particles by using nearest neighbour weighted (inverse
!!$!  distance weighted) interpolation between nearby particles whose temperatures have properly converged
!!$   IF(PERCENTAGE_TEMPERATURE_CONVERGED.GE.100 .OR.TEMPERATURE_ITERATION.EQ.TEMPERATURE_ITERATIONS) THEN
!!$      CALL USE_NEARBY_TEMPERATURES(NPART,PARTICLE)
!!$   END IF

   WRITE(6,*)
   WRITE(6,*) 'Writing output files...'
   WRITE(6,*)

   FILE_SUFFIX = '.out'

!  Write the final cloud properties to the output file
   CALL WRITE_PROPERTIES(FILE_PREFIX,FILE_SUFFIX,NPART,PARTICLE)

!  Write the chemical network analysis results to the output file
   CALL ANALYSE_CHEMISTRY(FILE_PREFIX,FILE_SUFFIX,NPART,NSPEC,NREAC, &
                          END_TIME,SPECIES,REACTANT,PRODUCT,PARTICLE, &
                        & (/0.0D0,1.0D0,1.0D1,1.0D2,1.0D3,3.0D3,1.0D4/), &
                        & (/"H  ","H2 ","C  ","CO ","H+ ","O+ ","e- "/))

!  Write the final reaction rates to the output file
   CALL WRITE_REACTION_RATES(FILE_PREFIX,FILE_SUFFIX,NPART,NREAC,REACTANT,PARTICLE)

!  Write the final abundances to the output file
   CALL WRITE_ABUNDANCES(FILE_PREFIX,FILE_SUFFIX,NPART,NSPEC,SPECIES,PARTICLE)

!  Write the final population densities to the output file
   CALL WRITE_POPULATIONS(FILE_PREFIX,FILE_SUFFIX,NPART,NCOOL,COOLANT,PARTICLE)

!  Write the final line opacities to the output file
   DO J=0,NRAYS-1
      CALL WRITE_OPACITIES(FILE_PREFIX,FILE_SUFFIX,NPART,NCOOL,NRAYS,J,COOLANT,PARTICLE)
   END DO

!  Write the final line emissivities to the output file
   CALL WRITE_EMISSIVITIES(FILE_PREFIX,FILE_SUFFIX,NPART,NCOOL,COOLANT,PARTICLE)

!  Write the final cooling rates to the output file
   CALL WRITE_COOLING_RATES(FILE_PREFIX,FILE_SUFFIX,NPART,NCOOL,COOLANT,PARTICLE)

!  Write the final heating rates to the output file
   CALL WRITE_HEATING_RATES(FILE_PREFIX,FILE_SUFFIX,NPART,NHEAT,PARTICLE)

!  Write the final thermal balance properties to the output file
   CALL WRITE_THERMAL_BALANCE(FILE_PREFIX,FILE_SUFFIX,NPART,PARTICLE)

!  Write the final convergence status flags to the output file
   CALL WRITE_CONVERGENCE_STATUS(FILE_PREFIX,FILE_SUFFIX,NPART,NCOOL,COOLANT,PARTICLE)

#ifdef OPENMP
   CPU_END = OMP_GET_WTIME()
#else
   CALL CPU_TIME(CPU_END)
#endif

   WRITE(6,*)
   WRITE(6,*) 'Simulation time = ',TRIM(NUM2STR(CPU_END-CPU_START)),' seconds.'
   WRITE(6,*) 'Abundance iteration time  = ',TRIM(NUM2STR(CPU_CHEM)),' seconds.'
   WRITE(6,*) 'Population iteration time = ',TRIM(NUM2STR(CPU_POP)),' seconds.'
   WRITE(6,*) 'Finished!'
   WRITE(6,*)

!  Close the log files
   CALL CLOSE_LOGFILES

!  Deallocate the arrays/pointers
   DEALLOCATE(SPECIES,INITIAL_ABUNDANCE,MOLECULAR_MASS)
   DEALLOCATE(REACTANT,PRODUCT,ALPHA,BETA,GAMMA,RTMIN,RTMAX,DUPLICATE)
   DEALLOCATE(DUMMY_INDEX,DUMMY_ABUNDANCE,DUMMY_RATE,DUMMY_DENSITY,DUMMY_TEMPERATURE)
   DEALLOCATE(PERCENTAGE_POPULATION_CONVERGED)
   DEALLOCATE(CROSS_SECTION,COOLANT,PARTICLE)

END PROGRAM UCL_PDR
!=======================================================================
