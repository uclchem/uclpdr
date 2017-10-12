!=======================================================================
!
!  Print the splash screen text and author/contributor list to screen.
!
!-----------------------------------------------------------------------
SUBROUTINE PRINT_SPLASHSCREEN

  IMPLICIT NONE

   WRITE(6,*)
   WRITE(6,*) '=============================================================================='
   WRITE(6,*) ' **       **      *****     **             *******    ******      *******     '
   WRITE(6,*) ' **       **    ***   ***   **             **    ***  **   ***    **    ***   '
   WRITE(6,*) ' **       **   **       **  **             **      ** **     **   **      **  '
   WRITE(6,*) ' **       **  **            **   *******   **    ***  **      **  **    ***   '
   WRITE(6,*) ' **       **  **            **   *******   *******    **      **  ********    '
   WRITE(6,*) ' **       **  **            **             **         **      **  **    ***   '
   WRITE(6,*) ' **       **   **       **  **             **         **     **   **     **   '
   WRITE(6,*) '  ***   ***     ***   ***   *********      **         **   ***    **      **  '
   WRITE(6,*) '    *****         *****      ********      **         ******      **       ** '
   WRITE(6,*) '=============================================================================='
   WRITE(6,*) '**********************   Coders: T.A.Bell, T.G.Bisbas   **********************'
   WRITE(6,*) '***************   Collaborators: S.Viti, J.Yates, M.Barlow   *****************'
   WRITE(6,*) '*************************        Version 2.0        **************************'
   WRITE(6,*) '=============================================================================='
   WRITE(6,*)

END SUBROUTINE PRINT_SPLASHSCREEN
!=======================================================================

!=======================================================================
!
!  Write the physical properties for each particle to the specified
!  output file.
!
!-----------------------------------------------------------------------
SUBROUTINE WRITE_PROPERTIES(PREFIX,SUFFIX,NPART,PARTICLE)

   USE HEALPIX_TYPES
   USE PARTICLE_MODULE

   IMPLICIT NONE

   CHARACTER(LEN=*),    INTENT(IN) :: PREFIX,SUFFIX
   INTEGER(KIND=I4B),   INTENT(IN) :: NPART
   TYPE(PARTICLE_TYPE), INTENT(IN) :: PARTICLE(:)

   CHARACTER(LEN=256) :: FILENAME
   INTEGER(KIND=I4B)  :: I,P

!  Specify the output file name
   FILENAME=TRIM(ADJUSTL(PREFIX))//TRIM(ADJUSTL('.prop'))//TRIM(ADJUSTL(SUFFIX))
   WRITE(6,*) 'Writing file [',TRIM(FILENAME),']'

!  Open and write to the output file
   OPEN(UNIT=1,FILE=FILENAME,STATUS='REPLACE')
   WRITE(1,"('Particle #  n_H (cm^-3)  n_g (cm^-3)   T_gas (K)   T_dust (K)   chi(Draine)  F_X (X-ray)')")
   DO P=1,NPART
      WRITE(1,"(I10,6ES13.5)") P,PARTICLE(P)%GAS_DENSITY,PARTICLE(P)%DUST_DENSITY, &
                               & PARTICLE(P)%GAS_TEMPERATURE,PARTICLE(P)%DUST_TEMPERATURE, &
                               & PARTICLE(P)%FUV_FLUX,PARTICLE(P)%XRAY_FLUX
   END DO
   CLOSE(1)

   RETURN
END SUBROUTINE WRITE_PROPERTIES
!=======================================================================

!=======================================================================
!
!  Write the visual extinction along each ray for each particle to the
!  specified output file.
!
!-----------------------------------------------------------------------
SUBROUTINE WRITE_EXTINCTION(PREFIX,SUFFIX,NPART,NRAYS,PARTICLE)

   USE HEALPIX_TYPES
   USE PARTICLE_MODULE

   IMPLICIT NONE

   CHARACTER(LEN=*),    INTENT(IN) :: PREFIX,SUFFIX
   INTEGER(KIND=I4B),   INTENT(IN) :: NPART,NRAYS
   TYPE(PARTICLE_TYPE), INTENT(IN) :: PARTICLE(:)

   CHARACTER(LEN=256) :: FILENAME
   CHARACTER(LEN=11), ALLOCATABLE :: AV_LABELS(:)

   INTEGER(KIND=I4B)  :: I,P

!  Create the extinction label for each ray
   ALLOCATE(AV_LABELS(1:NRAYS))
   DO I=1,NRAYS
      WRITE(AV_LABELS(I),"(I3)") I
      AV_LABELS(I)='A_V('//TRIM(ADJUSTL(AV_LABELS(I)))//') mag'
      AV_LABELS(I)=REPEAT(' ',(LEN(AV_LABELS(I))-LEN_TRIM(AV_LABELS(I)))/2) &
                 & //TRIM(ADJUSTL(AV_LABELS(I)))
   END DO

!  Specify the output file name
   FILENAME=TRIM(ADJUSTL(PREFIX))//TRIM(ADJUSTL('.av'))//TRIM(ADJUSTL(SUFFIX))
   WRITE(6,*) 'Writing file [',TRIM(FILENAME),']'

!  Open and write to the output file
   OPEN(UNIT=1,FILE=FILENAME,STATUS='REPLACE')
   WRITE(1,"('Particle #    x (cm)       y (cm)       z (cm)   ',500(2X,A11))") AV_LABELS
   DO P=1,NPART
      WRITE(1,"(I10,SP,3ES13.5,SS,500ES13.5)") P,(PARTICLE(P)%COORDINATES(I),I=1,3),(PARTICLE(P)%AV(I),I=0,NRAYS-1)
   END DO
   CLOSE(1)

   RETURN
END SUBROUTINE WRITE_EXTINCTION
!=======================================================================

!=======================================================================
!
!  Write the reaction rate coefficients for each particle to the
!  specified output file.
!
!-----------------------------------------------------------------------
SUBROUTINE WRITE_REACTION_RATES(PREFIX,SUFFIX,NPART,NREAC,REACTANT,PARTICLE)

   USE HEALPIX_TYPES
   USE PARTICLE_MODULE

   IMPLICIT NONE

   CHARACTER(LEN=*),    INTENT(IN) :: PREFIX,SUFFIX
   INTEGER(KIND=I4B),   INTENT(IN) :: NPART,NREAC
   CHARACTER(LEN=*),    INTENT(IN) :: REACTANT(:,:)
   TYPE(PARTICLE_TYPE), INTENT(IN) :: PARTICLE(:)

   CHARACTER(LEN=256) :: FILENAME
   CHARACTER(LEN=11), ALLOCATABLE :: REACTION_LABELS(:)

   INTEGER(KIND=I4B)  :: I,P

!  Create the reaction label for each rate
   ALLOCATE(REACTION_LABELS(1:NREAC))
   DO I=1,NREAC
      REACTION_LABELS(I)=TRIM(ADJUSTL(REACTANT(I,1)))//','//TRIM(ADJUSTL(REACTANT(I,2)))
      REACTION_LABELS(I)=REPEAT(' ',(LEN(REACTION_LABELS(I))-LEN_TRIM(REACTION_LABELS(I)))/2) &
                       & //TRIM(ADJUSTL(REACTION_LABELS(I)))
   END DO

!  Specify the output file name
   FILENAME=TRIM(ADJUSTL(PREFIX))//TRIM(ADJUSTL('.reac'))//TRIM(ADJUSTL(SUFFIX))
   WRITE(6,*) 'Writing file [',TRIM(FILENAME),']'

!  Open and write to the output file
   OPEN(UNIT=1,FILE=FILENAME,STATUS='REPLACE')
   WRITE(1,"('Particle #',9999(2X,A11))") REACTION_LABELS
   DO P=1,NPART
      WRITE(1,"(I10,9999ES13.5)") P,(PARTICLE(P)%RATE(I),I=1,NREAC)
   END DO
   CLOSE(1)

   DEALLOCATE(REACTION_LABELS)

   RETURN
END SUBROUTINE WRITE_REACTION_RATES
!=======================================================================

!=======================================================================
!
!  Write the fractional abundances for each particle to the specified
!  output file.
!
!-----------------------------------------------------------------------
SUBROUTINE WRITE_ABUNDANCES(PREFIX,SUFFIX,NPART,NSPEC,SPECIES,PARTICLE)

   USE HEALPIX_TYPES
   USE PARTICLE_MODULE

   IMPLICIT NONE

   CHARACTER(LEN=*),    INTENT(IN) :: PREFIX,SUFFIX
   INTEGER(KIND=I4B),   INTENT(IN) :: NPART,NSPEC
   CHARACTER(LEN=*),    INTENT(IN) :: SPECIES(:)
   TYPE(PARTICLE_TYPE), INTENT(IN) :: PARTICLE(:)

   CHARACTER(LEN=256) :: FILENAME
   CHARACTER(LEN=11), ALLOCATABLE :: ABUNDANCE_LABELS(:)

   INTEGER(KIND=I4B)  :: I,P

!  Create the abundance label for each species
   ALLOCATE(ABUNDANCE_LABELS(1:NSPEC))
   DO I=1,NSPEC
      ABUNDANCE_LABELS(I)=REPEAT(' ',(LEN(ABUNDANCE_LABELS(I))-LEN_TRIM(SPECIES(I))-3)/2) &
                        & //'X('//TRIM(ADJUSTL(SPECIES(I)))//')'
   END DO

!  Specify the output file name
   FILENAME=TRIM(ADJUSTL(PREFIX))//TRIM(ADJUSTL('.abun'))//TRIM(ADJUSTL(SUFFIX))
   WRITE(6,*) 'Writing file [',TRIM(FILENAME),']'

!  Open and write to the output file
   OPEN(UNIT=1,FILE=FILENAME,STATUS='REPLACE')
   WRITE(1,"('Particle #',999(2X,A11))") ABUNDANCE_LABELS
   DO P=1,NPART
      WRITE(1,"(I10,999ES13.5)") P,(PARTICLE(P)%ABUNDANCE(I),I=1,NSPEC)
   END DO
   CLOSE(1)

   DEALLOCATE(ABUNDANCE_LABELS)

   RETURN
END SUBROUTINE WRITE_ABUNDANCES
!=======================================================================

!=======================================================================
!
!  Write the population densities for each coolant and for each particle
!  to the specified output file.
!
!-----------------------------------------------------------------------
SUBROUTINE WRITE_POPULATIONS(PREFIX,SUFFIX,NPART,NCOOL,COOLANT,PARTICLE)

   USE HEALPIX_TYPES
   USE COOLANT_MODULE
   USE PARTICLE_MODULE

   IMPLICIT NONE

   CHARACTER(LEN=*),    INTENT(IN) :: PREFIX,SUFFIX
   INTEGER(KIND=I4B),   INTENT(IN) :: NPART,NCOOL
   TYPE(COOLANT_TYPE),  INTENT(IN) :: COOLANT(:)
   TYPE(PARTICLE_TYPE), INTENT(IN) :: PARTICLE(:)

   CHARACTER(LEN=256) :: FILENAME
   CHARACTER(LEN=11), ALLOCATABLE :: POPULATION_LABELS(:)

   INTEGER(KIND=I4B)  :: I,N,P

!  Create the population labels for each coolant
   ALLOCATE(POPULATION_LABELS(1:SUM((/(COOLANT(N)%NLEVEL,N=1,NCOOL)/))))
   P=1
   DO N=1,NCOOL
      DO I=1,COOLANT(N)%NLEVEL
         WRITE(POPULATION_LABELS(P),"(I3)") I-1
         POPULATION_LABELS(P)='n('//TRIM(ADJUSTL(COOLANT(N)%NAME))//','//TRIM(ADJUSTL(POPULATION_LABELS(P)))//')'
         POPULATION_LABELS(P)=REPEAT(' ',(LEN(POPULATION_LABELS(P))-LEN_TRIM(POPULATION_LABELS(P)))/2) &
                            & //TRIM(ADJUSTL(POPULATION_LABELS(P)))
         P=P+1
      END DO
   END DO

!  Specify the output file name
   FILENAME=TRIM(ADJUSTL(PREFIX))//TRIM(ADJUSTL('.pop'))//TRIM(ADJUSTL(SUFFIX))
   WRITE(6,*) 'Writing file [',TRIM(FILENAME),']'

!  Open and write to the output file
   OPEN(UNIT=1,FILE=FILENAME,STATUS='REPLACE')
   WRITE(1,"('Particle #',999(2X,A11))") POPULATION_LABELS
   DO P=1,NPART
      WRITE(1,"(I10,999ES13.5)") P,(PARTICLE(P)%COOLANT(N)%POPULATION,N=1,NCOOL)
   END DO
   CLOSE(1)

   DEALLOCATE(POPULATION_LABELS)

   RETURN
END SUBROUTINE WRITE_POPULATIONS
!=======================================================================

!=======================================================================
!
!  Write the optical depths along a given ray for the emission lines of
!  each coolant and for each particle to the specified output file.
!
!-----------------------------------------------------------------------
SUBROUTINE WRITE_OPACITIES(PREFIX,SUFFIX,NPART,NCOOL,NRAYS,RAY,COOLANT,PARTICLE)

   USE HEALPIX_TYPES
   USE COOLANT_MODULE
   USE PARTICLE_MODULE
   USE FUNCTIONS_MODULE

   IMPLICIT NONE

   CHARACTER(LEN=*),    INTENT(IN) :: PREFIX,SUFFIX
   INTEGER(KIND=I4B),   INTENT(IN) :: NPART,NCOOL
   INTEGER(KIND=I4B),   INTENT(IN) :: NRAYS,RAY
   TYPE(COOLANT_TYPE),  INTENT(IN) :: COOLANT(:)
   TYPE(PARTICLE_TYPE), INTENT(IN) :: PARTICLE(:)

   CHARACTER(LEN=256) :: FILENAME
   CHARACTER(LEN=11), ALLOCATABLE :: LINE_LABELS(:)

   INTEGER(KIND=I4B)  :: I,J,N,P

!  Create the transition label for each emission line
   ALLOCATE(LINE_LABELS(1:COUNT((/(COOLANT(N)%A_COEFF,N=1,NCOOL)/).GT.0)))
   P=1
   DO N=1,NCOOL ! Loop over coolants
      DO I=1,COOLANT(N)%NLEVEL ! Loop over levels (i)
         DO J=1,COOLANT(N)%NLEVEL ! Loop over levels (j)
            IF(COOLANT(N)%A_COEFF(I,J).EQ.0) CYCLE

!           Assume all coolants with < 6 levels are atoms with fine-structure lines
            IF(COOLANT(N)%NLEVEL.LT.6) THEN

!              Label fine-structure lines with their wavelength in micron (or in Angstrom, if appropriate)
               IF(COOLANT(N)%FREQUENCY(I,J).GE.3.0D14) THEN
                  WRITE(LINE_LABELS(P),"(I4)") NINT(C/COOLANT(N)%FREQUENCY(I,J)*1.0D8) ! Wavelength in Angstrom (nearest int)
                  LINE_LABELS(P)=TRIM(ADJUSTL(LINE_LABELS(P)))//'A'
               ELSE IF(COOLANT(N)%FREQUENCY(I,J).GE.3.0D13) THEN
                  WRITE(LINE_LABELS(P),"(F5.2)") C/COOLANT(N)%FREQUENCY(I,J)*1.0D4 ! Wavelength in micron (2 dp accuracy)
                  LINE_LABELS(P)=TRIM(ADJUSTL(LINE_LABELS(P)))//'um'
               ELSE IF(COOLANT(N)%FREQUENCY(I,J).GE.1.0D13) THEN
                  WRITE(LINE_LABELS(P),"(F5.1)") C/COOLANT(N)%FREQUENCY(I,J)*1.0D4 ! Wavelength in micron (1 dp accuracy)
                  LINE_LABELS(P)=TRIM(ADJUSTL(LINE_LABELS(P)))//'um'
               ELSE
                  WRITE(LINE_LABELS(P),"(I4)") NINT(C/COOLANT(N)%FREQUENCY(I,J)*1.0D4) ! Wavelength in micron (nearest int)
                  LINE_LABELS(P)=TRIM(ADJUSTL(LINE_LABELS(P)))//'um'
               END IF

!              Handle labelling of both neutral and ionized atomic species
               LINE_LABELS(P)='['//COOLANT(N)%NAME(1:VERIFY(COOLANT(N)%NAME,'+ ',.TRUE.))  &
                               & //REPEAT('I',COUNT_SUBSTRING(COOLANT(N)%NAME,'+'))//'I] ' &
                               & //TRIM(ADJUSTL(LINE_LABELS(P)))

!           Assume all coolants with more levels are molecules with pure rotational or ro-vibrational lines
            ELSE

!              If the line frequency corresponds to a wavelength in the micron range
!              or shorter, then label the line with its wavelength in micron instead
               IF(COOLANT(N)%FREQUENCY(I,J).GE.3.0D14) THEN
                  WRITE(LINE_LABELS(P),"(I4)") NINT(C/COOLANT(N)%FREQUENCY(I,J)*1.0D8) ! Wavelength in Angstrom (nearest int)
                  LINE_LABELS(P)=TRIM(ADJUSTL(COOLANT(N)%NAME))//' '//TRIM(ADJUSTL(LINE_LABELS(P)))//'A'
               ELSE IF(COOLANT(N)%FREQUENCY(I,J).GE.3.0D13) THEN
                  WRITE(LINE_LABELS(P),"(F5.2)") C/COOLANT(N)%FREQUENCY(I,J)*1.0D4 ! Wavelength in micron (2 dp accuracy)
                  LINE_LABELS(P)=TRIM(ADJUSTL(COOLANT(N)%NAME))//' '//TRIM(ADJUSTL(LINE_LABELS(P)))//'um'
               ELSE IF(COOLANT(N)%FREQUENCY(I,J).GE.1.0D13) THEN
                  WRITE(LINE_LABELS(P),"(F5.1)") C/COOLANT(N)%FREQUENCY(I,J)*1.0D4 ! Wavelength in micron (1 dp accuracy)
                  LINE_LABELS(P)=TRIM(ADJUSTL(COOLANT(N)%NAME))//' '//TRIM(ADJUSTL(LINE_LABELS(P)))//'um'
               ELSE IF(COOLANT(N)%FREQUENCY(I,J).GE.6.0D12) THEN
                  WRITE(LINE_LABELS(P),"(I4)") NINT(C/COOLANT(N)%FREQUENCY(I,J)*1.0D4) ! Wavelength in micron (nearest int)
                  LINE_LABELS(P)=TRIM(ADJUSTL(COOLANT(N)%NAME))//' '//TRIM(ADJUSTL(LINE_LABELS(P)))//'um'
               ELSE
                  IF(J.GT.100) THEN
                     WRITE(LINE_LABELS(P),"(I4,'-',I3)") I-1,J-1
                  ELSE IF(J.GT.10) THEN
                     WRITE(LINE_LABELS(P),"(I4,'-',I2)") I-1,J-1
                  ELSE
                     WRITE(LINE_LABELS(P),"(I4,'-',I1)") I-1,J-1
                  END IF
                  LINE_LABELS(P)=TRIM(ADJUSTL(COOLANT(N)%NAME))//'(' &
                               & //TRIM(ADJUSTL(LINE_LABELS(P)))//')'
                  LINE_LABELS(P)=REPEAT(' ',(LEN(LINE_LABELS(P))-LEN_TRIM(LINE_LABELS(P)))/2) &
                               & //TRIM(ADJUSTL(LINE_LABELS(P)))
               END IF

            END IF
            P=P+1
         END DO ! End of loop over levels (j)
      END DO ! End of loop over levels (i)
   END DO ! End of loop over coolants

!  Specify the output file name
   WRITE(FILENAME,"(I4)") RAY+1
   FILENAME=REPEAT('0',1+FLOOR(LOG10(REAL(NRAYS)))-LEN_TRIM(ADJUSTL(FILENAME)))//TRIM(ADJUSTL(FILENAME))
   FILENAME=TRIM(ADJUSTL(PREFIX))//TRIM(ADJUSTL('.tau_'))//TRIM(ADJUSTL(FILENAME))//TRIM(ADJUSTL(SUFFIX))
   WRITE(6,*) 'Writing file [',TRIM(FILENAME),']'

!  Open and write to the output file
   OPEN(UNIT=1,FILE=FILENAME,STATUS='REPLACE')
   WRITE(1,"('Particle #',999(2X,A11))") LINE_LABELS
   DO P=1,NPART
      WRITE(1,"(I10)",ADVANCE='NO') P
      DO N=1,NCOOL
         DO I=1,COOLANT(N)%NLEVEL
            DO J=1,COOLANT(N)%NLEVEL
               IF(COOLANT(N)%A_COEFF(I,J).EQ.0) CYCLE
               WRITE(1,"(999ES13.5)",ADVANCE='NO') PARTICLE(P)%COOLANT(N)%OPACITY(RAY,I,J)
            END DO
         END DO
      END DO
      WRITE(1,*)
   END DO
   CLOSE(1)

   DEALLOCATE(LINE_LABELS)

   RETURN
END SUBROUTINE WRITE_OPACITIES
!=======================================================================

!=======================================================================
!
!  Write the local emissivities for the emission lines of each coolant
!  and for each particle to the specified output file.
!
!-----------------------------------------------------------------------
SUBROUTINE WRITE_EMISSIVITIES(PREFIX,SUFFIX,NPART,NCOOL,COOLANT,PARTICLE)

   USE HEALPIX_TYPES
   USE COOLANT_MODULE
   USE PARTICLE_MODULE
   USE FUNCTIONS_MODULE

   IMPLICIT NONE

   CHARACTER(LEN=*),    INTENT(IN) :: PREFIX,SUFFIX
   INTEGER(KIND=I4B),   INTENT(IN) :: NPART,NCOOL
   TYPE(COOLANT_TYPE),  INTENT(IN) :: COOLANT(:)
   TYPE(PARTICLE_TYPE), INTENT(IN) :: PARTICLE(:)

   CHARACTER(LEN=256) :: FILENAME
   CHARACTER(LEN=11), ALLOCATABLE :: LINE_LABELS(:)

   INTEGER(KIND=I4B)  :: I,J,N,P

!  Create the transition label for each emission line
   ALLOCATE(LINE_LABELS(1:COUNT((/(COOLANT(N)%A_COEFF,N=1,NCOOL)/).GT.0)))
   P=1
   DO N=1,NCOOL ! Loop over coolants
      DO I=1,COOLANT(N)%NLEVEL ! Loop over levels (i)
         DO J=1,COOLANT(N)%NLEVEL ! Loop over levels (j)
            IF(COOLANT(N)%A_COEFF(I,J).EQ.0) CYCLE

!           Assume all coolants with < 6 levels are atoms with fine-structure lines
            IF(COOLANT(N)%NLEVEL.LT.6) THEN

!              Label fine-structure lines with their wavelength in micron (or in Angstrom, if appropriate)
               IF(COOLANT(N)%FREQUENCY(I,J).GE.3.0D14) THEN
                  WRITE(LINE_LABELS(P),"(I4)") NINT(C/COOLANT(N)%FREQUENCY(I,J)*1.0D8) ! Wavelength in Angstrom (nearest int)
                  LINE_LABELS(P)=TRIM(ADJUSTL(LINE_LABELS(P)))//'A'
               ELSE IF(COOLANT(N)%FREQUENCY(I,J).GE.3.0D13) THEN
                  WRITE(LINE_LABELS(P),"(F5.2)") C/COOLANT(N)%FREQUENCY(I,J)*1.0D4 ! Wavelength in micron (2 dp accuracy)
                  LINE_LABELS(P)=TRIM(ADJUSTL(LINE_LABELS(P)))//'um'
               ELSE IF(COOLANT(N)%FREQUENCY(I,J).GE.1.0D13) THEN
                  WRITE(LINE_LABELS(P),"(F5.1)") C/COOLANT(N)%FREQUENCY(I,J)*1.0D4 ! Wavelength in micron (1 dp accuracy)
                  LINE_LABELS(P)=TRIM(ADJUSTL(LINE_LABELS(P)))//'um'
               ELSE
                  WRITE(LINE_LABELS(P),"(I4)") NINT(C/COOLANT(N)%FREQUENCY(I,J)*1.0D4) ! Wavelength in micron (nearest int)
                  LINE_LABELS(P)=TRIM(ADJUSTL(LINE_LABELS(P)))//'um'
               END IF

!              Handle labelling of both neutral and ionized atomic species
               LINE_LABELS(P)='['//COOLANT(N)%NAME(1:VERIFY(COOLANT(N)%NAME,'+ ',.TRUE.))  &
                               & //REPEAT('I',COUNT_SUBSTRING(COOLANT(N)%NAME,'+'))//'I] ' &
                               & //TRIM(ADJUSTL(LINE_LABELS(P)))

!           Assume all coolants with more levels are molecules with pure rotational or ro-vibrational lines
            ELSE

!              If the line frequency corresponds to a wavelength in the micron range
!              or shorter, then label the line with its wavelength in micron instead
               IF(COOLANT(N)%FREQUENCY(I,J).GE.3.0D14) THEN
                  WRITE(LINE_LABELS(P),"(I4)") NINT(C/COOLANT(N)%FREQUENCY(I,J)*1.0D8) ! Wavelength in Angstrom (nearest int)
                  LINE_LABELS(P)=TRIM(ADJUSTL(COOLANT(N)%NAME))//' '//TRIM(ADJUSTL(LINE_LABELS(P)))//'A'
               ELSE IF(COOLANT(N)%FREQUENCY(I,J).GE.3.0D13) THEN
                  WRITE(LINE_LABELS(P),"(F5.2)") C/COOLANT(N)%FREQUENCY(I,J)*1.0D4 ! Wavelength in micron (2 dp accuracy)
                  LINE_LABELS(P)=TRIM(ADJUSTL(COOLANT(N)%NAME))//' '//TRIM(ADJUSTL(LINE_LABELS(P)))//'um'
               ELSE IF(COOLANT(N)%FREQUENCY(I,J).GE.1.0D13) THEN
                  WRITE(LINE_LABELS(P),"(F5.1)") C/COOLANT(N)%FREQUENCY(I,J)*1.0D4 ! Wavelength in micron (1 dp accuracy)
                  LINE_LABELS(P)=TRIM(ADJUSTL(COOLANT(N)%NAME))//' '//TRIM(ADJUSTL(LINE_LABELS(P)))//'um'
               ELSE IF(COOLANT(N)%FREQUENCY(I,J).GE.6.0D12) THEN
                  WRITE(LINE_LABELS(P),"(I4)") NINT(C/COOLANT(N)%FREQUENCY(I,J)*1.0D4) ! Wavelength in micron (nearest int)
                  LINE_LABELS(P)=TRIM(ADJUSTL(COOLANT(N)%NAME))//' '//TRIM(ADJUSTL(LINE_LABELS(P)))//'um'
               ELSE
                  IF(J.GT.100) THEN
                     WRITE(LINE_LABELS(P),"(I4,'-',I3)") I-1,J-1
                  ELSE IF(J.GT.10) THEN
                     WRITE(LINE_LABELS(P),"(I4,'-',I2)") I-1,J-1
                  ELSE
                     WRITE(LINE_LABELS(P),"(I4,'-',I1)") I-1,J-1
                  END IF
                  LINE_LABELS(P)=TRIM(ADJUSTL(COOLANT(N)%NAME))//'(' &
                               & //TRIM(ADJUSTL(LINE_LABELS(P)))//')'
                  LINE_LABELS(P)=REPEAT(' ',(LEN(LINE_LABELS(P))-LEN_TRIM(LINE_LABELS(P)))/2) &
                               & //TRIM(ADJUSTL(LINE_LABELS(P)))
               END IF

            END IF
            P=P+1
         END DO ! End of loop over levels (j)
      END DO ! End of loop over levels (i)
   END DO ! End of loop over coolants

!  Specify the output file name
   FILENAME=TRIM(ADJUSTL(PREFIX))//TRIM(ADJUSTL('.emis'))//TRIM(ADJUSTL(SUFFIX))
   WRITE(6,*) 'Writing file [',TRIM(FILENAME),']'

!  Open and write to the output file
   OPEN(UNIT=1,FILE=FILENAME,STATUS='REPLACE')
   WRITE(1,"('Particle #',999(2X,A11))") LINE_LABELS
   DO P=1,NPART
      WRITE(1,"(I10)",ADVANCE='NO') P
      DO N=1,NCOOL
         DO I=1,COOLANT(N)%NLEVEL
            DO J=1,COOLANT(N)%NLEVEL
               IF(COOLANT(N)%A_COEFF(I,J).EQ.0) CYCLE
               WRITE(1,"(999ES13.5)",ADVANCE='NO') PARTICLE(P)%COOLANT(N)%EMISSIVITY(I,J)
            END DO
         END DO
      END DO
      WRITE(1,*)
   END DO
   CLOSE(1)

   DEALLOCATE(LINE_LABELS)

   RETURN
END SUBROUTINE WRITE_EMISSIVITIES
!=======================================================================

!=======================================================================
!
!  Write the local cooling rate for each coolant species and for each
!  particle to the specified output file.
!
!-----------------------------------------------------------------------
SUBROUTINE WRITE_COOLING_RATES(PREFIX,SUFFIX,NPART,NCOOL,COOLANT,PARTICLE)

   USE HEALPIX_TYPES
   USE COOLANT_MODULE
   USE PARTICLE_MODULE
   USE FUNCTIONS_MODULE

   IMPLICIT NONE

   CHARACTER(LEN=*),    INTENT(IN) :: PREFIX,SUFFIX
   INTEGER(KIND=I4B),   INTENT(IN) :: NPART,NCOOL
   TYPE(COOLANT_TYPE),  INTENT(IN) :: COOLANT(:)
   TYPE(PARTICLE_TYPE), INTENT(IN) :: PARTICLE(:)

   CHARACTER(LEN=256) :: FILENAME
   CHARACTER(LEN=11), ALLOCATABLE :: COOLING_LABELS(:)

   INTEGER(KIND=I4B)  :: N,P

!  Create the label for each coolant species
   ALLOCATE(COOLING_LABELS(1:NCOOL+1))
   DO N=1,NCOOL
!     Assume all coolants with < 6 levels are atoms with fine-structure lines
      IF(COOLANT(N)%NLEVEL.LT.6) THEN
!        Handle labelling of both neutral and ionized atomic species
         COOLING_LABELS(N)='['//COOLANT(N)%NAME(1:VERIFY(COOLANT(N)%NAME,'+ ',.TRUE.))  &
                            & //REPEAT('I',COUNT_SUBSTRING(COOLANT(N)%NAME,'+'))//'I]'

!     Assume all coolants with more levels are molecules with pure rotational or ro-vibrational lines
      ELSE
         COOLING_LABELS(N)=TRIM(ADJUSTL(COOLANT(N)%NAME))
      END IF
      COOLING_LABELS(N)=REPEAT(' ',(LEN(COOLING_LABELS(N))-LEN_TRIM(COOLING_LABELS(N)))/2) &
                     & //TRIM(ADJUSTL(COOLING_LABELS(N)))
   END DO
   COOLING_LABELS(NCOOL+1)='   Total   '

!  Specify the output file name
   FILENAME=TRIM(ADJUSTL(PREFIX))//TRIM(ADJUSTL('.cool'))//TRIM(ADJUSTL(SUFFIX))
   WRITE(6,*) 'Writing file [',TRIM(FILENAME),']'

!  Open and write to the output file
   OPEN(UNIT=1,FILE=FILENAME,STATUS='REPLACE')
   WRITE(1,"('Particle #',500(2X,A11))") COOLING_LABELS
   DO P=1,NPART
      WRITE(1,"(I10,500ES13.5)") P,PARTICLE(P)%COOLING_RATE,SUM(PARTICLE(P)%COOLING_RATE)
   END DO
   CLOSE(1)

   DEALLOCATE(COOLING_LABELS)

   RETURN
END SUBROUTINE WRITE_COOLING_RATES
!=======================================================================

!=======================================================================
!
!  Write the local heating rate for each heating mechanism and for each
!  particle to the specified output file.
!
!-----------------------------------------------------------------------
SUBROUTINE WRITE_HEATING_RATES(PREFIX,SUFFIX,NPART,NHEAT,PARTICLE)

   USE HEALPIX_TYPES
   USE PARTICLE_MODULE

   IMPLICIT NONE

   CHARACTER(LEN=*),    INTENT(IN) :: PREFIX,SUFFIX
   INTEGER(KIND=I4B),   INTENT(IN) :: NPART,NHEAT
   TYPE(PARTICLE_TYPE), INTENT(IN) :: PARTICLE(:)

   CHARACTER(LEN=256) :: FILENAME
   CHARACTER(LEN=11), ALLOCATABLE :: HEATING_LABELS(:)

   INTEGER(KIND=I4B)  :: I,P

!  Create the label for each heating mechanism
   ALLOCATE(HEATING_LABELS(1:NHEAT+1))
   HEATING_LABELS(1) =' Dust P.E. '
   HEATING_LABELS(2) ='  H2 Form  '
   HEATING_LABELS(3) ='H2* Pumping'
   HEATING_LABELS(4) ='H2 Photodis'
   HEATING_LABELS(5) ='CI Photoion'
   HEATING_LABELS(6) ='Cosmic-Rays'
   HEATING_LABELS(7) ='Turbulence '
   HEATING_LABELS(8) =' Chemistry '
   HEATING_LABELS(9) =' Gas-Grain '
   HEATING_LABELS(10)='  Coulomb  '
   HEATING_LABELS(11)='   Total   '

!  Specify the output file name
   FILENAME=TRIM(ADJUSTL(PREFIX))//TRIM(ADJUSTL('.heat'))//TRIM(ADJUSTL(SUFFIX))
   WRITE(6,*) 'Writing file [',TRIM(FILENAME),']'

!  Open and write to the output file
   OPEN(UNIT=1,FILE=FILENAME,STATUS='REPLACE')
   WRITE(1,"('Particle #',500(2X,A11))") HEATING_LABELS
   DO P=1,NPART
      WRITE(1,"(I10,500ES13.5)") P,PARTICLE(P)%HEATING_RATE,SUM(PARTICLE(P)%HEATING_RATE)
   END DO
   CLOSE(1)

   DEALLOCATE(HEATING_LABELS)

   RETURN
END SUBROUTINE WRITE_HEATING_RATES
!=======================================================================

!=======================================================================
!
!  Write the gas temperatures and thermal balance iteration parameters
!  for each particle to the specified output file.
!
!-----------------------------------------------------------------------
SUBROUTINE WRITE_THERMAL_BALANCE(PREFIX,SUFFIX,NPART,PARTICLE)

   USE HEALPIX_TYPES
   USE PARTICLE_MODULE

   IMPLICIT NONE

   CHARACTER(LEN=*),    INTENT(IN) :: PREFIX,SUFFIX
   INTEGER(KIND=I4B),   INTENT(IN) :: NPART
   TYPE(PARTICLE_TYPE), INTENT(IN) :: PARTICLE(:)

   CHARACTER(LEN=256) :: FILENAME
   INTEGER(KIND=I4B)  :: P

!  Specify the output file name
   FILENAME=TRIM(ADJUSTL(PREFIX))//TRIM(ADJUSTL('.temp'))//TRIM(ADJUSTL(SUFFIX))
   WRITE(6,*) 'Writing file [',TRIM(FILENAME),']'

!  Open and write to the output file
   OPEN(UNIT=1,FILE=FILENAME,STATUS='REPLACE')
   WRITE(1,"('Particle #   T_old (K)    T_new (K)     Heating      Cooling    Difference    F_ratio  Converged?')")
   DO P=1,NPART
      IF(PARTICLE(P)%TEMPERATURE_CONVERGED) THEN
         WRITE(1,"(I10,4ES13.5,SP,ES13.5,SS,F9.3,'%  ---Yes----')") &
            &     P,PARTICLE(P)%PREVIOUS_TEMPERATURE,PARTICLE(P)%GAS_TEMPERATURE, &
            &     SUM(PARTICLE(P)%HEATING_RATE),SUM(PARTICLE(P)%COOLING_RATE),  &
            &     SUM(PARTICLE(P)%HEATING_RATE)-SUM(PARTICLE(P)%COOLING_RATE),  &
            & ABS(SUM(PARTICLE(P)%HEATING_RATE)-SUM(PARTICLE(P)%COOLING_RATE))/ &
            & ABS(SUM(PARTICLE(P)%HEATING_RATE)+SUM(PARTICLE(P)%COOLING_RATE))*200
      ELSE
         WRITE(1,"(I10,4ES13.5,SP,ES13.5,SS,F9.3,'%  <<  No  >>')") &
            &     P,PARTICLE(P)%PREVIOUS_TEMPERATURE,PARTICLE(P)%GAS_TEMPERATURE, &
            &     SUM(PARTICLE(P)%HEATING_RATE),SUM(PARTICLE(P)%COOLING_RATE),  &
            &     SUM(PARTICLE(P)%HEATING_RATE)-SUM(PARTICLE(P)%COOLING_RATE),  &
            & ABS(SUM(PARTICLE(P)%HEATING_RATE)-SUM(PARTICLE(P)%COOLING_RATE))/ &
            & ABS(SUM(PARTICLE(P)%HEATING_RATE)+SUM(PARTICLE(P)%COOLING_RATE))*200
      ENDIF
   END DO
   CLOSE(1)

   RETURN
END SUBROUTINE WRITE_THERMAL_BALANCE
!=======================================================================

!=======================================================================
!
!  Write the convergence status of the level populations of each coolant
!  and for the thermal balance for each particle to the specified output
!  file.
!
!-----------------------------------------------------------------------
SUBROUTINE WRITE_CONVERGENCE_STATUS(PREFIX,SUFFIX,NPART,NCOOL,COOLANT,PARTICLE)

   USE HEALPIX_TYPES
   USE COOLANT_MODULE
   USE PARTICLE_MODULE
   USE FUNCTIONS_MODULE

   IMPLICIT NONE

   CHARACTER(LEN=*),    INTENT(IN) :: PREFIX,SUFFIX
   INTEGER(KIND=I4B),   INTENT(IN) :: NPART,NCOOL
   TYPE(COOLANT_TYPE),  INTENT(IN) :: COOLANT(:)
   TYPE(PARTICLE_TYPE), INTENT(IN) :: PARTICLE(:)

   CHARACTER(LEN=256) :: FILENAME
   CHARACTER(LEN=11), ALLOCATABLE :: COOLING_LABELS(:)

   INTEGER(KIND=I4B)  :: N,P

!  Create the label for each coolant species
   ALLOCATE(COOLING_LABELS(1:NCOOL+1))
   DO N=1,NCOOL
!     Assume all coolants with < 6 levels are atoms with fine-structure lines
      IF(COOLANT(N)%NLEVEL.LT.6) THEN
!        Handle labelling of both neutral and ionized atomic species
         COOLING_LABELS(N)='['//COOLANT(N)%NAME(1:VERIFY(COOLANT(N)%NAME,'+ ',.TRUE.))  &
                            & //REPEAT('I',COUNT_SUBSTRING(COOLANT(N)%NAME,'+'))//'I]'

!     Assume all coolants with more levels are molecules with pure rotational or ro-vibrational lines
      ELSE
         COOLING_LABELS(N)=TRIM(ADJUSTL(COOLANT(N)%NAME))
      END IF
      COOLING_LABELS(N)=REPEAT(' ',(LEN(COOLING_LABELS(N))-LEN_TRIM(COOLING_LABELS(N)))/2) &
                     & //TRIM(ADJUSTL(COOLING_LABELS(N)))
   END DO
   COOLING_LABELS(NCOOL+1)='   T_gas   '

!  Specify the output file name
   FILENAME=TRIM(ADJUSTL(PREFIX))//TRIM(ADJUSTL('.conv'))//TRIM(ADJUSTL(SUFFIX))
   WRITE(6,*) 'Writing file [',TRIM(FILENAME),']'

!  Open and write to the output file
   OPEN(UNIT=1,FILE=FILENAME,STATUS='REPLACE')
   WRITE(1,"('Particle #',500(2X,A11))") COOLING_LABELS
   DO P=1,NPART
      WRITE(1,"(I10)",ADVANCE='NO') P
      DO N=1,NCOOL
         IF(PARTICLE(P)%COOLANT(N)%CONVERGED) THEN
            WRITE(1,"('  ----Yes----')",ADVANCE='NO')
         ELSE
            WRITE(1,"('      No     ')",ADVANCE='NO')
         END IF
      END DO
      IF(PARTICLE(P)%TEMPERATURE_CONVERGED) THEN
         WRITE(1,"('  ----Yes----')",ADVANCE='NO')
      ELSE
         WRITE(1,"('      No     ')",ADVANCE='NO')
      END IF
      WRITE(1,*)
   END DO
   CLOSE(1)

   DEALLOCATE(COOLING_LABELS)

   RETURN
END SUBROUTINE WRITE_CONVERGENCE_STATUS
!=======================================================================
