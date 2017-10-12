!=======================================================================
!
!  Read in the coolant species, their energy level structure, transition
!  properties and collisional rate coefficients. The specified files are
!  assumed to contain entries in the LAMDA/RADEX format, allowing files
!  to be downloaded directly from the online database.
!
!-----------------------------------------------------------------------
SUBROUTINE READ_COOLANTS(NCOOL,COOLANT)

   USE HEALPIX_TYPES
   USE COOLANT_MODULE

   IMPLICIT NONE

   INTEGER(KIND=I4B),  INTENT(IN)    :: NCOOL
   TYPE(COOLANT_TYPE), INTENT(INOUT) :: COOLANT(*)

   INTEGER(KIND=I4B) :: I,J,K,L,M,N,INDEX,IER
   INTEGER(KIND=I4B) :: NLEVEL,NLINE,NTEMP,NPARTNER,NCOLL,PARTNER_ID,MAX_NTEMP

   DO N=1,NCOOL ! Loop over coolants

!     Open the input file
      OPEN(UNIT=1,FILE='Datafiles/Collisional-Rates/'//COOLANT(N)%FILENAME,IOSTAT=IER,ACTION='READ',STATUS='OLD')
      READ(1,*,IOSTAT=IER) ! Skip the first comment line

!     Produce an error message if the file does not exist (or cannot be opened for whatever reason)
      IF(IER.NE.0) THEN
         WRITE(6,*) 'ERROR! Cannot open coolant data file ',TRIM(COOLANT(N)%FILENAME),' for input'
         WRITE(6,*)
         CLOSE(1)
         STOP
      END IF

      READ(1,*,IOSTAT=IER) COOLANT(N)%NAME ! Read the name of the coolant
      READ(1,*,IOSTAT=IER)
      READ(1,*,IOSTAT=IER) COOLANT(N)%MOLECULAR_MASS ! Read the molecular mass
      READ(1,*,IOSTAT=IER)
      COOLANT(N)%INDEX=0 ! Initialize the coolant species index (assigned later)

!     Read the number of levels and allocate the energy, statistical weight,
!     Einstein A & B coefficient and transition frequency arrays accordingly
      READ(1,*,IOSTAT=IER) NLEVEL
      READ(1,*,IOSTAT=IER)
      IF(NLEVEL.LT.2) THEN
         WRITE(6,*) 'ERROR! Incorrect number of energy levels in coolant data file ',TRIM(COOLANT(N)%FILENAME),' (NLEVEL=',NLEVEL,')'
         CLOSE(1)
         STOP
      END IF
      COOLANT(N)%NLEVEL=NLEVEL
      ALLOCATE(COOLANT(N)%ENERGY(1:NLEVEL),COOLANT(N)%WEIGHT(1:NLEVEL))
      ALLOCATE(COOLANT(N)%A_COEFF(1:NLEVEL,1:NLEVEL))
      ALLOCATE(COOLANT(N)%B_COEFF(1:NLEVEL,1:NLEVEL))
      ALLOCATE(COOLANT(N)%FREQUENCY(1:NLEVEL,1:NLEVEL))

!     Initialize the level energies, statistical weights,
!     Einstein coefficients and transition frequencies
      COOLANT(N)%ENERGY=0.0D0
      COOLANT(N)%WEIGHT=0.0D0
      COOLANT(N)%A_COEFF=0.0D0
      COOLANT(N)%B_COEFF=0.0D0
      COOLANT(N)%FREQUENCY=0.0D0

!     Read the energy (cm^-1) and statistical weight of each level
      DO L=1,NLEVEL ! Loop over levels
         READ(1,*,IOSTAT=IER) I,COOLANT(N)%ENERGY(I),COOLANT(N)%WEIGHT(I)
         COOLANT(N)%ENERGY(I)=COOLANT(N)%ENERGY(I)*C*HP ! Convert from cm^-1 to erg
      END DO ! End of loop over levels
      READ(1,*,IOSTAT=IER)

!     Read the Einstein A coefficient (s^-1) and frequency (GHz) of each radiative transition
      READ(1,*,IOSTAT=IER) NLINE
      READ(1,*,IOSTAT=IER)
      DO L=1,NLINE ! Loop over radiative transitions
         READ(1,*,IOSTAT=IER) INDEX,I,J,COOLANT(N)%A_COEFF(I,J),COOLANT(N)%FREQUENCY(I,J)
         COOLANT(N)%FREQUENCY(I,J)=COOLANT(N)%FREQUENCY(I,J)*1.0D9 ! Convert from GHz to Hz
!        Calculate the Einstein B coefficient using B_ij = A_ij/(2.h.nu^3/c^2)
         COOLANT(N)%B_COEFF(I,J)=COOLANT(N)%A_COEFF(I,J)/(2*HP*(COOLANT(N)%FREQUENCY(I,J)**3)/(C**2))
!        Calculate the Einstein B coefficient for the reverse transition from detailed balance
         COOLANT(N)%B_COEFF(J,I)=COOLANT(N)%B_COEFF(I,J)*(COOLANT(N)%WEIGHT(I)/COOLANT(N)%WEIGHT(J))
      END DO ! End of loop over radiative transitions
      READ(1,*,IOSTAT=IER)

!     Calculate the transition frequencies between all levels (even if forbidden)
      DO I=1,NLEVEL
         DO J=1,NLEVEL
!           Check that the calculated and measured frequencies differ by <0.1%
!           Produce an error message if the difference between them is greater
            IF(COOLANT(N)%FREQUENCY(I,J).NE.0.0D0) THEN
               IF(ABS(COOLANT(N)%FREQUENCY(I,J)-ABS(COOLANT(N)%ENERGY(I)-COOLANT(N)%ENERGY(J))/HP) &
                   & /COOLANT(N)%FREQUENCY(I,J).GT.1.0D-3) THEN
                  WRITE(6,*) 'ERROR! Calculated frequency differs from measured frequency by >0.1%'
                  WRITE(6,"(1PD12.5,'Hz vs',1PD12.5,'Hz')") ABS(COOLANT(N)%ENERGY(I)-COOLANT(N)%ENERGY(J))/HP, &
                                                          & COOLANT(N)%FREQUENCY(I,J)
                  CLOSE(1)
                  STOP
               END IF
            ELSE
               COOLANT(N)%FREQUENCY(I,J)=ABS(COOLANT(N)%ENERGY(I)-COOLANT(N)%ENERGY(J))/HP
            END IF
         END DO
      END DO

!     Allocate and initialize the collisional rate coefficient arrays
!     allowing a maximum of 1000 temperature values per collision partner
      MAX_NTEMP=1000
      COOLANT(N)%NTEMP=MAX_NTEMP
      ALLOCATE(COOLANT(N)%TEMPERATURE(1:7,1:MAX_NTEMP))
      ALLOCATE(COOLANT(N)%C_COEFF(1:7,1:NLEVEL,1:NLEVEL,1:MAX_NTEMP))
      COOLANT(N)%TEMPERATURE=0.0D0
      COOLANT(N)%C_COEFF=0.0D0

!     Read the collisional rate coefficients (cm^3 s^-1) for each collision partner
      MAX_NTEMP=0
      READ(1,*,IOSTAT=IER) NPARTNER
      DO L=1,NPARTNER ! Loop over collision partners
         READ(1,*,IOSTAT=IER)
         READ(1,*,IOSTAT=IER) PARTNER_ID
         IF(PARTNER_ID.LT.1 .OR. PARTNER_ID.GT.7) THEN
            WRITE(6,*) 'ERROR! Unrecognized collision partner ID in coolant data file ',TRIM(COOLANT(N)%FILENAME),' (ID=',PARTNER_ID,')'
            CLOSE(1)
            STOP
         END IF
         READ(1,*,IOSTAT=IER)
         READ(1,*,IOSTAT=IER) NCOLL
         READ(1,*,IOSTAT=IER)
         READ(1,*,IOSTAT=IER) NTEMP ; MAX_NTEMP=MAX(MAX_NTEMP,NTEMP)
         IF(NTEMP.GT.1000) THEN
            WRITE(6,*) 'ERROR! Number of temperature values exceeds limit in coolant data file ',TRIM(COOLANT(N)%FILENAME),' (NTEMP=',NTEMP,')'
            CLOSE(1)
            STOP
         END IF
         READ(1,*,IOSTAT=IER)
         READ(1,*,IOSTAT=IER) (COOLANT(N)%TEMPERATURE(PARTNER_ID,K),K=1,NTEMP)
         READ(1,*,IOSTAT=IER)
         DO M=1,NCOLL ! Loop over collisional transitions
            READ(1,*,IOSTAT=IER) INDEX,I,J,(COOLANT(N)%C_COEFF(PARTNER_ID,I,J,K),K=1,NTEMP)
!           Calculate the reverse (excitation) rate coefficients from
!           detailed balance: C_ji = C_ij*gi/gj*exp(-(Ei-Ej)/kT)
            DO K=1,NTEMP ! Loop over temperatures
               IF(COOLANT(N)%C_COEFF(PARTNER_ID,I,J,K).NE.0.0D0 .AND. &
                & COOLANT(N)%C_COEFF(PARTNER_ID,J,I,K).EQ.0.0D0) THEN
                  COOLANT(N)%C_COEFF(PARTNER_ID,J,I,K)=COOLANT(N)%C_COEFF(PARTNER_ID,I,J,K) &
                                                    & *(COOLANT(N)%WEIGHT(I)/COOLANT(N)%WEIGHT(J)) &
                                                    & *EXP(-(COOLANT(N)%ENERGY(I)-COOLANT(N)%ENERGY(J)) &
                                                    &      /(KB*COOLANT(N)%TEMPERATURE(PARTNER_ID,K)))
               END IF
            END DO ! End of loop over temperatures
         END DO ! End of loop over collisional transitions
      END DO ! End of loop over collision partners

      COOLANT(N)%NTEMP=MAX_NTEMP

      CLOSE(1)

   END DO ! End of loop over coolants
   
   RETURN
END SUBROUTINE READ_COOLANTS
!=======================================================================
