!=======================================================================
!
!  Calculate the rate coefficients for all reactions at the specified
!  temperature and visual extinction A_V. The photodissociation of H2
!  and CO and the photoionization of CI and SI are treated separately
!  in detail (see the routines in photorates.f90). Multiple rates for
!  the same reaction (duplicates) are allowed in the ratefile and are
!  activated based on their minimum and maximum temperature specified
!  in that file. Negative gamma factors are ignored below the minimum
!  temperature at which the reaction rate is valid.
!
!  X-ray induced reaction rates are calculated following the detailed
!  treatment of Meijerink & Spaans (2005, A&A, 436, 397).
!
!-----------------------------------------------------------------------
SUBROUTINE CALCULATE_REACTION_RATES(NRAYS,NSPEC,NREAC,NXSEC,CROSS_SECTION, &
                                  & GAS_TEMPERATURE,DUST_TEMPERATURE,FUV_FIELD, &
                                  & XRAY_FIELD,FUV_SURFACE,AV,COLUMN_DENSITY, &
                                  & REACTANT,PRODUCT,ALPHA,BETA,GAMMA,RTMIN,RTMAX, &
                                  & DUPLICATE,RATE)

   USE HEALPIX_TYPES
   USE GLOBAL_MODULE
   USE FUNCTIONS_MODULE
   USE NUM2STR_FUNCTION
   USE CROSS_SECTION_MODULE

   IMPLICIT NONE

   INTEGER(KIND=I4B),        INTENT(IN)  :: NRAYS,NSPEC,NREAC,NXSEC
   TYPE(CROSS_SECTION_TYPE), INTENT(IN)  :: CROSS_SECTION(:)
   REAL(KIND=DP),            INTENT(IN)  :: GAS_TEMPERATURE,DUST_TEMPERATURE
   REAL(KIND=DP),            INTENT(IN)  :: FUV_FIELD,XRAY_FIELD
   REAL(KIND=DP),            INTENT(IN)  :: FUV_SURFACE(0:),AV(0:),COLUMN_DENSITY(0:,:)
   CHARACTER(LEN=10),        INTENT(IN)  :: REACTANT(:,:),PRODUCT(:,:)
   REAL(KIND=DP),            INTENT(IN)  :: ALPHA(:),BETA(:),GAMMA(:)
   REAL(KIND=DP),            INTENT(IN)  :: RTMIN(:),RTMAX(:)
   INTEGER(KIND=I4B),        INTENT(IN)  :: DUPLICATE(:)
   REAL(KIND=DP),            INTENT(OUT) :: RATE(:)

   INTEGER(KIND=I4B) :: I,J,K
   REAL(KIND=DP)     :: PHI_PAH,CION,STICKING,FLUX,YIELD

!  Initialize the rate coefficients
   RATE=0.0D0

!!$!  Calculate the photoreaction rates for species with available cross section data
!!$   DO I=1,NXSEC ! Loop over available cross sections
!!$      IF(CROSS_SECTION(I)%INDEX.NE.0) THEN
!!$         DO K=0,NRAYS-1 ! Loop over all rays
!!$            RATE(CROSS_SECTION(I)%INDEX) = RATE(CROSS_SECTION(I)%INDEX) + PHOTOREACTION_RATE(FUV_SURFACE(K),AV(K),CROSS_SECTION(I))
!!$         END DO ! End of loop over all rays
!!$      END IF
!!$   END DO ! End of loop over available cross sections

!  Loop over all reactions
   DO I=1,NREAC

!     Determine the type of reaction
      IF(REACTANT(I,2).EQ."PHOTON") GOTO 1
      IF(REACTANT(I,2).EQ."CRP   ") GOTO 2
      IF(REACTANT(I,2).EQ."CRPHOT") GOTO 3
      IF(REACTANT(I,2).EQ."XRAY  ") GOTO 4
      IF(REACTANT(I,2).EQ."XRSEC ") GOTO 5
      IF(REACTANT(I,2).EQ."XRLYA ") GOTO 6
      IF(REACTANT(I,2).EQ."XRPHOT") GOTO 6
      IF(REACTANT(I,2).EQ."FREEZE") GOTO 7
      IF(REACTANT(I,2).EQ."CRH   ") GOTO 8
      IF(REACTANT(I,2).EQ."PHOTD ") GOTO 9
      IF(REACTANT(I,2).EQ."THERM ") GOTO 10
      IF(REACTANT(I,2)(1:1).EQ."#") GOTO 11

!-----------------------------------------------------------------------

!     Thermal reactions:

!     The rate of H2 formation on grains is calculated separately
!     by the function H2_FORMATION_RATE (see function for details)
      IF((REACTANT(I,1).EQ."H  " .AND. REACTANT(I,2).EQ."H  " .AND. REACTANT(I,3).EQ."#  ") .AND. &
       & (PRODUCT(I,1).EQ."H2  " .AND. PRODUCT(I,2).EQ."#   ")) THEN
         RATE(I)=H2_FORMATION_RATE(GAS_TEMPERATURE,DUST_TEMPERATURE)
         GOTO 20
      END IF

!     Rates for reactions involving PAHs are calculated according to the
!     treatment of Wolfire et al. (2003, ApJ, 587, 278; 2008, ApJ, 680, 384)
      IF(ANY(REACTANT(I,:).EQ."PAH  ") .OR. ANY(REACTANT(I,:).EQ."PAH0 ") .OR. &
       & ANY(REACTANT(I,:).EQ."PAH+ ") .OR. ANY(REACTANT(I,:).EQ."PAH- ")) THEN
         PHI_PAH=0.4D0
         RATE(I)=ALPHA(I)*(GAS_TEMPERATURE/100.0D0)**BETA(I)*PHI_PAH
         GOTO 20
      END IF

!     Check for large negative gamma values that might cause discrepant
!     rates at low temperatures. Set these rates to zero when T < RTMIN
      IF(DUPLICATE(I).EQ.0) THEN
         IF(GAMMA(I).LT.-200.0D0 .AND. GAS_TEMPERATURE.LT.RTMIN(I)) THEN
            RATE(I)=0.0D0
         ELSE
            RATE(I)=ALPHA(I)*(GAS_TEMPERATURE/300.0D0)**BETA(I)*EXP(-(GAMMA(I)/GAS_TEMPERATURE))
         END IF
      ELSE IF(DUPLICATE(I).EQ.1) THEN
         J=I
         DO
            IF(GAS_TEMPERATURE.LE.RTMAX(J)) THEN
               IF(GAMMA(J).LT.-200.0D0 .AND. GAS_TEMPERATURE.LT.RTMIN(J)) THEN
                  RATE(J)=0.0D0
               ELSE
                  RATE(J)=ALPHA(J)*(GAS_TEMPERATURE/300.0D0)**BETA(J)*EXP(-(GAMMA(J)/GAS_TEMPERATURE))
               END IF
               EXIT
            ELSE IF(DUPLICATE(J+1).LT.DUPLICATE(J)) THEN
               IF(GAMMA(J).LT.-200.0D0 .AND. GAS_TEMPERATURE.LT.RTMIN(J)) THEN
                  RATE(J)=0.0D0
               ELSE
                  RATE(J)=ALPHA(J)*(GAS_TEMPERATURE/300.0D0)**BETA(J)*EXP(-(GAMMA(J)/GAS_TEMPERATURE))
               END IF
               EXIT
            ELSE
               RATE(J)=0.0D0
               J=J+1
            END IF
         END DO
      END IF
      GOTO 20

!-----------------------------------------------------------------------

!     FUV photoreactions:

!     The rate of H2 photodissociation is calculated separately by the function
!     H2_PHOTODISSOCIATION_RATE (see photoreaction_rates.f90 for details)
 1    IF(REACTANT(I,1).EQ."H2 " .AND. REACTANT(I,3).EQ."   ") THEN
         DO K=0,NRAYS-1 ! Loop over all rays
            RATE(I)=RATE(I) + H2_PHOTODISSOCIATION_RATE(ALPHA(I),FUV_SURFACE(K),AV(K),COLUMN_DENSITY(K,NH2))
         END DO
         GOTO 20
      END IF

!     The rate of HD photodissociation is calculated separately by the function
!     H2_PHOTODISSOCIATION_RATE (see photoreaction_rates.f90 for details)
      IF(REACTANT(I,1).EQ."HD " .AND. REACTANT(I,3).EQ."   ") THEN
         DO K=0,NRAYS-1 ! Loop over all rays
            RATE(I)=RATE(I) + H2_PHOTODISSOCIATION_RATE(ALPHA(I),FUV_SURFACE(K),AV(K),COLUMN_DENSITY(K,NHD))
         END DO
         GOTO 20
      END IF

!     The rate of CO photodissociation is calculated separately by the function
!     CO_PHOTODISSOCIATION_RATE (see photoreaction_rates.f90 for details)
      IF(REACTANT(I,1).EQ."CO " .AND. REACTANT(I,3).EQ."   " .AND. ANY(PRODUCT(I,:).EQ."C ") .AND. ANY(PRODUCT(I,:).EQ."O ")) THEN
         DO K=0,NRAYS-1 ! Loop over all rays
            RATE(I)=RATE(I) + CO_PHOTODISSOCIATION_RATE(ALPHA(I),FUV_SURFACE(K),AV(K),COLUMN_DENSITY(K,NCO),COLUMN_DENSITY(K,NH2))
         END DO
         GOTO 20
      END IF

!     The rate of CI photoionization is calculated separately by the function
!     CI_PHOTOIONIZATION_RATE (see photoreaction_rates.f90 for details)
      IF(REACTANT(I,1).EQ."C  " .AND. REACTANT(I,3).EQ."   ") THEN
         DO K=0,NRAYS-1 ! Loop over all rays
            RATE(I)=RATE(I) + CI_PHOTOIONIZATION_RATE(ALPHA(I),FUV_SURFACE(K),AV(K),GAMMA(I),COLUMN_DENSITY(K,NC),COLUMN_DENSITY(K,NH2),GAS_TEMPERATURE)
         END DO
         GOTO 20
      END IF

!     The rate of SI photoionization is calculated separately by the function
!     SI_PHOTOIONIZATION_RATE (see photoreaction_rates.f90 for details)
      IF(REACTANT(I,1).EQ."S  " .AND. REACTANT(I,3).EQ."   ") THEN
         DO K=0,NRAYS-1 ! Loop over all rays
            RATE(I)=RATE(I) + SI_PHOTOIONIZATION_RATE(ALPHA(I),FUV_SURFACE(K),AV(K),GAMMA(I),COLUMN_DENSITY(K,NS))
         END DO
         GOTO 20
      END IF

      IF(DUPLICATE(I).EQ.0) THEN
         IF(RATE(I).EQ.0.0D0) THEN
            DO K=0,NRAYS-1 ! Loop over all rays
                RATE(I)=RATE(I) + ALPHA(I)*FUV_SURFACE(K)*EXP(-(GAMMA(I)*AV(K)))
            END DO
         END IF
      ELSE IF(DUPLICATE(I).EQ.1) THEN
         J=I
         DO
            IF(GAS_TEMPERATURE.LE.RTMAX(J)) THEN
               IF(RATE(J).EQ.0.0D0) THEN
                  DO K=0,NRAYS-1 ! Loop over all rays
                     RATE(J)=RATE(J) + ALPHA(J)*FUV_SURFACE(K)*EXP(-(GAMMA(J)*AV(K)))
                  END DO
               END IF
               EXIT
            ELSE IF(DUPLICATE(J+1).LT.DUPLICATE(J)) THEN
               IF(RATE(J).EQ.0.0D0) THEN
                  DO K=0,NRAYS-1 ! Loop over all rays
                     RATE(J)=RATE(J) + ALPHA(J)*FUV_SURFACE(K)*EXP(-(GAMMA(J)*AV(K)))
                  END DO
               END IF
               EXIT
            ELSE
               RATE(J)=0.0D0
               J=J+1
            END IF
         END DO
      END IF
      GOTO 20

!-----------------------------------------------------------------------

!     Cosmic-ray induced ionization:

 2    IF(DUPLICATE(I).EQ.0) THEN
         RATE(I)=ALPHA(I)*ZETA
      ELSE IF(DUPLICATE(I).EQ.1) THEN
         J=I
         DO
            IF(GAS_TEMPERATURE.LE.RTMAX(J)) THEN
               RATE(J)=ALPHA(J)*ZETA
               EXIT
            ELSE IF(DUPLICATE(J+1).LT.DUPLICATE(J)) THEN
               RATE(J)=ALPHA(J)*ZETA
               EXIT
            ELSE
               RATE(J)=0.0D0
               J=J+1
            END IF
         END DO
      END IF
      GOTO 20

!-----------------------------------------------------------------------

!     Photoreactions due to cosmic-ray induced secondary photons:

 3    IF(DUPLICATE(I).EQ.0) THEN
         RATE(I)=ALPHA(I)*ZETA*(GAS_TEMPERATURE/300.0D0)**BETA(I)*GAMMA(I)/(1.0D0-OMEGA)
      ELSE IF(DUPLICATE(I).EQ.1) THEN
         J=I
         DO
            IF(GAS_TEMPERATURE.LE.RTMAX(J)) THEN
               RATE(J)=ALPHA(J)*ZETA*(GAS_TEMPERATURE/300.0D0)**BETA(J)*GAMMA(J)/(1.0D0-OMEGA)
               EXIT
            ELSE IF(DUPLICATE(J+1).LT.DUPLICATE(J)) THEN
               RATE(J)=ALPHA(J)*ZETA*(GAS_TEMPERATURE/300.0D0)**BETA(J)*GAMMA(J)/(1.0D0-OMEGA)
               EXIT
            ELSE
               RATE(J)=0.0D0
               J=J+1
            END IF
         END DO
      END IF
      GOTO 20

!-----------------------------------------------------------------------

!     X-ray photoreactions:

 4    IF(DUPLICATE(I).EQ.0) THEN
         RATE(I)=ALPHA(I)
      ELSE IF(DUPLICATE(I).EQ.1) THEN
         J=I
         DO
            IF(GAS_TEMPERATURE.LE.RTMAX(J)) THEN
               RATE(J)=ALPHA(J)
               EXIT
            ELSE IF(DUPLICATE(J+1).LT.DUPLICATE(J)) THEN
               RATE(J)=ALPHA(J)
               EXIT
            ELSE
               RATE(J)=0.0D0
               J=J+1
            END IF
         END DO
      END IF
      GOTO 20


!-----------------------------------------------------------------------

!     X-ray induced secondary ionization:

 5    IF(DUPLICATE(I).EQ.0) THEN
         RATE(I)=ALPHA(I)*XRAY_FIELD*(1.0D0/eV) ! Convert from erg to eV
      ELSE IF(DUPLICATE(I).EQ.1) THEN
         J=I
         DO
            IF(GAS_TEMPERATURE.LE.RTMAX(J)) THEN
               RATE(J)=ALPHA(J)*XRAY_FIELD*(1.0D0/eV)
               EXIT
            ELSE IF(DUPLICATE(J+1).LT.DUPLICATE(J)) THEN
               RATE(J)=ALPHA(J)*XRAY_FIELD*(1.0D0/eV)
               EXIT
            ELSE
               RATE(J)=0.0D0
               J=J+1
            END IF
         END DO
      END IF
      GOTO 20

!-----------------------------------------------------------------------

!     Photoreactions due to X-ray induced secondary photons (both Lyman-alpha and Lyman-Werner):
!     For Lyman-alpha photons:  k_i = epsilon_Lya * P_i,Lya * x_H * zeta_H / (1 - omega)  [s^-1]
!     For Lyman-Werner photons: k_i = epsilon_LyW * P_i,LyW * x_H2* zeta_H2/ (1 - omega)  [s^-1]
!
!     The value of epsilon_Lya or epsilon_LyW is stored in the alpha parameter of each reaction
!     The value of P_i,Lya or P_i,LyW is stored in the gamma parameter of each reaction

 6    IF(DUPLICATE(I).EQ.0) THEN
         RATE(I)=ALPHA(I)*XRAY_FIELD*(1.0D0/eV)*(GAS_TEMPERATURE/300.0D0)**BETA(I)*GAMMA(I)/(1.0D0-OMEGA)
      ELSE IF(DUPLICATE(I).EQ.1) THEN
         J=I
         DO
            IF(GAS_TEMPERATURE.LE.RTMAX(J)) THEN
               RATE(J)=ALPHA(J)*XRAY_FIELD*(1.0D0/eV)*(GAS_TEMPERATURE/300.0D0)**BETA(J)*GAMMA(J)/(1.0D0-OMEGA)
               EXIT
            ELSE IF(DUPLICATE(J+1).LT.DUPLICATE(J)) THEN
               RATE(J)=ALPHA(J)*XRAY_FIELD*(1.0D0/eV)*(GAS_TEMPERATURE/300.0D0)**BETA(J)*GAMMA(J)/(1.0D0-OMEGA)
               EXIT
            ELSE
               RATE(J)=0.0D0
               J=J+1
            END IF
         END DO
      END IF
      GOTO 20

!-----------------------------------------------------------------------

!     Freeze-out of neutral species and singly charged positive ions:

 7    IF(BETA(I).EQ.0.0D0) THEN
         CION=1.0D0
      ELSE IF(BETA(I).EQ.1.0D0) THEN
         CION=1.0D0+16.71D-4/(GRAIN_RADIUS*GAS_TEMPERATURE)
      ELSE
         CION=0.0D0
      END IF
      STICKING=0.3D0
      RATE(I)=ALPHA(I)*4.57D4*2.4D-22*SQRT(GAS_TEMPERATURE/GAMMA(I))*CION*STICKING
      GOTO 20

!-----------------------------------------------------------------------

!     Desorption due to cosmic-ray heating:

!!$!     Treatment of Hasegawa & Herbst (1993, MNRAS, 261, 83, Equation 15)
!!$ 8    RATE(I)=ALPHA(I)*ZETA

!     Treatment of Roberts et al. (2007, MNRAS, 382, 773, Equation 3)
 8    IF(GAMMA(I).LE.1210.0D0) THEN
         YIELD=1.0D5 ! Number of adsorbed molecules released per cosmic-ray impact
      ELSE
         YIELD=0.0D0
      END IF
      FLUX=2.06D-3 ! Flux of iron nuclei cosmic rays (in cm^-2 s^-1)
      RATE(I)=FLUX*ZETA*2.4D-22*YIELD
      GOTO 20

!-----------------------------------------------------------------------

!     Photodesorption:

 9    IF(GAS_TEMPERATURE.LT.50.0D0) THEN
         YIELD=3.5D-3
      ELSE IF(GAS_TEMPERATURE.LT.85.0D0) THEN
         YIELD=4.0D-3
      ELSE IF(GAS_TEMPERATURE.LT.100.0D0) THEN
         YIELD=5.5D-3
      ELSE
         YIELD=7.5D-3
      END IF
      FLUX=1.71D8 ! Flux of FUV photons in the unattenuated Draine field (in photons cm^-2 s^-1)
      DO K=0,NRAYS-1 ! Loop over all rays
         RATE(I)=RATE(I) + FLUX*FUV_SURFACE(K)*EXP(-(GAMMA(I)*AV(K)))*2.4D-22*YIELD
      END DO
      GOTO 20

!-----------------------------------------------------------------------

!     Thermal desorption:

!     Treatment of Hasegawa, Herbst & Leung (1992, ApJS, 82, 167, Equations 2 & 3)
 10   RATE(I)=SQRT(2.0D0*1.5D15*KB/(PI**2*AU)*ALPHA(I)/GAMMA(I))*EXP(-(ALPHA(I)/DUST_TEMPERATURE))
      GOTO 20

!-----------------------------------------------------------------------

!     Grain mantle reactions:

 11   RATE(I)=ALPHA(I)
      GOTO 20

!-----------------------------------------------------------------------

!     Check that the rate is physical (0<RATE(I)<1) and produce an error
!     message if not. Impose a lower cut-off on all rate coefficients to
!     prevent the problem becoming too stiff. Rates lower than 1E-99 are
!     set to zero. Grain-surface reactions and desorption mechanisms are
!     allowed rates greater than 1.
 20   IF(RATE(I).LT.0.0D0) THEN
        WRITE(6,*) 'ERROR! Negative rate coefficient for reaction #',TRIM(NUM2STR(I))
        WRITE(6,*) 'Rate = ',TRIM(NUM2STR(RATE(I)))
        WRITE(6,*)
        STOP
      END IF
      IF(RATE(I).GT.1.0D0 .AND. REACTANT(I,1)(1:1).NE."G" .AND. REACTANT(I,1)(1:1).NE."#") THEN
        WRITE(10,*) 'WARNING! Rate coefficient is too large for reaction #',TRIM(NUM2STR(I))
        WRITE(10,*) 'Rate = ',TRIM(NUM2STR(RATE(I)))
        RATE(I)=1.0D0
      END IF
      IF(RATE(I).LT.1.0D-99) RATE(I)=0.0D0

!  End of loop over reactions
   END DO

   RETURN
END SUBROUTINE CALCULATE_REACTION_RATES
!=======================================================================
