!=======================================================================
!
!  Calculate the various heating rates for the current particle.
!
!-----------------------------------------------------------------------
SUBROUTINE CALCULATE_HEATING_RATES(NSPEC,NREAC,PARTICLE)

   USE HEALPIX_TYPES
   USE PARTICLE_MODULE
   USE GLOBAL_MODULE

   IMPLICIT NONE

   INTEGER(KIND=I4B),   INTENT(IN)    :: NSPEC,NREAC
   TYPE(PARTICLE_TYPE), INTENT(INOUT) :: PARTICLE

   REAL(KIND=DP) :: TH85_PHOTOELECTRIC_HEATING_RATE,BT94_PHOTOELECTRIC_HEATING_RATE,WD01_PHOTOELECTRIC_HEATING_RATE
   REAL(KIND=DP) :: CARBON_IONIZATION_HEATING_RATE,H2_FORMATION_HEATING_RATE,H2_DISSOCIATION_HEATING_RATE,H2_FUV_PUMPING_HEATING_RATE
   REAL(KIND=DP) :: COSMIC_RAY_HEATING_RATE,TURBULENT_HEATING_RATE,CHEMICAL_HEATING_RATE,GAS_GRAIN_HEATING_RATE,COULOMB_HEATING_RATE

   REAL(KIND=DP) :: GAS_DENSITY,GAS_TEMPERATURE,DUST_TEMPERATURE
   REAL(KIND=DP), ALLOCATABLE :: DENSITY(:),ABUNDANCE(:),REACTION_RATE(:)

!  Tielens & Hollenbach (1985) treatment of PE heating
   INTEGER(KIND=I4B) :: ITERATION
   REAL(KIND=DP) :: HABING_FIELD
   REAL(KIND=DP) :: X,XX,XK,XD,GAMMA,DELTA
   REAL(KIND=DP) :: DELTA_D,DELTA_UV,Y,HNU_D,HNU_H

!  Bakes & Tielens (1994) treatment of PE heating/cooling
   REAL(KIND=DP) :: EPSILON,ALPHA,BETA
   REAL(KIND=DP) :: PAH_HEATING_RATE,PAH_COOLING_RATE

!  Wolfire et al. (2003) scaling factor for PAH collision rates
   REAL(KIND=DP) :: PHI_PAH

!  Weingartner & Draine (2001) treatment of PE heating
   REAL(KIND=DP) :: C0,C1,C2,C3,C4,C5,C6

!  H2* FUV pumping heating
   REAL(KIND=DP) :: NCRIT_H2

!  Rollig et al. (2006) treatment of H2 vibrational heating/cooling
   REAL(KIND=DP) :: H2_VIBRATIONAL_HEATING,H2_VIBRATIONAL_COOLING
   REAL(KIND=DP) :: DELTA_E_10,A_COEFF_10,C_COEFF_10,R_PHOTO_1
   REAL(KIND=DP) :: DELTA_E_EFF,R_PUMP_EFF,R_PHOTO_EFF,A_COEFF_EFF,C_COEFF_EFF

!  Turbulent heating
   REAL(KIND=DP) :: L_TURB

!  Gas-grain coupling heating/cooling
   REAL(KIND=DP) :: N_GRAIN,C_GRAIN,ACCOMMODATION

!  Coulomb heating
   REAL(KIND=DP) :: R,X_PRIME,ETA_H2_He,ETA_H_He,ETA,H_X

!  Store the particle physical properties
   GAS_DENSITY=PARTICLE%GAS_DENSITY
   GAS_TEMPERATURE=PARTICLE%GAS_TEMPERATURE
   DUST_TEMPERATURE=PARTICLE%DUST_TEMPERATURE

   ALLOCATE(DENSITY(1:NSPEC))
   ALLOCATE(ABUNDANCE(1:NSPEC))
   ALLOCATE(REACTION_RATE(1:NREAC))
   DENSITY=PARTICLE%ABUNDANCE*PARTICLE%GAS_DENSITY
   ABUNDANCE=PARTICLE%ABUNDANCE
   REACTION_RATE=PARTICLE%RATE

!  Convert the FUV flux (in Draine units) to the Habing equivalent
   HABING_FIELD=1.71D0*PARTICLE%FUV_FLUX

!-----------------------------------------------------------------------
!  Grain photoelectric heating (large grains only; r ~ 100 Å)
!
!  Use the treatment of Tielens & Hollenbach (1985, ApJ, 291, 722)
!  which in turn follows de Jong (1977, 1980)
!
!  The charge of a dust grain can be found by equating the rate of
!  photo-ejection of electrons from the dust grain to the rate of
!  recombination of electrons with the dust grain (Spitzer)
!
!  The various parameter values are taken from Table 2 of the paper
!-----------------------------------------------------------------------

   DELTA_D=1.0D0
   DELTA_UV=1.8D0
   Y=0.1D0
   HNU_D=6.0D0
   HNU_H=13.6D0

   XK=KB*GAS_TEMPERATURE/(HNU_H*eV)
   XD=HNU_D/HNU_H
   GAMMA=2.9D-4*Y*SQRT(GAS_TEMPERATURE)*HABING_FIELD/DENSITY(nelect)
   DELTA=XK-XD+GAMMA

!  Iterate to determine X by finding the zero of the function F
   X=0.5D0
   DO ITERATION=1,100
      XX=X-(F(X,DELTA,GAMMA)/FF(X,DELTA))
      IF(ABS(XX-X).LT.1.0D-2) EXIT
      X=XX
   END DO
   X=XX

   IF(ITERATION.GE.100) THEN
      WRITE(10,*)'WARNING! Grain parameter X not found in PE heating'
      WRITE(10,*)'Using final value from interation loop: X =',X
   END IF

   TH85_PHOTOELECTRIC_HEATING_RATE=2.7D-25*DELTA_UV*DELTA_D*GAS_DENSITY*Y*HABING_FIELD &
                                 & *(((1.0D0-X)**2)/X + XK*((X**2)-1.0D0)/(X**2))

!  Assume the PE heating rate scales linearly with metallicity
   TH85_PHOTOELECTRIC_HEATING_RATE=TH85_PHOTOELECTRIC_HEATING_RATE*METALLICITY

!-----------------------------------------------------------------------
!  Grain + PAH photoelectric heating (MRN size distribution; r = 3-100 Å)
!
!  Use the treatment of Bakes & Tielens (1994, ApJ, 427, 822) with the
!  modifications suggested by Wolfire et al. (2003, ApJ, 587, 278) to
!  account for the revised PAH abundance estimate from Spitzer data.
!
!  See also:
!  Wolfire et al. (1995, ApJ, 443, 152)
!  Le Page, Snow & Bierbaum (2001, ApJS, 132, 233)
!-----------------------------------------------------------------------

!  Adopt the PAH rate scaling factor of Wolfire et al. (2008, ApJ, 680, 384)
!  Setting this factor to 1.0 gives the standard Bakes & Tielens expression
   PHI_PAH=0.4D0

   ALPHA=0.944D0
   BETA=0.735D0/GAS_TEMPERATURE**0.068
   DELTA=HABING_FIELD*SQRT(GAS_TEMPERATURE)/(DENSITY(nelect)*PHI_PAH)
   EPSILON=4.87D-2/(1.0D0+4.0D-3*DELTA**0.73) + 3.65D-2*(GAS_TEMPERATURE/1.0D4)**0.7/(1.0D0+2.0D-4*DELTA)

   PAH_HEATING_RATE=1.30D-24*EPSILON*HABING_FIELD*GAS_DENSITY
   PAH_COOLING_RATE=4.65D-30*GAS_TEMPERATURE**ALPHA*(DELTA**BETA)*DENSITY(nelect)*PHI_PAH*GAS_DENSITY

   BT94_PHOTOELECTRIC_HEATING_RATE=PAH_HEATING_RATE - PAH_COOLING_RATE

!  Assume the PE heating rate scales linearly with metallicity
   BT94_PHOTOELECTRIC_HEATING_RATE=BT94_PHOTOELECTRIC_HEATING_RATE*METALLICITY

!-----------------------------------------------------------------------
!  Grain + PAH photoelectric heating (with graphitic and silicate grains)

!  Weingartner & Draine (2001, ApJS, 134, 263)
!
!  Includes photoelectric heating due to PAHs, VSGs and larger grains
!  Assumes a gas-to-dust mass ratio of 100:1
!-----------------------------------------------------------------------

   C0=5.72D+0
   C1=3.45D-2
   C2=7.08D-3
   C3=1.98D-2
   C4=4.95D-1
   C5=6.92D-1
   C6=5.20D-1

   WD01_PHOTOELECTRIC_HEATING_RATE=1.0D-26*(HABING_FIELD*GAS_DENSITY)*(C0+C1*GAS_TEMPERATURE**C4) &
           & /(1.0D0+C2*(HABING_FIELD*SQRT(GAS_TEMPERATURE)/DENSITY(nelect))**C5  &
           & *(1.0D0+C3*(HABING_FIELD*SQRT(GAS_TEMPERATURE)/DENSITY(nelect))**C6))

!  Assume the PE heating rate scales linearly with metallicity
   WD01_PHOTOELECTRIC_HEATING_RATE=WD01_PHOTOELECTRIC_HEATING_RATE*METALLICITY

!-----------------------------------------------------------------------
!  H2 formation heating
!
!  Assume that 1.5 eV is liberated as heat during H2 formation
!
!  See: Hollenbach & Tielens (Review of Modern Physics, 1999, 71, 173)
!  Use the H2 formation rate determined by subroutine H2_FORMATION_RATE
!  and stored as REACTION_RATE(nRGR) (cm^3 s^-1)
!-----------------------------------------------------------------------

   H2_FORMATION_HEATING_RATE=(1.5*eV)*REACTION_RATE(nRGR)*DENSITY(nH)*GAS_DENSITY

!-----------------------------------------------------------------------
!  H2 FUV pumping heating
!
!  On average, 2.2 eV released per vibrationally excited H2* molecule
!
!  Use the treatment of Hollenbach & McKee (1979, ApJ)

!  Use the H2 photodissociation rate determined by the subroutine
!  CALCULATE_REACTION_RATES and stored as REACTION_RATE(nRH2) (s^-1)
!
!  Use the H2 critical density expression from Hollenbach & McKee (1979)
!  NOTE: The equation for the collisional de-excitation rate coefficient
!  for the v=2-1 and v=1-0 transitions by collisions with H2 was wrongly
!  stated in the Hollenbach & McKee (1979) paper (equation 6.29), but is
!  corrected in Hollenbach & McKee (1989, equation 2.8) and used below.
!-----------------------------------------------------------------------

   NCRIT_H2=1.0D6/SQRT(GAS_TEMPERATURE)/(1.6D0*ABUNDANCE(nH)*EXP(-((400.0D0/GAS_TEMPERATURE)**2)) &
                                     & + 1.4D0*ABUNDANCE(nH2)*EXP(-(18100.0D0/(GAS_TEMPERATURE+1200.0D0))))

   H2_FUV_PUMPING_HEATING_RATE=(2.2*eV)*9.0D0*REACTION_RATE(nRH2)*DENSITY(nH2)/(1.0D0+NCRIT_H2/GAS_DENSITY)

!  If vibrationally excited H2 (H2*) is included in the chemical network,
!  then use the treatment of Tielens & Hollenbach (1985, ApJ, 291, 722)
   IF(nH2v.NE.0) THEN
      H2_FUV_PUMPING_HEATING_RATE=(DENSITY(nH)*1.0D-12*SQRT(GAS_TEMPERATURE)*EXP(-1000.0D0/GAS_TEMPERATURE) &
                               & +DENSITY(nH2)*1.4D-12*SQRT(GAS_TEMPERATURE)*EXP(-18100.0D0/(GAS_TEMPERATURE+1200.0D0))) &
                                 *(2.6*eV)*DENSITY(nH2v)
   END IF

!-----------------------------------------------------------------------
!  H2 vibrational heating/cooling
!
!  Treat the vibrationally excited levels of H2 as a single pseudo level
!  with effective rates of spontaneous emission, collisional excitation,
!  FUV pumping and photodissociation that describe the behaviour of all
!  the vibrational levels combined.
!
!  Use the treatment of Rollig et al. (2006, A&A, 451, 917)
!-----------------------------------------------------------------------

   DELTA_E_10=6587.0 ! Energy gap (K) between the v=1 and v=0 levels of H2
   A_COEFF_10=8.6D-7 ! Einstein A-coefficient (s^-1) for emission from the v=1 to v=0 level of H2
   C_COEFF_10=5.4D-13*SQRT(GAS_TEMPERATURE) ! Collisional rate coefficient (cm^3 s^-1) for v=0 to v=1
   R_PHOTO_1=REACTION_RATE(nRH2) ! Photodissociation rate (s^-1) for the v=1 level of H2

   DELTA_E_EFF=23500.0 ! Characteristic vibrational level energy (K)
   A_COEFF_EFF=1.9D-6  ! Effective Einstein A-coefficient (s^-1)
   C_COEFF_EFF=C_COEFF_10 ! Effective collisional rate coefficient (cm^3 s^-1)
!!$   R_PUMP_EFF=2.9D-10 ! Effective FUV pumping rate (s^-1)
!!$   R_PHOTO_EFF=4.7D-10 ! Effective photodissociation rate (s^-1)
   R_PUMP_EFF=11.2*REACTION_RATE(nRH2) ! Effective vibrational pumping rate (s^-1)
   R_PHOTO_EFF=18.0*REACTION_RATE(nRH2) ! Effective photodissociation rate (s^-1)

   H2_VIBRATIONAL_COOLING=KB*DELTA_E_10*C_COEFF_10*GAS_DENSITY*EXP(-DELTA_E_10/GAS_TEMPERATURE)*DENSITY(nH2) &
                       & *(A_COEFF_10+R_PHOTO_1)/(C_COEFF_10*GAS_DENSITY+A_COEFF_10+R_PHOTO_1)

   H2_VIBRATIONAL_HEATING=DENSITY(nH2)*(R_PUMP_EFF*KB*DELTA_E_EFF) &
                       & /(1.0D0+(A_COEFF_EFF+R_PHOTO_EFF)/(C_COEFF_EFF*GAS_DENSITY))

!   H2_FUV_PUMPING_HEATING_RATE=H2_VIBRATIONAL_HEATING-H2_VIBRATIONAL_COOLING

!-----------------------------------------------------------------------
!  H2 photodissociation heating
!
!  On average, 0.4 eV of kinetic energy per photodissociated molecule
!
!  Use the H2 photodissociation rate determined by the subroutine
!  CALCULATE_REACTION_RATES and stored as REACTION_RATE(nRH2) (s^-1)
!-----------------------------------------------------------------------

   H2_DISSOCIATION_HEATING_RATE=(0.4*eV)*REACTION_RATE(nRH2)*DENSITY(nH2)

!-----------------------------------------------------------------------
!  Carbon photoionization heating
!
!  On average, 1 eV of kinetic energy released per carbon ionization
!  Use the carbon photoionization rate determined by the subroutine
!  CALCULATE_REACTION_RATES and stored as REACTION_RATE(nRCI) (s^-1)
!-----------------------------------------------------------------------

   CARBON_IONIZATION_HEATING_RATE=(1.0*eV)*REACTION_RATE(nRCI)*DENSITY(nC)

!-----------------------------------------------------------------------
!  Cosmic-ray ionization heating
!
!  20.0 eV of kinetic energy deposited per H2 ionization,
!  based on the estimate of Goldsmith (2001, ApJ, 557,736)
!
!  See also:
!  Clavel et al. (1978, A&A, 65, 435)
!  Tielens & Hollenbach (1985, ApJ, 291, 722)
!  Shull & van Steenberg (1985, ApJ, 298, 268)
!  Kamp & van Zadelhoff (2001)
!-----------------------------------------------------------------------

   COSMIC_RAY_HEATING_RATE=(20.0*eV)*(1.3D-17*ZETA)*DENSITY(nH2)

!-----------------------------------------------------------------------
!  Supersonic turbulent decay heating
!
!  Most relevant for the inner parsecs of galaxies (Black)
!  Black, in Interstellar Processes, 1987, p731
!  See also: Rodriguez-Fernandez et al., 2001, A&A, 365, 174
!
!  V_TURB = turbulent velocity (km/s); Galactic center ~ 15 km/s
!  L_TURB = turbulent scale length (pc); typically 5 pc
!-----------------------------------------------------------------------

   L_TURB=5.0D0
   TURBULENT_HEATING_RATE=3.5D-28*((V_TURB/1.0D5)**3)*(1.0D0/L_TURB)*GAS_DENSITY

!-----------------------------------------------------------------------
!  Exothermic chemical reaction heating
!
!  See:
!  Clavel et al. (1978, A&A,  65, 435)
!  Meijerink & Spaans (2005, A&A, 436, 397)
!  Glassgold & Langer (1973, ApJ, 179, L147)
!
!  Recombination reactions:
!     H2+ (10.9 eV); H3+ (9.23+4.76 eV); H3O+ (1.16+5.63+6.27 eV); HCO+ (7.51 eV)
!
!  Ion-neutral reactions:
!     H2+ + H (0.94 eV); He+ + H2 (6.51 eV); He+ + CO (2.22 eV)
!
!  For each reaction, the heating rate is given by: n(1) * n(2) * K * E
!  where n(1) and n(2) are the number densities, K the rate coefficient
!  (cm^3 s^-1), and E the energy released (erg).
!-----------------------------------------------------------------------

   CHEMICAL_HEATING_RATE = &
    & + DENSITY(nH2x) *DENSITY(nelect)*(REACTION_RATE(nR_H2x_1) *(10.9*eV)) & ! H2+ + e-
    & + DENSITY(nH2x) *DENSITY(nH)    *(REACTION_RATE(nR_H2x_2) *(0.94*eV)) & ! H2+ + H
    & + DENSITY(nH3x) *DENSITY(nelect)*(REACTION_RATE(nR_H3x_1) *(9.23*eV)+REACTION_RATE(nR_H3x_2)*(4.76*eV)) & ! H3+ + e-
    & + DENSITY(nH3Ox)*DENSITY(nelect)*(REACTION_RATE(nR_H3Ox_1)*(1.16*eV)+REACTION_RATE(nR_H3Ox_2)*(5.63*eV)+REACTION_RATE(nR_H3Ox_3)*(6.27*eV)) & ! H3O+ + e-
    & + DENSITY(nHCOx)*DENSITY(nelect)*(REACTION_RATE(nR_HCOx_1)*(7.51*eV)) & ! HCO+ + e-
    & + DENSITY(nHex) *DENSITY(nH2)   *(REACTION_RATE(nR_Hex_1) *(6.51*eV)+REACTION_RATE(nR_Hex_2)*(6.51*eV)) & ! He+ + H2
    & + DENSITY(nHex) *DENSITY(nCO)   *(REACTION_RATE(nR_Hex_3) *(2.22*eV)+REACTION_RATE(nR_Hex_4)*(2.22*eV))   ! He+ + CO

!-----------------------------------------------------------------------
!  Gas-grain collisional heating
!
!  Use the treatment of Burke & Hollenbach (1983, ApJ, 265, 223)
!
!  Other relevant references:
!  Hollenbach & McKee (1979, ApJS, 41, 555)
!  Tielens & Hollenbach (1985, ApJ, 291, 722)
!  Goldsmith (2001, ApJ, 557, 736)
!  Young et al. (2004, ApJ, 614, 252)
!
!  This process is insignificant for the energy balance of the dust,
!  but can influence the gas temperature. If the dust temperature is
!  lower than the gas temperature, this becomes a cooling mechanism.
!-----------------------------------------------------------------------

   N_GRAIN=PARTICLE%DUST_DENSITY
   C_GRAIN=PI*GRAIN_RADIUS**2

!!$!  Accommodation fitting formula of Groenewegen (1994, A&A, 290, 531)
!!$   ACCOMMODATION=0.35D0*EXP(-SQRT((DUST_TEMPERATURE+GAS_TEMPERATURE)/5.0D2))+0.1D0

!  Accommodation coefficient of Burke & Hollenbach (1983, ApJ, 265, 223)
   ACCOMMODATION=0.37D0*(1.0D0-0.8D0*EXP(-75.0D0/GAS_TEMPERATURE))

   GAS_GRAIN_HEATING_RATE=N_GRAIN*C_GRAIN*GAS_DENSITY*SQRT(8*KB*GAS_TEMPERATURE/(PI*MH)) &
                       & *ACCOMMODATION*(2*KB*DUST_TEMPERATURE-2*KB*GAS_TEMPERATURE)

!-----------------------------------------------------------------------
!  Coulomb heating
!
!  Use the treatment of Meijerink & Spaans (2005, A&A, 436, 397)
!
!  Other relevant references:
!  Dalgarno et al. (1999, ApJS, 125, 237)
!
!  This is an X-ray heating mechanism. When X-rays are absorbed, fast
!  electrons are produced. These fast electrons lose part of their
!  energy through Coulomb interactions with thermal electrons.
!-----------------------------------------------------------------------

   R=DENSITY(nH2)/DENSITY(nH) ! n(H2)/n(H) ratio
   X_PRIME=1.83D0*ABUNDANCE(nelect)/(1.0D0+0.83D0*ABUNDANCE(nelect)) ! Correction to the electron abundance for a pure H2-He mixture

   ETA_H2_He=1.0D0+(0.055D0-1.0D0)/(1.0D0+2.17D0*X_PRIME**0.366) ! Heating efficiency for a pure H2-He mixture
   ETA_H_He =1.0D0+(0.117D0-1.0D0)/(1.0D0+7.95D0*ABUNDANCE(nelect)**0.678) ! Heating efficiency for a pure H-He mixture
   ETA=(10.0D0*R*ETA_H2_He+ETA_H_He)/(10.0D0*R+1.0D0) ! Total heating efficiency for mixed atomic and molecular gas
   H_X=PARTICLE%XRAY_ENERGY_DEPOSITION_RATE ! X-ray energy deposition rate per hydrogen nucleus (erg s^-1)

   COULOMB_HEATING_RATE=ETA*GAS_DENSITY*H_X

!-----------------------------------------------------------------------

!  Store the various heating rates
!   PARTICLE%HEATING_RATE(1)=TH85_PHOTOELECTRIC_HEATING_RATE
!   PARTICLE%HEATING_RATE(1)=WD01_PHOTOELECTRIC_HEATING_RATE
   PARTICLE%HEATING_RATE(1)=BT94_PHOTOELECTRIC_HEATING_RATE
   PARTICLE%HEATING_RATE(2)=H2_FORMATION_HEATING_RATE
   PARTICLE%HEATING_RATE(3)=H2_FUV_PUMPING_HEATING_RATE
   PARTICLE%HEATING_RATE(4)=H2_DISSOCIATION_HEATING_RATE
   PARTICLE%HEATING_RATE(5)=CARBON_IONIZATION_HEATING_RATE
   PARTICLE%HEATING_RATE(6)=COSMIC_RAY_HEATING_RATE
   PARTICLE%HEATING_RATE(7)=TURBULENT_HEATING_RATE
   PARTICLE%HEATING_RATE(8)=CHEMICAL_HEATING_RATE
   PARTICLE%HEATING_RATE(9)=GAS_GRAIN_HEATING_RATE
   PARTICLE%HEATING_RATE(10)=COULOMB_HEATING_RATE

   RETURN

CONTAINS ! Dust photoelectric heating functions

!=======================================================================
!  X is the grain charge parameter and is the solution to F(X)=0
!-----------------------------------------------------------------------
   FUNCTION F(X,DELTA,GAMMA)
      USE HEALPIX_TYPES
      IMPLICIT NONE
      REAL(KIND=DP) :: F
      REAL(KIND=DP), INTENT(IN) :: X,DELTA,GAMMA
      F=(X**3)+DELTA*(X**2)-GAMMA
   END FUNCTION F

!=======================================================================
!  FF(X) is the derivative of F(X) with respect to X
!-----------------------------------------------------------------------
   FUNCTION FF(X,DELTA)
      USE HEALPIX_TYPES
      IMPLICIT NONE
      REAL(KIND=DP) :: FF
      REAL(KIND=DP), INTENT(IN) :: X,DELTA
      FF=3*(X**2)+DELTA*(2*X)
   END FUNCTION FF

END SUBROUTINE CALCULATE_HEATING_RATES
!=======================================================================
