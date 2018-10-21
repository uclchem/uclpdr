!=======================================================================
MODULE POPULATION_MODULE

   USE HEALPIX_TYPES
   IMPLICIT NONE

!  Specify the properties that define the level populations and emission of a species
   TYPE POPULATION_TYPE

      CHARACTER(LEN=10) :: NAME ! Species name
      INTEGER(KIND=I4B) :: NLEVEL ! Number of levels in the system

      REAL(KIND=DP) :: DENSITY ! Total number density of the species (cm^-3)
      REAL(KIND=DP) :: LINEWIDTH ! Doppler line width of the species (cm s^-1)

      REAL(KIND=DP), POINTER :: POPULATION(:) ! Population density (cm^-3) of each level
      REAL(KIND=DP), POINTER :: PREVIOUS_POPULATION(:) ! Population density calculated at the previous iteration step
      REAL(KIND=DP), POINTER :: EMISSIVITY(:,:) ! Local emissivity (erg cm^-3 s^-1) of each transition
      REAL(KIND=DP), POINTER :: OPACITY(:,:,:) ! Optical depth of each transition along each HEALPix ray to the PDR surface (or simulation boundary)

#ifdef USE_ALI
      REAL(KIND=DP), POINTER :: LAMBDA(:,:) ! Lambda operator value for each transition
#endif

      LOGICAL :: CONVERGED ! Flag indicating whether the level populations have converged

   END TYPE POPULATION_TYPE

END MODULE POPULATION_MODULE
!=======================================================================

!=======================================================================
MODULE PARTICLE_MODULE

   USE HEALPIX_TYPES
   USE POPULATION_MODULE
   IMPLICIT NONE

!  Specify the physical properties that define each particle in the model cloud
   TYPE PARTICLE_TYPE

      INTEGER(KIND=I4B) :: TYPE ! Particle type (ionized=1, PDR=2, dark cloud=3)

      REAL(KIND=DP) :: COORDINATES(1:3) ! Particle coordinates (x,y,z; pc)
      REAL(KIND=DP) :: GAS_DENSITY,DUST_DENSITY ! Gas and dust number densities (cm^-3)
      REAL(KIND=DP) :: GAS_TEMPERATURE,DUST_TEMPERATURE ! Gas and dust temperatures (K)
      REAL(KIND=DP) :: FUV_FLUX,XRAY_FLUX ! Local FUV and X-ray fluxes (chi, Draines; F_X, erg cm^-2 s^-1)
      REAL(KIND=DP) :: XRAY_ENERGY_DEPOSITION_RATE ! X-ray photon energy deposition rate (H_X, erg s^-1)

      REAL(KIND=DP), POINTER :: FUV_WAVELENGTHS(:),FUV_FLUXES(:) ! Local wavelength-dependent FUV fluxes (erg cm^-2 s^-1 Å^-1)
      REAL(KIND=DP), POINTER :: XRAY_ENERGIES(:),XRAY_FLUXES(:) ! Local energy-dependent X-ray fluxes (erg cm^-2 s^-1 keV^-1)

      INTEGER(KIND=I4B), POINTER :: EVALUATION_POINT(:,:) ! List of evaluation points along each HEALPix ray to the PDR surface (or simulation boundary)
      INTEGER(KIND=I4B), POINTER :: NPOINT(:) ! Number of evaluation points along each ray

      REAL(KIND=DP), POINTER :: TOTAL_COLUMN(:),AV(:) ! Total hydrogen column density (N_H, cm^-2) and visual extinction (A_V, mag) along each HEALPix ray to the PDR surface (or simulation boundary)
      REAL(KIND=DP), POINTER :: FUV_SURFACE(:),XRAY_SURFACE(:) ! Unattenuated FUV and X-ray fluxes (chi_0, Draines; F_X, erg cm^-2 s^-1) at the PDR surface intersected by each HEALPix ray
      REAL(KIND=DP), POINTER :: RATE(:),ABUNDANCE(:) ! Rate coefficient of each reaction and fractional abundance of each species (with respect to the total hydrogen abundance)
      REAL(KIND=DP), POINTER :: COLUMN_DENSITY(:,:) ! Column density (cm^-2) of each species along each HEALPix ray to the PDR surface (or simulation boundary)
      REAL(KIND=DP), POINTER :: COOLING_RATE(:),HEATING_RATE(:) ! Cooling and heating rate (erg cm^-3 s^-1) of each mechanism considered

      TYPE(POPULATION_TYPE), POINTER :: COOLANT(:) ! Coolant species whose level populations and emission properties will be calculated

      REAL(KIND=DP), POINTER :: PREVIOUS_ABUNDANCE(:) ! Fractional abundance calculated at the previous iteration step
      LOGICAL       :: CHEMISTRY_CONVERGED ! Flag indicating whether the chemical abundances have converged

      REAL(KIND=DP) :: PREVIOUS_DIFFERENCE ! Heating-cooling imbalance at the previous iteration step
      REAL(KIND=DP) :: PREVIOUS_TEMPERATURE ! Gas temperature calculated at the previous iteration step
      REAL(KIND=DP) :: TLOW,THIGH ! Temperature brackets determined by the thermal balance iteration routine
      LOGICAL       :: BINARY_CHOP_SEARCH ! Flag indicating whether the binary chop search method is to be used
      LOGICAL       :: BRACKET_EXPANDED ! Flag indicating whether the binary chop search bracket has been expanded
      LOGICAL       :: TEMPERATURE_CONVERGED ! Flag indicating whether the gas temperature has converged

   END TYPE PARTICLE_TYPE

END MODULE PARTICLE_MODULE
!=======================================================================

!=======================================================================
MODULE COOLANT_MODULE

   USE HEALPIX_TYPES
   IMPLICIT NONE

!  Specify the properties that define each coolant species
   TYPE COOLANT_TYPE

      CHARACTER(LEN=10) :: NAME ! Coolant species name
      CHARACTER(LEN=256):: FILENAME ! Name of the coolant data file

      INTEGER(KIND=I4B) :: INDEX  ! Index number of the coolant species
      INTEGER(KIND=I4B) :: NLEVEL ! Number of levels in the system
      INTEGER(KIND=I4B) :: NTEMP  ! Number of temperatures for which collisional rates are available

      REAL(KIND=DP) :: MOLECULAR_MASS ! Molecular mass of the coolant species

      REAL(KIND=DP), POINTER :: ENERGY(:),WEIGHT(:) ! Energy (K) and statistical weight of each level
      REAL(KIND=DP), POINTER :: A_COEFF(:,:),B_COEFF(:,:) ! Einstein A and B coefficients for each transition between levels
      REAL(KIND=DP), POINTER :: FREQUENCY(:,:) ! Frequency (Hz) of each transition between levels
      REAL(KIND=DP), POINTER :: TEMPERATURE(:,:) ! Temperatures (K) at which collisional rates are available for each collision partner
      REAL(KIND=DP), POINTER :: C_COEFF(:,:,:,:) ! Collisional rate coefficient (cm^3 s^-1) for each transition at each specified temperature for each collision partner

      LOGICAL :: CONVERGED ! Flag indicating whether the level populations of all particles have converged

   END TYPE COOLANT_TYPE

END MODULE COOLANT_MODULE
!=======================================================================

!=======================================================================
MODULE CROSS_SECTION_MODULE

   USE HEALPIX_TYPES
   IMPLICIT NONE

!  Specify the properties that define each photoreaction cross section
   TYPE CROSS_SECTION_TYPE

      CHARACTER(LEN=10) :: SPECIES  ! Name of the species
      CHARACTER(LEN=256):: FILENAME ! Name of the cross-section data file

      INTEGER(KIND=I4B) :: INDEX ! Index number of the relevant photoreaction
      INTEGER(KIND=I4B) :: NLINE ! Number of discrete line absorption cross sections
      INTEGER(KIND=I4B) :: NCONT ! Number of wavelengths for which continuous cross section values are available

      REAL(KIND=DP) :: SCALING_FACTOR ! Scaling factor to account for the probability of dissociation into multiple channels

      REAL(KIND=DP), POINTER :: LINE_WAVELENGTH(:),LINE_CROSS_SECTION(:) ! Wavelength (Å) and cross section (cm^2 Å) of each line absorption
      REAL(KIND=DP), POINTER :: CONTINUUM_WAVELENGTH(:),CONTINUUM_CROSS_SECTION(:) ! Wavelength (Å) and cross section (cm^2) values for continuous absorption

   END TYPE CROSS_SECTION_TYPE

END MODULE CROSS_SECTION_MODULE
!=======================================================================

!=======================================================================
MODULE MAIN_MODULE

   USE HEALPIX_TYPES
   USE PARTICLE_MODULE
   USE COOLANT_MODULE
   USE CROSS_SECTION_MODULE

   IMPLICIT NONE

   INTEGER(KIND=I4B) :: NPART ! Number of particles in the model cloud
   INTEGER(KIND=I4B) :: NRAYS ! Number of rays in the HEALPix structure
   INTEGER(KIND=I4B) :: NSPEC ! Number of species in the chemical network
   INTEGER(KIND=I4B) :: NREAC ! Number of reactions in the chemical network
   INTEGER(KIND=I4B) :: NXSEC ! Number of photoreactions with cross sections
   INTEGER(KIND=I4B) :: NCOOL ! Number of species contributing to the cooling
   INTEGER(KIND=I4B) :: NHEAT ! Number of mechanisms contributing to the heating

   INTEGER(KIND=I4B) :: CHEMISTRY_ITERATIONS ! Number of chemical abundance iterations
   INTEGER(KIND=I4B) :: POPULATION_ITERATIONS ! Number of level population iterations
   INTEGER(KIND=I4B) :: TEMPERATURE_ITERATIONS ! Number of thermal balance iterations

   REAL(KIND=DP) :: TMIN,TMAX,TDIFF,FCRIT ! Thermal balance convergence parameters
   REAL(KIND=DP) :: GAS_TEMPERATURE_GUESS,DUST_TEMPERATURE_GUESS ! Initial guesses for the gas and dust temperature
   REAL(KIND=DP) :: ABUNDANCE_LIMIT,CHEMISTRY_CONVERGENCE_CRITERION ! Chemistry convergence parameters
   REAL(KIND=DP) :: POPULATION_LIMIT,POPULATION_CONVERGENCE_CRITERION ! Population convergence parameters
   REAL(KIND=DP) :: POPULATION_ITERATION_LIMIT,POPULATION_PERCENTAGE_LIMIT ! Population convergence parameters

   TYPE(PARTICLE_TYPE), ALLOCATABLE :: PARTICLE(:) ! Array of particles in the model cloud
   TYPE(COOLANT_TYPE),  ALLOCATABLE :: COOLANT(:)  ! Array of coolants in the model cloud
   TYPE(CROSS_SECTION_TYPE), ALLOCATABLE :: CROSS_SECTION(:) ! Array of photoreaction cross sections

END MODULE MAIN_MODULE
!=======================================================================

!=======================================================================
MODULE HEALPIX_MODULE

   USE HEALPIX_TYPES
   IMPLICIT NONE

   INTEGER(KIND=I4B) :: HEALPIX_LEVEL ! HEALPix level of refinement
   REAL(KIND=DP) :: THETA_CRITICAL ! Evaluation point search angle

END MODULE HEALPIX_MODULE
!=======================================================================

!=======================================================================
MODULE CHEMISTRY_MODULE

   USE ISO_C_BINDING
   USE HEALPIX_TYPES
   IMPLICIT NONE

   CHARACTER(LEN=10), ALLOCATABLE :: SPECIES(:) ! Names of the species in the chemical network
   REAL(KIND=DP),     ALLOCATABLE :: INITIAL_ABUNDANCE(:),MOLECULAR_MASS(:) ! Initial abundance and molecular mass of each species

   CHARACTER(LEN=10), ALLOCATABLE :: REACTANT(:,:),PRODUCT(:,:) ! Names of the reactants and products involved in each chemical reaction
   REAL(KIND=DP),     ALLOCATABLE :: ALPHA(:),BETA(:),GAMMA(:),RTMIN(:),RTMAX(:) ! Arrhenius equation parameters for each chemical reaction and the minimum/maximum temperatures within which they are valid
   INTEGER(KIND=I4B), ALLOCATABLE :: DUPLICATE(:) ! Duplicate reaction indices (0: new reaction; 1: duplicate of the previous reaction; 2: second duplicate of a previous reaction; etc.)

   REAL(KIND=DP),BIND(C,name="chemistry_module_mp_relative_abundance_tolerance_") :: RELATIVE_ABUNDANCE_TOLERANCE
   REAL(KIND=DP),BIND(C,name="chemistry_module_mp_absolute_abundance_tolerance_") :: ABSOLUTE_ABUNDANCE_TOLERANCE ! Relative (RTOL) and absolute (ATOL) error tolerances for the calculated abundances

END MODULE CHEMISTRY_MODULE
!=======================================================================

!=======================================================================
MODULE GLOBAL_MODULE

   USE HEALPIX_TYPES
   USE ISO_C_BINDING
   IMPLICIT NONE

!  Species and reaction indices (assigned within the subroutines READ_SPECIES and READ_REACTIONS)
   INTEGER(KIND=I4B),BIND(C,name="global_module_mp_nelect_") :: nelect
   INTEGER(KIND=I4B) :: nH,nHx,nD,nDx,nH2,nH2x,nH3x,nHD,nH2Dx,nHe,nHex,nC,nCx,nN,nNx,nO,nOx, &
                      & nF,nFx,nNa,nNax,nMg,nMgx,nSi,nSix,nS,nSx,nCl,nClx,nCa,nCax,nCaxx,nFe,nFex,  &
                      & nCH,nCHx,nCH2,nCH2x,nCH3x,nOH,nH2O,nH2Ox,nH3Ox,nNH,nNH2,nNH3,nCO,nHCOx,nCS, &
                      & nH2v,nPAH,nPAHx,nPAHm
   INTEGER(KIND=I4B) :: nRGR,nRH2,nRHD,nRCO,nRCI,nRSI
   INTEGER(KIND=I4B) :: nR_H2x_1,nR_H2x_2,nR_H3x_1,nR_H3x_2, &
                      & nR_Hex_1,nR_Hex_2,nR_Hex_3,nR_Hex_4, &
                      & nR_H3Ox_1,nR_H3Ox_2,nR_H3Ox_3,nR_HCOx_1

!  Specify various global properties used throughout the code
   REAL(KIND=DP),BIND(C,name="global_module_mp_start_time_") :: START_TIME
   REAL(KIND=DP),BIND(C,name="global_module_mp_end_time_") :: END_TIME ! Start and end times for the chemical evolution (yr)
   REAL(KIND=DP) :: ZETA ! Cosmic-ray ionization rate (s^-1)
   REAL(KIND=DP) :: T_CMB ! Cosmic microwave background temperature (K)
   REAL(KIND=DP) :: V_TURB ! Microturbulent velocity (cm s^-1)
   REAL(KIND=DP) :: METALLICITY ! Metallicity (Z_solar)
   REAL(KIND=DP) :: AV_FAC ! A_V/N_H conversion factor (mag cm^2)
   REAL(KIND=DP) :: UV_FAC ! Attenuation scaling factor between visible and FUV wavelengths
   REAL(KIND=DP) :: SIGMA  ! Photoelectric absorption cross section per hydrogen nucleus (cm^2)
   REAL(KIND=DP) :: OMEGA ! Dust grain albedo
   REAL(KIND=DP) :: GRAIN_RADIUS ! Dust grain radius (cm)

END MODULE GLOBAL_MODULE
!=======================================================================

!=======================================================================
MODULE SHIELDING_MODULE

   USE HEALPIX_TYPES
   IMPLICIT NONE

   LOGICAL :: START=.TRUE.

!  H2 line shielding data from Lee et al. (1996, A&A, 311, 690, Table 10)
   INTEGER(KIND=I4B),             SAVE :: NUM_COL=105
   REAL(KIND=DP), DIMENSION(105), SAVE :: COL_GRID=(/ &
      & 0.000D+00,3.690D+11,3.715D+12,3.948D+13,1.233D+14, &
      & 2.536D+14,4.342D+14,6.653D+14,6.689D+14,9.075D+14, &
      & 1.234D+15,1.631D+15,2.105D+15,2.363D+15,2.899D+15, &
      & 3.207D+15,3.848D+15,4.636D+15,5.547D+15,6.604D+15, &
      & 7.855D+15,9.368D+15,1.122D+16,1.352D+16,1.643D+16, &
      & 2.017D+16,2.515D+16,3.190D+16,4.128D+16,5.439D+16, &
      & 7.315D+16,1.009D+17,1.432D+17,2.092D+17,3.123D+17, &
      & 4.738D+17,5.388D+17,8.935D+17,1.381D+18,2.164D+18, &
      & 3.330D+18,5.024D+18,7.404D+18,9.029D+18,1.316D+19, &
      & 1.813D+19,2.453D+19,3.248D+19,4.216D+19,5.370D+19, &
      & 6.722D+19,8.277D+19,9.894D+19,1.186D+20,1.404D+20, &
      & 1.644D+20,1.908D+20,2.197D+20,2.510D+20,2.849D+20, &
      & 3.214D+20,3.604D+20,4.019D+20,4.456D+20,4.915D+20, &
      & 5.393D+20,5.886D+20,6.392D+20,6.909D+20,7.433D+20, &
      & 7.965D+20,8.505D+20,9.056D+20,9.627D+20,1.011D+21, &
      & 1.068D+21,1.125D+21,1.185D+21,1.250D+21,1.327D+21, &
      & 1.428D+21,1.578D+21,1.851D+21,2.128D+21,2.298D+21, &
      & 2.389D+21,2.459D+21,2.519D+21,2.571D+21,2.618D+21, &
      & 2.707D+21,2.790D+21,2.887D+21,3.001D+21,3.139D+21, &
      & 3.303D+21,3.497D+21,3.722D+21,3.983D+21,4.283D+21, &
      & 4.644D+21,5.127D+21,5.945D+21,8.205D+21,1.015D+22/)
   REAL(KIND=DP), DIMENSION(105), SAVE :: SH2_GRID=(/ &
      & 1.000D+00,9.983D-01,9.853D-01,8.761D-01,7.199D-01, &
      & 5.728D-01,4.455D-01,3.431D-01,3.418D-01,2.732D-01, &
      & 2.110D-01,1.619D-01,1.236D-01,1.084D-01,8.447D-02, &
      & 7.410D-02,5.774D-02,4.416D-02,3.390D-02,2.625D-02, &
      & 2.048D-02,1.606D-02,1.264D-02,9.987D-03,7.937D-03, &
      & 6.343D-03,5.088D-03,4.089D-03,3.283D-03,2.640D-03, &
      & 2.130D-03,1.725D-03,1.397D-03,1.129D-03,9.097D-04, &
      & 7.340D-04,6.883D-04,5.377D-04,4.352D-04,3.475D-04, &
      & 2.771D-04,2.205D-04,1.753D-04,1.549D-04,1.210D-04, &
      & 9.666D-05,7.705D-05,6.148D-05,4.904D-05,3.909D-05, &
      & 3.112D-05,2.473D-05,1.997D-05,1.578D-05,1.244D-05, &
      & 9.769D-06,7.634D-06,5.932D-06,4.581D-06,3.515D-06, &
      & 2.679D-06,2.029D-06,1.527D-06,1.144D-06,8.523D-07, &
      & 6.332D-07,4.693D-07,3.475D-07,2.574D-07,1.907D-07, &
      & 1.413D-07,1.047D-07,7.739D-08,5.677D-08,4.386D-08, &
      & 3.227D-08,2.385D-08,1.750D-08,1.248D-08,8.389D-09, &
      & 5.026D-09,2.382D-09,6.259D-10,1.653D-10,7.399D-11, &
      & 4.824D-11,3.474D-11,2.633D-11,2.069D-11,1.663D-11, &
      & 1.099D-11,7.506D-12,4.825D-12,2.864D-12,1.534D-12, &
      & 7.324D-13,3.087D-13,1.135D-13,3.591D-14,9.689D-15, &
      & 2.045D-15,2.618D-16,8.918D-18,3.041D-21,1.739D-23/)
   REAL(KIND=DP), DIMENSION(105), SAVE :: SH2_DERIV

!  12CO line shielding data from van Dishoeck & Black (1988, ApJ, 334, 771, Table 5)
   INTEGER(KIND=I4B),             SAVE :: NUM_NCO=8,NUM_NH2=6
   REAL(KIND=DP), DIMENSION(8),   SAVE :: NCO_GRID=(/12.0D0,13.0D0,14.0D0,15.0D0,16.0D0,17.0D0,18.0D0,19.0D0/)
   REAL(KIND=DP), DIMENSION(6),   SAVE :: NH2_GRID=(/18.0D0,19.0D0,20.0D0,21.0D0,22.0D0,23.0D0/)
   REAL(KIND=DP), DIMENSION(8,6), SAVE :: SCO_GRID=RESHAPE((/ &
      &  0.000D+00,-1.408D-02,-1.099D-01,-4.400D-01,-1.154D+00,-1.888D+00,-2.760D+00,-4.001D+00, &
      & -8.539D-02,-1.015D-01,-2.104D-01,-5.608D-01,-1.272D+00,-1.973D+00,-2.818D+00,-4.055D+00, &
      & -1.451D-01,-1.612D-01,-2.708D-01,-6.273D-01,-1.355D+00,-2.057D+00,-2.902D+00,-4.122D+00, &
      & -4.559D-01,-4.666D-01,-5.432D-01,-8.665D-01,-1.602D+00,-2.303D+00,-3.146D+00,-4.421D+00, &
      & -1.303D+00,-1.312D+00,-1.367D+00,-1.676D+00,-2.305D+00,-3.034D+00,-3.758D+00,-5.077D+00, &
      & -3.883D+00,-3.888D+00,-3.936D+00,-4.197D+00,-4.739D+00,-5.165D+00,-5.441D+00,-6.446D+00/), (/8,6/))
   REAL(KIND=DP), DIMENSION(8,6), SAVE :: SCO_DERIV

!  Ratio of tau(λ)/tau(V) based on the extinction curve from Savage & Mathis (1979, ARA&A, 17, 73, Table 2)
   INTEGER(KIND=I4B),             SAVE :: NUM_LAMBDA=30
   REAL(KIND=DP), DIMENSION(30),  SAVE :: LAMBDA_GRID=(/ &
      &  910.0D0, 950.0D0,1000.0D0,1050.0D0,1110.0D0, &
      & 1180.0D0,1250.0D0,1390.0D0,1490.0D0,1600.0D0, &
      & 1700.0D0,1800.0D0,1900.0D0,2000.0D0,2100.0D0, &
      & 2190.0D0,2300.0D0,2400.0D0,2500.0D0,2740.0D0, &
      & 3440.0D0,4000.0D0,4400.0D0,5500.0D0,7000.0D0, &
      & 9000.0D0,12500.0D0,22000.0D0,34000.0D0,1.0D9/)
   REAL(KIND=DP), DIMENSION(30),  SAVE :: XLAMBDA_GRID=(/ &
      & 5.76D0,5.18D0,4.65D0,4.16D0,3.73D0, &
      & 3.40D0,3.11D0,2.74D0,2.63D0,2.62D0, &
      & 2.54D0,2.50D0,2.58D0,2.78D0,3.01D0, &
      & 3.12D0,2.86D0,2.58D0,2.35D0,2.00D0, &
      & 1.58D0,1.42D0,1.32D0,1.00D0,0.75D0, &
      & 0.48D0,0.28D0,0.12D0,0.05D0,0.00D0/)
   REAL(KIND=DP), DIMENSION(30),  SAVE :: XLAMBDA_DERIV

END MODULE SHIELDING_MODULE
!=======================================================================

!=======================================================================
MODULE FUNCTIONS_MODULE

   INTERFACE

      FUNCTION COUNT_LINES(FILENAME) RESULT(LINES)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         INTEGER(KIND=I4B) :: LINES
         CHARACTER(LEN=*), INTENT(IN) :: FILENAME
      END FUNCTION COUNT_LINES

      FUNCTION COUNT_SUBSTRING(STRING,SUBSTRING) RESULT(COUNT)
         IMPLICIT NONE
         INTEGER :: COUNT
         CHARACTER(LEN=*), INTENT(IN) :: STRING,SUBSTRING
      END FUNCTION COUNT_SUBSTRING

      FUNCTION ORTHO_PARA_RATIO(TEMPERATURE)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         REAL(KIND=DP) :: ORTHO_PARA_RATIO
         REAL(KIND=DP), INTENT(IN) :: TEMPERATURE
      END FUNCTION ORTHO_PARA_RATIO

      FUNCTION PHOTOREACTION_RATE(CHI,AV,CROSS_SECTION) RESULT(TOTAL_RATE)
         USE HEALPIX_TYPES
         USE CROSS_SECTION_MODULE
         IMPLICIT NONE
         REAL(KIND=DP),            INTENT(IN) :: CHI,AV
         TYPE(CROSS_SECTION_TYPE), INTENT(IN) :: CROSS_SECTION
         REAL(KIND=DP) :: TOTAL_RATE
      END FUNCTION PHOTOREACTION_RATE

      FUNCTION H2_FORMATION_RATE(GAS_TEMPERATURE,GRAIN_TEMPERATURE) RESULT(RATE)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         REAL(KIND=DP) :: RATE
         REAL(KIND=DP), INTENT(IN) :: GAS_TEMPERATURE,GRAIN_TEMPERATURE
      END FUNCTION H2_FORMATION_RATE

      FUNCTION H2_PHOTODISSOCIATION_RATE(K0,G0,AV,NH2) RESULT(RATE)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         REAL(KIND=DP) :: RATE
         REAL(KIND=DP), INTENT(IN) :: K0,G0,AV,NH2
      END FUNCTION H2_PHOTODISSOCIATION_RATE

      FUNCTION CO_PHOTODISSOCIATION_RATE(K0,G0,AV,NCO,NH2) RESULT(RATE)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         REAL(KIND=DP) :: RATE
         REAL(KIND=DP), INTENT(IN) :: K0,G0,AV,NCO,NH2
      END FUNCTION CO_PHOTODISSOCIATION_RATE

      FUNCTION CI_PHOTOIONIZATION_RATE(K0,G0,AV,KAV,NCI,NH2,TGAS) RESULT(RATE)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         REAL(KIND=DP) :: RATE
         REAL(KIND=DP), INTENT(IN) :: K0,G0,AV,KAV,NCI,NH2,TGAS
      END FUNCTION CI_PHOTOIONIZATION_RATE

      FUNCTION SI_PHOTOIONIZATION_RATE(K0,G0,AV,KAV,NSI) RESULT(RATE)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         REAL(KIND=DP) :: RATE
         REAL(KIND=DP), INTENT(IN) :: K0,G0,AV,KAV,NSI
      END FUNCTION SI_PHOTOIONIZATION_RATE

      FUNCTION H2SHIELD_FGK(NH2,DOPW,RADW) RESULT(SHIELDING_FACTOR)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         REAL(KIND=DP) :: SHIELDING_FACTOR
         REAL(KIND=DP), INTENT(IN) :: NH2,DOPW,RADW
      END FUNCTION H2SHIELD_FGK

      FUNCTION H2SHIELD_LEE(NH2) RESULT(SHIELDING_FACTOR)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         REAL(KIND=DP) :: SHIELDING_FACTOR
         REAL(KIND=DP), INTENT(IN) :: NH2
      END FUNCTION H2SHIELD_LEE

      FUNCTION COSHIELD(NCO,NH2) RESULT(SHIELDING_FACTOR)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         REAL(KIND=DP) :: SHIELDING_FACTOR
         REAL(KIND=DP), INTENT(IN) :: NCO,NH2
      END FUNCTION COSHIELD

      FUNCTION SCATTER(AV,LAMBDA)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         REAL(KIND=DP) :: SCATTER
         REAL(KIND=DP), INTENT(IN) :: AV,LAMBDA
      END FUNCTION SCATTER

      FUNCTION XLAMBDA(LAMBDA)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         REAL(KIND=DP) :: XLAMBDA
         REAL(KIND=DP), INTENT(IN) :: LAMBDA
      END FUNCTION XLAMBDA

      FUNCTION LBAR(NCO,NH2)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         REAL(KIND=DP) :: LBAR
         REAL(KIND=DP), INTENT(IN) :: NCO,NH2
      END FUNCTION LBAR

      FUNCTION ESCAPE_PROBABILITY(NRAYS,TAU) RESULT(BETA)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         REAL(KIND=DP) :: BETA
         INTEGER(KIND=I4B), INTENT(IN) :: NRAYS
         REAL(KIND=DP),     INTENT(IN) :: TAU(0:NRAYS-1)
      END FUNCTION ESCAPE_PROBABILITY

      FUNCTION INTEGRATE(F,X,N) RESULT(F_DX)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         REAL(KIND=DP) :: F_DX
         INTEGER(KIND=I4B), INTENT(IN) :: N
         REAL(KIND=DP),     INTENT(IN) :: F(1:N),X(1:N)
      END FUNCTION INTEGRATE

   END INTERFACE

END MODULE FUNCTIONS_MODULE
!=======================================================================

!=======================================================================
MODULE SUBROUTINES_MODULE

   INTERFACE

      SUBROUTINE READ_AVAILABLE_CROSS_SECTIONS(FILENAME,REACTANT,PRODUCT, &
                                             & NREAC,NXSEC,CROSS_SECTION)
         USE HEALPIX_TYPES
         USE CROSS_SECTION_MODULE
         IMPLICIT NONE
         CHARACTER(LEN=*),         INTENT(IN)  :: FILENAME
         INTEGER(KIND=I4B),        INTENT(IN)  :: NREAC,NXSEC
         CHARACTER(LEN=10),        INTENT(IN)  :: REACTANT(:,:),PRODUCT(:,:)
         TYPE(CROSS_SECTION_TYPE), INTENT(OUT) :: CROSS_SECTION(:)
      END SUBROUTINE READ_AVAILABLE_CROSS_SECTIONS

      SUBROUTINE CALCULATE_REACTION_RATES(NRAYS,NSPEC,NREAC,NXSEC,CROSS_SECTION, &
                                        & GAS_TEMPERATURE,DUST_TEMPERATURE,FUV_FIELD, &
                                        & XRAY_FIELD,FUV_SURFACE,AV,COLUMN_DENSITY, &
                                        & REACTANT,PRODUCT,ALPHA,BETA,GAMMA,RTMIN,RTMAX, &
                                        & DUPLICATE,RATE)
         USE HEALPIX_TYPES
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
      END SUBROUTINE CALCULATE_REACTION_RATES

      SUBROUTINE CALCULATE_XRAY_IONIZATION_RATES(ENERGY,FLUX,RATE, &
                                               & REACTANT,NREAC)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         INTEGER(KIND=I4B), INTENT(IN)    :: NREAC
         CHARACTER(LEN=10), INTENT(IN)    :: REACTANT(:,:)
         REAL(KIND=DP),     INTENT(IN)    :: ENERGY(:),FLUX(:)
         REAL(KIND=DP),     INTENT(INOUT) :: RATE(:)
      END SUBROUTINE CALCULATE_XRAY_IONIZATION_RATES

      SUBROUTINE CALCULATE_COLLISIONAL_RATES(COOLANT,DENSITY,TEMPERATURE, &
                                           & ABUNDANCE,COLLISIONAL_RATE)
         USE HEALPIX_TYPES
         USE COOLANT_MODULE
         IMPLICIT NONE
         TYPE(COOLANT_TYPE), INTENT(IN)  :: COOLANT
         REAL(KIND=DP),      INTENT(IN)  :: DENSITY,TEMPERATURE
         REAL(KIND=DP),      INTENT(IN)  :: ABUNDANCE(:)
         REAL(KIND=DP),      INTENT(OUT) :: COLLISIONAL_RATE(:,:)
      END SUBROUTINE CALCULATE_COLLISIONAL_RATES

      SUBROUTINE UPDATE_ABUNDANCES(TEMPERATURE_ITERATION,CHEMISTRY_ITERATION,NPART,PARTICLE,INITIAL_ABUNDANCE, &
                                 & NCHEM,DUMMY_INDEX,DUMMY_ABUNDANCE,DUMMY_RATE,DUMMY_DENSITY,DUMMY_TEMPERATURE)
         USE HEALPIX_TYPES
         USE PARTICLE_MODULE
         IMPLICIT NONE
         INTEGER(KIND=I4B),   INTENT(IN)  :: TEMPERATURE_ITERATION,CHEMISTRY_ITERATION,NPART
         TYPE(PARTICLE_TYPE), INTENT(IN)  :: PARTICLE(*)
         REAL(KIND=DP),       INTENT(IN)  :: INITIAL_ABUNDANCE(:)
         INTEGER(KIND=I4B),   INTENT(OUT) :: NCHEM,DUMMY_INDEX(:)
         REAL(KIND=DP),       INTENT(OUT) :: DUMMY_ABUNDANCE(:,:),DUMMY_RATE(:,:)
         REAL(KIND=DP),       INTENT(OUT) :: DUMMY_DENSITY(:),DUMMY_TEMPERATURE(:)
      END SUBROUTINE UPDATE_ABUNDANCES

      SUBROUTINE CALCULATE_LTE_POPULATIONS(NLEVEL,ENERGY,WEIGHT,DENSITY,TEMPERATURE,POPULATION)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         INTEGER(KIND=I4B), INTENT(IN)  :: NLEVEL
         REAL(KIND=DP),     INTENT(IN)  :: ENERGY(:),WEIGHT(:)
         REAL(KIND=DP),     INTENT(IN)  :: DENSITY,TEMPERATURE
         REAL(KIND=DP),     INTENT(OUT) :: POPULATION(:)
      END SUBROUTINE CALCULATE_LTE_POPULATIONS

      SUBROUTINE GAUSS_JORDAN(N,A,B)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         INTEGER(KIND=I4B), INTENT(IN)    :: N
         REAL(KIND=DP),     INTENT(INOUT) :: A(:,:)
         REAL(KIND=DP),     INTENT(INOUT) :: B(:)
      END SUBROUTINE GAUSS_JORDAN

      SUBROUTINE CHECK_CHEMISTRY_CONVERGENCE(NPART,NSPEC,PARTICLE,PERCENTAGE_CONVERGED)
         USE HEALPIX_TYPES
         USE PARTICLE_MODULE
         IMPLICIT NONE
         INTEGER(KIND=I4B),   INTENT(IN)    :: NPART,NSPEC
         TYPE(PARTICLE_TYPE), INTENT(INOUT) :: PARTICLE(:)
         REAL(KIND=DP),       INTENT(OUT)   :: PERCENTAGE_CONVERGED
      END SUBROUTINE CHECK_CHEMISTRY_CONVERGENCE

      SUBROUTINE CHECK_POPULATION_CONVERGENCE(NPART,NCOOL,COOLANT,PARTICLE,PERCENTAGE_CONVERGED)
         USE HEALPIX_TYPES
         USE COOLANT_MODULE
         USE PARTICLE_MODULE
         IMPLICIT NONE
         INTEGER(KIND=I4B),   INTENT(IN)    :: NPART,NCOOL
         TYPE(COOLANT_TYPE),  INTENT(INOUT) :: COOLANT(:)
         TYPE(PARTICLE_TYPE), INTENT(INOUT) :: PARTICLE(:)
         REAL(KIND=DP),       INTENT(OUT)   :: PERCENTAGE_CONVERGED(:)
      END SUBROUTINE CHECK_POPULATION_CONVERGENCE

      SUBROUTINE ANALYSE_CHEMISTRY(PREFIX,SUFFIX,NPART,NSPEC,NREAC,TIME, &
                                 & SPECIES,REACTANT,PRODUCT,PARTICLE, &
                                 & AV_LIST,SPECIES_LIST)
         USE HEALPIX_TYPES
         USE PARTICLE_MODULE
         IMPLICIT NONE
         CHARACTER(LEN=*),    INTENT(IN) :: PREFIX,SUFFIX
         INTEGER(KIND=I4B),   INTENT(IN) :: NPART,NSPEC,NREAC
         REAL(KIND=DP),       INTENT(IN) :: TIME,AV_LIST(:)
         CHARACTER(LEN=*),    INTENT(IN) :: SPECIES(:),SPECIES_LIST(:)
         CHARACTER(LEN=*),    INTENT(IN) :: REACTANT(:,:),PRODUCT(:,:)
         TYPE(PARTICLE_TYPE), INTENT(IN) :: PARTICLE(:)
      END SUBROUTINE ANALYSE_CHEMISTRY

   END INTERFACE

END MODULE SUBROUTINES_MODULE
!=======================================================================

!=======================================================================
MODULE OUTPUT_MODULE

   INTERFACE

      SUBROUTINE WRITE_PROPERTIES(PREFIX,SUFFIX,NPART,PARTICLE)
         USE HEALPIX_TYPES
         USE PARTICLE_MODULE
         IMPLICIT NONE
         CHARACTER(LEN=*),    INTENT(IN) :: PREFIX,SUFFIX
         INTEGER(KIND=I4B),   INTENT(IN) :: NPART
         TYPE(PARTICLE_TYPE), INTENT(IN) :: PARTICLE(:)
      END SUBROUTINE WRITE_PROPERTIES

      SUBROUTINE WRITE_EXTINCTION(PREFIX,SUFFIX,NPART,NRAYS,PARTICLE)
         USE HEALPIX_TYPES
         USE PARTICLE_MODULE
         IMPLICIT NONE
         CHARACTER(LEN=*),    INTENT(IN) :: PREFIX,SUFFIX
         INTEGER(KIND=I4B),   INTENT(IN) :: NPART,NRAYS
         TYPE(PARTICLE_TYPE), INTENT(IN) :: PARTICLE(:)
      END SUBROUTINE WRITE_EXTINCTION

      SUBROUTINE WRITE_REACTION_RATES(PREFIX,SUFFIX,NPART,NREAC,REACTANT,PARTICLE)
         USE HEALPIX_TYPES
         USE PARTICLE_MODULE
         IMPLICIT NONE
         CHARACTER(LEN=*),    INTENT(IN) :: PREFIX,SUFFIX
         INTEGER(KIND=I4B),   INTENT(IN) :: NPART,NREAC
         CHARACTER(LEN=*),    INTENT(IN) :: REACTANT(:,:)
         TYPE(PARTICLE_TYPE), INTENT(IN) :: PARTICLE(:)
      END SUBROUTINE WRITE_REACTION_RATES

      SUBROUTINE WRITE_ABUNDANCES(PREFIX,SUFFIX,NPART,NSPEC,SPECIES,PARTICLE)
         USE HEALPIX_TYPES
         USE PARTICLE_MODULE
         IMPLICIT NONE
         CHARACTER(LEN=*),    INTENT(IN) :: PREFIX,SUFFIX
         INTEGER(KIND=I4B),   INTENT(IN) :: NPART,NSPEC
         CHARACTER(LEN=*),    INTENT(IN) :: SPECIES(:)
         TYPE(PARTICLE_TYPE), INTENT(IN) :: PARTICLE(:)
      END SUBROUTINE WRITE_ABUNDANCES

      SUBROUTINE WRITE_POPULATIONS(PREFIX,SUFFIX,NPART,NCOOL,COOLANT,PARTICLE)
         USE HEALPIX_TYPES
         USE COOLANT_MODULE
         USE PARTICLE_MODULE
         IMPLICIT NONE
         CHARACTER(LEN=*),    INTENT(IN) :: PREFIX,SUFFIX
         INTEGER(KIND=I4B),   INTENT(IN) :: NPART,NCOOL
         TYPE(COOLANT_TYPE),  INTENT(IN) :: COOLANT(:)
         TYPE(PARTICLE_TYPE), INTENT(IN) :: PARTICLE(:)
      END SUBROUTINE WRITE_POPULATIONS

      SUBROUTINE WRITE_OPACITIES(PREFIX,SUFFIX,NPART,NCOOL,NRAYS,RAY,COOLANT,PARTICLE)
         USE HEALPIX_TYPES
         USE COOLANT_MODULE
         USE PARTICLE_MODULE
         IMPLICIT NONE
         CHARACTER(LEN=*),    INTENT(IN) :: PREFIX,SUFFIX
         INTEGER(KIND=I4B),   INTENT(IN) :: NPART,NCOOL
         INTEGER(KIND=I4B),   INTENT(IN) :: NRAYS,RAY
         TYPE(COOLANT_TYPE),  INTENT(IN) :: COOLANT(:)
         TYPE(PARTICLE_TYPE), INTENT(IN) :: PARTICLE(:)
      END SUBROUTINE WRITE_OPACITIES

      SUBROUTINE WRITE_EMISSIVITIES(PREFIX,SUFFIX,NPART,NCOOL,COOLANT,PARTICLE)
         USE HEALPIX_TYPES
         USE COOLANT_MODULE
         USE PARTICLE_MODULE
         IMPLICIT NONE
         CHARACTER(LEN=*),    INTENT(IN) :: PREFIX,SUFFIX
         INTEGER(KIND=I4B),   INTENT(IN) :: NPART,NCOOL
         TYPE(COOLANT_TYPE),  INTENT(IN) :: COOLANT(:)
         TYPE(PARTICLE_TYPE), INTENT(IN) :: PARTICLE(:)
      END SUBROUTINE WRITE_EMISSIVITIES

      SUBROUTINE WRITE_COOLING_RATES(PREFIX,SUFFIX,NPART,NCOOL,COOLANT,PARTICLE)
         USE HEALPIX_TYPES
         USE COOLANT_MODULE
         USE PARTICLE_MODULE
         IMPLICIT NONE
         CHARACTER(LEN=*),    INTENT(IN) :: PREFIX,SUFFIX
         INTEGER(KIND=I4B),   INTENT(IN) :: NPART,NCOOL
         TYPE(COOLANT_TYPE),  INTENT(IN) :: COOLANT(:)
         TYPE(PARTICLE_TYPE), INTENT(IN) :: PARTICLE(:)
      END SUBROUTINE WRITE_COOLING_RATES

      SUBROUTINE WRITE_HEATING_RATES(PREFIX,SUFFIX,NPART,NHEAT,PARTICLE)
         USE HEALPIX_TYPES
         USE PARTICLE_MODULE
         IMPLICIT NONE
         CHARACTER(LEN=*),    INTENT(IN) :: PREFIX,SUFFIX
         INTEGER(KIND=I4B),   INTENT(IN) :: NPART,NHEAT
         TYPE(PARTICLE_TYPE), INTENT(IN) :: PARTICLE(:)
      END SUBROUTINE WRITE_HEATING_RATES

      SUBROUTINE WRITE_THERMAL_BALANCE(PREFIX,SUFFIX,NPART,PARTICLE)
         USE HEALPIX_TYPES
         USE PARTICLE_MODULE
         IMPLICIT NONE
         CHARACTER(LEN=*),    INTENT(IN) :: PREFIX,SUFFIX
         INTEGER(KIND=I4B),   INTENT(IN) :: NPART
         TYPE(PARTICLE_TYPE), INTENT(IN) :: PARTICLE(:)
      END SUBROUTINE WRITE_THERMAL_BALANCE

      SUBROUTINE WRITE_CONVERGENCE_STATUS(PREFIX,SUFFIX,NPART,NCOOL,COOLANT,PARTICLE)
         USE HEALPIX_TYPES
         USE COOLANT_MODULE
         USE PARTICLE_MODULE
         IMPLICIT NONE
         CHARACTER(LEN=*),    INTENT(IN) :: PREFIX,SUFFIX
         INTEGER(KIND=I4B),   INTENT(IN) :: NPART,NCOOL
         TYPE(COOLANT_TYPE),  INTENT(IN) :: COOLANT(:)
         TYPE(PARTICLE_TYPE), INTENT(IN) :: PARTICLE(:)
      END SUBROUTINE WRITE_CONVERGENCE_STATUS

   END INTERFACE

END MODULE OUTPUT_MODULE

MODULE TRANSITION_MATRIX_MODULE

   INTERFACE

      SUBROUTINE CONSTRUCT_TRANSITION_MATRIX(NRAYS,N,COOLANT,PARTICLE,TRANSITION_MATRIX)
         USE HEALPIX_TYPES
         USE COOLANT_MODULE
         USE PARTICLE_MODULE
         IMPLICIT NONE
         INTEGER(KIND=I4B),   INTENT(IN)    :: NRAYS
         INTEGER(KIND=I4B),   INTENT(IN)    :: N
         TYPE(COOLANT_TYPE),  INTENT(IN)    :: COOLANT
         TYPE(PARTICLE_TYPE), INTENT(INOUT) :: PARTICLE
         REAL(KIND=DP),       INTENT(OUT)   :: TRANSITION_MATRIX(:,:)
      END SUBROUTINE CONSTRUCT_TRANSITION_MATRIX

   END INTERFACE

END MODULE TRANSITION_MATRIX_MODULE
!=======================================================================
