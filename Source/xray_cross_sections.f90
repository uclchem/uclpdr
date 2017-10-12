!=======================================================================
MODULE XRAY_FUNCTIONS

   INTERFACE

      FUNCTION XRAY_FLUX(ENERGY,SIGMA,NPOINT,T_X,F_X,N_H) RESULT(FLUX)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
         REAL(KIND=DP),     INTENT(IN)  :: ENERGY(1:NPOINT)
         REAL(KIND=DP),     INTENT(IN)  :: SIGMA(1:NPOINT)
         REAL(KIND=DP),     INTENT(IN)  :: T_X,F_X,N_H
         REAL(KIND=DP) :: FLUX(1:NPOINT)
      END FUNCTION XRAY_FLUX

      FUNCTION ENERGY_DEPOSITION_RATE(ENERGY,FLUX,SIGMA,NPOINT) RESULT(H_X)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         INTEGER(KIND=I4B), INTENT(IN) :: NPOINT
         REAL(KIND=DP),     INTENT(IN) :: ENERGY(1:NPOINT)
         REAL(KIND=DP),     INTENT(IN) :: FLUX(1:NPOINT)
         REAL(KIND=DP),     INTENT(IN) :: SIGMA(1:NPOINT)
         REAL(KIND=DP) :: H_X
      END FUNCTION ENERGY_DEPOSITION_RATE

      FUNCTION XRAY_IONIZATION_RATE(ENERGY,FLUX,SPECIES,NPOINT) RESULT(RATE)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         INTEGER(KIND=I4B), INTENT(IN) :: NPOINT
         CHARACTER(LEN=10), INTENT(IN) :: SPECIES
         REAL(KIND=DP),     INTENT(IN) :: ENERGY(1:NPOINT)
         REAL(KIND=DP),     INTENT(IN) :: FLUX(1:NPOINT)
         REAL(KIND=DP) :: RATE
      END FUNCTION XRAY_IONIZATION_RATE

      FUNCTION TOTAL_XRAY_CROSS_SECTION(ENERGY,NPOINT) RESULT(SIGMA)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
         REAL(KIND=DP),     INTENT(IN)  :: ENERGY(1:NPOINT)
         REAL(KIND=DP) :: SIGMA(1:NPOINT)
      END FUNCTION TOTAL_XRAY_CROSS_SECTION

      FUNCTION HYDROGEN_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT) RESULT(SIGMA)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
         REAL(KIND=DP),     INTENT(IN)  :: ENERGY(1:NPOINT)
         REAL(KIND=DP) :: SIGMA(1:NPOINT)
      END FUNCTION HYDROGEN_PHOTOIONIZATION_CROSS_SECTION

      FUNCTION HELIUM_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT) RESULT(SIGMA)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
         REAL(KIND=DP),     INTENT(IN)  :: ENERGY(1:NPOINT)
         REAL(KIND=DP) :: SIGMA(1:NPOINT)
      END FUNCTION HELIUM_PHOTOIONIZATION_CROSS_SECTION

      FUNCTION CARBON_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT) RESULT(SIGMA)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
         REAL(KIND=DP),     INTENT(IN)  :: ENERGY(1:NPOINT)
         REAL(KIND=DP) :: SIGMA(1:NPOINT)
      END FUNCTION CARBON_PHOTOIONIZATION_CROSS_SECTION

      FUNCTION NITROGEN_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT) RESULT(SIGMA)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
         REAL(KIND=DP),     INTENT(IN)  :: ENERGY(1:NPOINT)
         REAL(KIND=DP) :: SIGMA(1:NPOINT)
      END FUNCTION NITROGEN_PHOTOIONIZATION_CROSS_SECTION

      FUNCTION OXYGEN_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT) RESULT(SIGMA)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
         REAL(KIND=DP),     INTENT(IN)  :: ENERGY(1:NPOINT)
         REAL(KIND=DP) :: SIGMA(1:NPOINT)
      END FUNCTION OXYGEN_PHOTOIONIZATION_CROSS_SECTION

      FUNCTION NEON_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT) RESULT(SIGMA)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
         REAL(KIND=DP),     INTENT(IN)  :: ENERGY(1:NPOINT)
         REAL(KIND=DP) :: SIGMA(1:NPOINT)
      END FUNCTION NEON_PHOTOIONIZATION_CROSS_SECTION

      FUNCTION SODIUM_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT) RESULT(SIGMA)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
         REAL(KIND=DP),     INTENT(IN)  :: ENERGY(1:NPOINT)
         REAL(KIND=DP) :: SIGMA(1:NPOINT)
      END FUNCTION SODIUM_PHOTOIONIZATION_CROSS_SECTION

      FUNCTION MAGNESIUM_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT) RESULT(SIGMA)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
         REAL(KIND=DP),     INTENT(IN)  :: ENERGY(1:NPOINT)
         REAL(KIND=DP) :: SIGMA(1:NPOINT)
      END FUNCTION MAGNESIUM_PHOTOIONIZATION_CROSS_SECTION

      FUNCTION ALUMINIUM_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT) RESULT(SIGMA)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
         REAL(KIND=DP),     INTENT(IN)  :: ENERGY(1:NPOINT)
         REAL(KIND=DP) :: SIGMA(1:NPOINT)
      END FUNCTION ALUMINIUM_PHOTOIONIZATION_CROSS_SECTION

      FUNCTION SILICON_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT) RESULT(SIGMA)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
         REAL(KIND=DP),     INTENT(IN)  :: ENERGY(1:NPOINT)
         REAL(KIND=DP) :: SIGMA(1:NPOINT)
      END FUNCTION SILICON_PHOTOIONIZATION_CROSS_SECTION

      FUNCTION PHOSPHORUS_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT) RESULT(SIGMA)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
         REAL(KIND=DP),     INTENT(IN)  :: ENERGY(1:NPOINT)
         REAL(KIND=DP) :: SIGMA(1:NPOINT)
      END FUNCTION PHOSPHORUS_PHOTOIONIZATION_CROSS_SECTION

      FUNCTION SULPHUR_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT) RESULT(SIGMA)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
         REAL(KIND=DP),     INTENT(IN)  :: ENERGY(1:NPOINT)
         REAL(KIND=DP) :: SIGMA(1:NPOINT)
      END FUNCTION SULPHUR_PHOTOIONIZATION_CROSS_SECTION

      FUNCTION CHLORINE_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT) RESULT(SIGMA)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
         REAL(KIND=DP),     INTENT(IN)  :: ENERGY(1:NPOINT)
         REAL(KIND=DP) :: SIGMA(1:NPOINT)
      END FUNCTION CHLORINE_PHOTOIONIZATION_CROSS_SECTION

      FUNCTION ARGON_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT) RESULT(SIGMA)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
         REAL(KIND=DP),     INTENT(IN)  :: ENERGY(1:NPOINT)
         REAL(KIND=DP) :: SIGMA(1:NPOINT)
      END FUNCTION ARGON_PHOTOIONIZATION_CROSS_SECTION

      FUNCTION CALCIUM_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT) RESULT(SIGMA)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
         REAL(KIND=DP),     INTENT(IN)  :: ENERGY(1:NPOINT)
         REAL(KIND=DP) :: SIGMA(1:NPOINT)
      END FUNCTION CALCIUM_PHOTOIONIZATION_CROSS_SECTION

      FUNCTION CHROMIUM_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT) RESULT(SIGMA)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
         REAL(KIND=DP),     INTENT(IN)  :: ENERGY(1:NPOINT)
         REAL(KIND=DP) :: SIGMA(1:NPOINT)
      END FUNCTION CHROMIUM_PHOTOIONIZATION_CROSS_SECTION

      FUNCTION IRON_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT) RESULT(SIGMA)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
         REAL(KIND=DP),     INTENT(IN)  :: ENERGY(1:NPOINT)
         REAL(KIND=DP) :: SIGMA(1:NPOINT)
      END FUNCTION IRON_PHOTOIONIZATION_CROSS_SECTION

      FUNCTION NICKEL_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT) RESULT(SIGMA)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
         REAL(KIND=DP),     INTENT(IN)  :: ENERGY(1:NPOINT)
         REAL(KIND=DP) :: SIGMA(1:NPOINT)
      END FUNCTION NICKEL_PHOTOIONIZATION_CROSS_SECTION

      FUNCTION INTEGRATE(F,X,N) RESULT(F_DX)
         USE HEALPIX_TYPES
         IMPLICIT NONE
         REAL(KIND=DP) :: F_DX
         INTEGER(KIND=I4B), INTENT(IN) :: N
         REAL(KIND=DP),     INTENT(IN) :: F(1:N),X(1:N)
      END FUNCTION INTEGRATE

   END INTERFACE

END MODULE XRAY_FUNCTIONS
!=======================================================================

!=======================================================================
!
!  Define the range of X-ray energies and determine the corresponding
!  photoionization cross sections, sigma (cm^2), based on the assumed
!  elemental abundances, then calculate the local X-ray flux (in ergs
!  cm^-2 s^-1) and integrated energy deposition rate, H_X (erg s^-1),
!  for each particle based on the incident fluxes along each ray.
!
!-----------------------------------------------------------------------
SUBROUTINE CALCULATE_XRAY_PROPERTIES(NPART,NRAYS,PARTICLE)

   USE HEALPIX_TYPES
   USE PARTICLE_MODULE
   USE XRAY_FUNCTIONS

   IMPLICIT NONE

   INTEGER(KIND=I4B),   INTENT(IN)    :: NPART,NRAYS
   TYPE(PARTICLE_TYPE), INTENT(INOUT) :: PARTICLE(*)

   INTEGER(KIND=I4B) :: I,J,P,NPOINT
   REAL(KIND=DP) :: E_STEP,E_MIN,E_MAX,XRAY_TEMPERATURE,H_X
   REAL(KIND=DP), ALLOCATABLE :: ENERGY(:),SIGMA(:),FLUX(:)

!  Define the effective temperature (K) of the X-ray spectrum
   XRAY_TEMPERATURE = 1.16D7

!  Define the array of X-ray energies (eV) at which cross sections and fluxes will be calculated
   E_STEP = 1.0D0 ! 1 eV
   E_MIN  = 1.0D3 ! 1 keV
   E_MAX  = 1.0D4 ! 10 keV
   NPOINT = INT((E_MAX - E_MIN)/E_STEP) + 1

   ALLOCATE(ENERGY(1:NPOINT),SIGMA(1:NPOINT),FLUX(1:NPOINT))

   DO I=1,NPOINT
      ENERGY(I) = E_MIN + E_STEP*(I-1)
   END DO

!  Calculate the photoionization absorption cross section at each energy
   SIGMA = TOTAL_XRAY_CROSS_SECTION(ENERGY,NPOINT)

   DO P=1,NPART ! Loop over all particles

!     Allocate the X-ray photon energy and flux arrays
      ALLOCATE(PARTICLE(P)%XRAY_ENERGIES(1:NPOINT))
      ALLOCATE(PARTICLE(P)%XRAY_FLUXES(1:NPOINT))

!     Store the X-ray photon energies
      PARTICLE(P)%XRAY_ENERGIES = ENERGY

!     Calculate the contribution to the local X-ray flux and the energy
!     deposition rate, H_X, from the incident radiation along each ray
      PARTICLE(P)%XRAY_FLUX = 0.0D0
      PARTICLE(P)%XRAY_FLUXES = 0.0D0
      PARTICLE(P)%XRAY_ENERGY_DEPOSITION_RATE = 0.0D0
      DO J=0,NRAYS-1
         FLUX = XRAY_FLUX(ENERGY,SIGMA,NPOINT,XRAY_TEMPERATURE,PARTICLE(P)%XRAY_SURFACE(J),PARTICLE(P)%TOTAL_COLUMN(J))
         H_X  = ENERGY_DEPOSITION_RATE(ENERGY,FLUX,SIGMA,NPOINT)
         PARTICLE(P)%XRAY_FLUX = PARTICLE(P)%XRAY_FLUX + INTEGRATE(FLUX,ENERGY,NPOINT)
         PARTICLE(P)%XRAY_FLUXES = PARTICLE(P)%XRAY_FLUXES + FLUX
         PARTICLE(P)%XRAY_ENERGY_DEPOSITION_RATE = PARTICLE(P)%XRAY_ENERGY_DEPOSITION_RATE + H_X
      END DO

!     Impose a lower cut-off of 1E-99 for the energy deposition rate H_X
      IF(PARTICLE(P)%XRAY_ENERGY_DEPOSITION_RATE.LT.1.0D-99) PARTICLE(P)%XRAY_ENERGY_DEPOSITION_RATE = 0.0D0

      IF(PARTICLE(P)%XRAY_FLUX.GE.1.0D-10) PARTICLE(P)%TYPE = 2

   END DO ! End of loop over particles

   RETURN
END SUBROUTINE CALCULATE_XRAY_PROPERTIES
!=======================================================================

!=======================================================================
!
!  Calculate the local X-ray flux (erg cm^-2 s^-1) at specified energies
!  E (eV), due to emission from a thermal plasma with unattenuated total
!  X-ray flux F_X (erg cm^-2 s^-1), temperature T_X (K), photoionization
!  absorption cross sections at the specified energies, σ_pa(E) (cm^2),
!  and total hydrogen column density N_H (cm^-2).
!
!-----------------------------------------------------------------------
FUNCTION XRAY_FLUX(ENERGY,SIGMA,NPOINT,T_X,F_X,N_H) RESULT(FLUX)

   USE HEALPIX_TYPES

   IMPLICIT NONE

   INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
   REAL(KIND=DP),     INTENT(IN)  :: ENERGY(1:NPOINT)
   REAL(KIND=DP),     INTENT(IN)  :: SIGMA(1:NPOINT)
   REAL(KIND=DP),     INTENT(IN)  :: T_X,F_X,N_H

   REAL(KIND=DP) :: FLUX(1:NPOINT)
   REAL(KIND=DP) :: NORMALIZATION

!  Initialize the X-ray fluxes
   FLUX = 0.0D0

!  Calculate the normalization factor (convert the temperature from K to eV by the factor kB/eV)
   NORMALIZATION = (KB/EV*T_X)*(EXP(-MINVAL(ENERGY)/(KB/EV*T_X)) - EXP(-MAXVAL(ENERGY)/(KB/EV*T_X)))

!  J(E) = F_X / N . exp(-E/kT) . exp(-σ(E).N_H) [erg cm^-2 s^-1]
   FLUX = F_X/NORMALIZATION*EXP(-ENERGY/(KB/EV*T_X))*EXP(-SIGMA*N_H)

   RETURN
END FUNCTION XRAY_FLUX
!=======================================================================

!=======================================================================
!
!  Calculate the integrated energy deposition rate per hydrogen nucleus
!  H_X (erg s^-1) due to X-ray photon absorption. The energy-dependent
!  X-ray flux (erg cm^-2 s^-1), photoelectric absorption cross section
!  σ_pa(E) (cm^2), and energy values at which they are both specified
!  must be supplied.
!
!-----------------------------------------------------------------------
FUNCTION ENERGY_DEPOSITION_RATE(ENERGY,FLUX,SIGMA,NPOINT) RESULT(H_X)

   USE HEALPIX_TYPES
   USE FUNCTIONS_MODULE

   IMPLICIT NONE

   INTEGER(KIND=I4B), INTENT(IN) :: NPOINT
   REAL(KIND=DP),     INTENT(IN) :: FLUX(1:NPOINT)
   REAL(KIND=DP),     INTENT(IN) :: SIGMA(1:NPOINT)
   REAL(KIND=DP),     INTENT(IN) :: ENERGY(1:NPOINT)

   REAL(KIND=DP) :: H_X

   H_X = INTEGRATE(SIGMA*FLUX,ENERGY,NPOINT)

   RETURN
END FUNCTION ENERGY_DEPOSITION_RATE
!=======================================================================

!=======================================================================
!
!  Integrate the supplied function f(x) using the trapezium rule for the
!  specified array of x values, [x_1, ..., x_N].
!
!-----------------------------------------------------------------------
FUNCTION INTEGRATE(F,X,N) RESULT(F_DX)

   USE HEALPIX_TYPES

   IMPLICIT NONE

   INTEGER(KIND=I4B), INTENT(IN) :: N
   REAL(KIND=DP),     INTENT(IN) :: F(1:N),X(1:N)

   INTEGER(KIND=I4B) :: I
   REAL(KIND=DP) :: F_DX

!  Initialize the output value
   F_DX = 0.0D0

   DO I=1,N-1 ! Loop over x-values
      F_DX = F_DX + 0.5D0*(F(I+1)+F(I))*(X(I+1)-X(I))
   END DO ! End of loop over x-values

   RETURN
END FUNCTION INTEGRATE
!=======================================================================

!=======================================================================
!
!  Calculate photoelectric absorption cross sections per hydrogen
!  nucleus, σ_pa (cm^2), using the treatment of Meijerink & Spaans
!  (2005, A&A, 436, 397, Appendix E), based on the analytical fits
!  to photoionization cross sections σ_i(E) from Verner & Yakovlev
!  (1995, A&AS, 109, 125) and assumed total (gas + dust) elemental
!  abundances.
!
!-----------------------------------------------------------------------
FUNCTION TOTAL_XRAY_CROSS_SECTION(ENERGY,NPOINT) RESULT(SIGMA)

   USE HEALPIX_TYPES
   USE XRAY_FUNCTIONS

   IMPLICIT NONE

   INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
   REAL(KIND=DP),     INTENT(IN)  :: ENERGY(1:NPOINT)

   REAL(KIND=DP) :: SIGMA(1:NPOINT)
   REAL(KIND=DP) :: TOTAL_ABUNDANCE

!  Initialize the cross sections
   SIGMA = 0.0D0

!  Calculate the photoelectric absorption cross section, σ_pa (cm^2), at each energy
!  by adding the contribution from the atomic ground electronic state of each element

!  Hydrogen
   TOTAL_ABUNDANCE = 1.0D+00
   SIGMA = SIGMA + TOTAL_ABUNDANCE*HYDROGEN_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)

!  Helium
   TOTAL_ABUNDANCE = 8.5D-02
   SIGMA = SIGMA + TOTAL_ABUNDANCE*HELIUM_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)

!  Carbon
   TOTAL_ABUNDANCE = 2.5D-04
   SIGMA = SIGMA + TOTAL_ABUNDANCE*CARBON_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)

!  Nitrogen
   TOTAL_ABUNDANCE = 7.2D-05
   SIGMA = SIGMA + TOTAL_ABUNDANCE*NITROGEN_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)

!  Oxygen
   TOTAL_ABUNDANCE = 4.7D-04
   SIGMA = SIGMA + TOTAL_ABUNDANCE*OXYGEN_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)

!  Neon
   TOTAL_ABUNDANCE = 6.9D-05
   SIGMA = SIGMA + TOTAL_ABUNDANCE*NEON_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)

!  Sodium
   TOTAL_ABUNDANCE = 1.5D-06
   SIGMA = SIGMA + TOTAL_ABUNDANCE*SODIUM_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)

!  Magnesium
   TOTAL_ABUNDANCE = 3.4D-06
   SIGMA = SIGMA + TOTAL_ABUNDANCE*MAGNESIUM_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)

!  Aluminium
   TOTAL_ABUNDANCE = 2.3D-06
   SIGMA = SIGMA + TOTAL_ABUNDANCE*ALUMINIUM_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)

!  Silicon
   TOTAL_ABUNDANCE = 3.4D-05
   SIGMA = SIGMA + TOTAL_ABUNDANCE*SILICON_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)

!  Phosphorus
   TOTAL_ABUNDANCE = 2.9D-07
   SIGMA = SIGMA + TOTAL_ABUNDANCE*PHOSPHORUS_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)

!  Sulphur
   TOTAL_ABUNDANCE = 1.4D-05
   SIGMA = SIGMA + TOTAL_ABUNDANCE*SULPHUR_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)

!  Chlorine
   TOTAL_ABUNDANCE = 2.4D-07
   SIGMA = SIGMA + TOTAL_ABUNDANCE*CHLORINE_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)

!  Argon
   TOTAL_ABUNDANCE = 1.5D-06
   SIGMA = SIGMA + TOTAL_ABUNDANCE*ARGON_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)

!  Calcium
   TOTAL_ABUNDANCE = 2.0D-06
   SIGMA = SIGMA + TOTAL_ABUNDANCE*CALCIUM_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)

!  Chromium
   TOTAL_ABUNDANCE = 4.4D-07
   SIGMA = SIGMA + TOTAL_ABUNDANCE*CHROMIUM_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)

!  Iron
   TOTAL_ABUNDANCE = 2.8D-05
   SIGMA = SIGMA + TOTAL_ABUNDANCE*IRON_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)

!  Nickel
   TOTAL_ABUNDANCE = 1.7D-06
   SIGMA = SIGMA + TOTAL_ABUNDANCE*NICKEL_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT)

   RETURN
END FUNCTION TOTAL_XRAY_CROSS_SECTION
!=======================================================================

!=======================================================================
!
!  A set of functions to calculate the photoionization cross sections,
!  σ_i (cm^2), at the specified photon energies, E(eV), for the atomic
!  ground electronic state of each element from H to Zn (Z ≤ 30) using
!  the analytical fits of Verner & Yakovlev (1995, A&AS, 109, 125).
!
!  Only the ground electronic state of each atom is considered.
!
!-----------------------------------------------------------------------
FUNCTION HYDROGEN_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT) RESULT(SIGMA)

   USE HEALPIX_TYPES

   IMPLICIT NONE

   INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
   REAL(KIND=DP),     INTENT(IN)  :: ENERGY(1:NPOINT)

   REAL(KIND=DP) :: SIGMA(1:NPOINT)

   INTEGER(KIND=I4B) :: I
   REAL(KIND=DP) :: E_th,E_0,SIGMA_0
   REAL(KIND=DP) :: L,P,Q,Y,Y_A,Y_W

!  Initialize the cross sections
   SIGMA = 0.0D0

!  Hydrogen (n=1,l=0)
   L       = 0.000D+00
   E_th    = 1.360D+01
   E_0     = 4.298D-01
   SIGMA_0 = 5.475D+04
   Y_A     = 3.288D+01
   Y_W     = 0.000D+00
   P       = 2.963D+00
   Q = 5.5D0+L-0.5D0*P

   DO I=1,NPOINT
      IF(ENERGY(I).GE.E_th) THEN
         Y = ENERGY(I)/E_0
         SIGMA(I) = (SIGMA_0*1.0D-18)*((Y-1.0D0)**2+Y_W**2)*Y**(-Q)*(1.0D0+SQRT(Y/Y_A))**(-P)
      END IF
   END DO

   RETURN
END FUNCTION HYDROGEN_PHOTOIONIZATION_CROSS_SECTION
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
FUNCTION HELIUM_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT) RESULT(SIGMA)

   USE HEALPIX_TYPES

   IMPLICIT NONE

   INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
   REAL(KIND=DP),     INTENT(IN)  :: ENERGY(1:NPOINT)

   REAL(KIND=DP) :: SIGMA(1:NPOINT)

   INTEGER(KIND=I4B) :: I
   REAL(KIND=DP) :: E_th,E_0,SIGMA_0
   REAL(KIND=DP) :: L,P,Q,Y,Y_A,Y_W

!  Initialize the cross sections
   SIGMA = 0.0D0

!  Helium (n=1,l=0)
   L       = 0.000D+00
   E_th    = 2.459D+01
   E_0     = 5.996D+00
   SIGMA_0 = 4.470D+03
   Y_A     = 2.199D+00
   Y_W     = 0.000D+00
   P       = 6.098D+00
   Q = 5.5D0+L-0.5D0*P

   DO I=1,NPOINT
      IF(ENERGY(I).GE.E_th) THEN
         Y = ENERGY(I)/E_0
         SIGMA(I) = (SIGMA_0*1.0D-18)*((Y-1.0D0)**2+Y_W**2)*Y**(-Q)*(1.0D0+SQRT(Y/Y_A))**(-P)
      END IF
   END DO

   RETURN
END FUNCTION HELIUM_PHOTOIONIZATION_CROSS_SECTION
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
FUNCTION CARBON_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT) RESULT(SIGMA)

   USE HEALPIX_TYPES

   IMPLICIT NONE

   INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
   REAL(KIND=DP),     INTENT(IN)  :: ENERGY(1:NPOINT)

   REAL(KIND=DP) :: SIGMA(1:NPOINT)

   INTEGER(KIND=I4B) :: I
   REAL(KIND=DP) :: E_th,E_0,SIGMA_0
   REAL(KIND=DP) :: L,P,Q,Y,Y_A,Y_W

!  Initialize the cross sections
   SIGMA = 0.0D0

!  Carbon (n=2,l=1)
   L       = 1.000D+00
   E_th    = 1.126D+01
   E_0     = 9.435D+00
   SIGMA_0 = 1.152D+03
   Y_A     = 5.687D+00
   Y_W     = 4.474D-01
   P       = 6.336D+00
   Q = 5.5D0+L-0.5D0*P

!!$!  Carbon (n=2,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 1.939D+01
!!$   E_0     = 1.026D+01
!!$   SIGMA_0 = 4.564D+03
!!$   Y_A     = 1.568D+00
!!$   Y_W     = 0.000D+00
!!$   P       = 1.085D+01
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Carbon (n=1,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 2.910D+02
!!$   E_0     = 8.655D+01
!!$   SIGMA_0 = 7.421D+01
!!$   Y_A     = 5.498D+01
!!$   Y_W     = 0.000D+00
!!$   P       = 1.503D+00
!!$   Q = 5.5D0+L-0.5D0*P

   DO I=1,NPOINT
      IF(ENERGY(I).GE.E_th) THEN
         Y = ENERGY(I)/E_0
         SIGMA(I) = (SIGMA_0*1.0D-18)*((Y-1.0D0)**2+Y_W**2)*Y**(-Q)*(1.0D0+SQRT(Y/Y_A))**(-P)
      END IF
   END DO

   RETURN
END FUNCTION CARBON_PHOTOIONIZATION_CROSS_SECTION
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
FUNCTION NITROGEN_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT) RESULT(SIGMA)

   USE HEALPIX_TYPES

   IMPLICIT NONE

   INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
   REAL(KIND=DP),     INTENT(IN)  :: ENERGY(1:NPOINT)

   REAL(KIND=DP) :: SIGMA(1:NPOINT)

   INTEGER(KIND=I4B) :: I
   REAL(KIND=DP) :: E_th,E_0,SIGMA_0
   REAL(KIND=DP) :: L,P,Q,Y,Y_A,Y_W

!  Initialize the cross sections
   SIGMA = 0.0D0

!  Nitrogen (n=2,l=1)
   L       = 1.000D+00
   E_th    = 1.453D+01
   E_0     = 1.164D+01
   SIGMA_0 = 1.029D+04
   Y_A     = 2.361D+00
   Y_W     = 4.239D-01
   P       = 8.821D+00
   Q = 5.5D0+L-0.5D0*P

!!$!  Nitrogen (n=2,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 2.541D+01
!!$   E_0     = 1.482D+01
!!$   SIGMA_0 = 7.722D+02
!!$   Y_A     = 2.306D+00
!!$   Y_W     = 0.000D+00
!!$   P       = 9.139D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Nitrogen (n=1,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 4.048D+02
!!$   E_0     = 1.270D+02
!!$   SIGMA_0 = 4.748D+01
!!$   Y_A     = 1.380D+02
!!$   Y_W     = 0.000D+00
!!$   P       = 1.252D+00
!!$   Q = 5.5D0+L-0.5D0*P

   DO I=1,NPOINT
      IF(ENERGY(I).GE.E_th) THEN
         Y = ENERGY(I)/E_0
         SIGMA(I) = (SIGMA_0*1.0D-18)*((Y-1.0D0)**2+Y_W**2)*Y**(-Q)*(1.0D0+SQRT(Y/Y_A))**(-P)
      END IF
   END DO

   RETURN
END FUNCTION NITROGEN_PHOTOIONIZATION_CROSS_SECTION
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
FUNCTION OXYGEN_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT) RESULT(SIGMA)

   USE HEALPIX_TYPES

   IMPLICIT NONE

   INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
   REAL(KIND=DP),     INTENT(IN)  :: ENERGY(1:NPOINT)

   REAL(KIND=DP) :: SIGMA(1:NPOINT)

   INTEGER(KIND=I4B) :: I
   REAL(KIND=DP) :: E_th,E_0,SIGMA_0
   REAL(KIND=DP) :: L,P,Q,Y,Y_A,Y_W

!  Initialize the cross sections
   SIGMA = 0.0D0

!  Oxygen (n=2,l=1)
   L       = 1.000D+00
   E_th    = 1.362D+01
   E_0     = 1.391D+01
   SIGMA_0 = 1.220D+05
   Y_A     = 1.364D+00
   Y_W     = 4.103D-01
   P       = 1.140D+01
   Q = 5.5D0+L-0.5D0*P

!!$!  Oxygen (n=2,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 2.848D+01
!!$   E_0     = 1.994D+01
!!$   SIGMA_0 = 2.415D+02
!!$   Y_A     = 3.241D+00
!!$   Y_W     = 0.000D+00
!!$   P       = 8.037D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Oxygen (n=1,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 5.380D+02
!!$   E_0     = 1.774D+02
!!$   SIGMA_0 = 3.237D+01
!!$   Y_A     = 3.812D+02
!!$   Y_W     = 0.000D+00
!!$   P       = 1.083D+00
!!$   Q = 5.5D0+L-0.5D0*P

   DO I=1,NPOINT
      IF(ENERGY(I).GE.E_th) THEN
         Y = ENERGY(I)/E_0
         SIGMA(I) = (SIGMA_0*1.0D-18)*((Y-1.0D0)**2+Y_W**2)*Y**(-Q)*(1.0D0+SQRT(Y/Y_A))**(-P)
      END IF
   END DO

   RETURN
END FUNCTION OXYGEN_PHOTOIONIZATION_CROSS_SECTION
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
FUNCTION NEON_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT) RESULT(SIGMA)

   USE HEALPIX_TYPES

   IMPLICIT NONE

   INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
   REAL(KIND=DP),     INTENT(IN)  :: ENERGY(1:NPOINT)

   REAL(KIND=DP) :: SIGMA(1:NPOINT)

   INTEGER(KIND=I4B) :: I
   REAL(KIND=DP) :: E_th,E_0,SIGMA_0
   REAL(KIND=DP) :: L,P,Q,Y,Y_A,Y_W

!  Initialize the cross sections
   SIGMA = 0.0D0

!  Neon (n=2,l=1)
   L       = 1.000D+00
   E_th    = 2.156D+01
   E_0     = 2.000D+01
   SIGMA_0 = 1.691D+04
   Y_A     = 2.442D+00
   Y_W     = 3.345D-01
   P       = 1.043D+01
   Q = 5.5D0+L-0.5D0*P

!!$!  Neon (n=2,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 4.847D+01
!!$   E_0     = 3.204D+01
!!$   SIGMA_0 = 5.615D+01
!!$   Y_A     = 5.808D+00
!!$   Y_W     = 0.000D+00
!!$   P       = 6.678D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Neon (n=1,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 8.701D+02
!!$   E_0     = 3.144D+02
!!$   SIGMA_0 = 1.664D+01
!!$   Y_A     = 2.042D+05
!!$   Y_W     = 0.000D+00
!!$   P       = 8.450D-01
!!$   Q = 5.5D0+L-0.5D0*P

   DO I=1,NPOINT
      IF(ENERGY(I).GE.E_th) THEN
         Y = ENERGY(I)/E_0
         SIGMA(I) = (SIGMA_0*1.0D-18)*((Y-1.0D0)**2+Y_W**2)*Y**(-Q)*(1.0D0+SQRT(Y/Y_A))**(-P)
      END IF
   END DO

   RETURN
END FUNCTION NEON_PHOTOIONIZATION_CROSS_SECTION
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
FUNCTION SODIUM_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT) RESULT(SIGMA)

   USE HEALPIX_TYPES

   IMPLICIT NONE

   INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
   REAL(KIND=DP),     INTENT(IN)  :: ENERGY(1:NPOINT)

   REAL(KIND=DP) :: SIGMA(1:NPOINT)

   INTEGER(KIND=I4B) :: I
   REAL(KIND=DP) :: E_th,E_0,SIGMA_0
   REAL(KIND=DP) :: L,P,Q,Y,Y_A,Y_W

!  Initialize the cross sections
   SIGMA = 0.0D0

!  Sodium (n=3,l=0)
   L       = 0.000D+00
   E_th    = 5.139D+00
   E_0     = 5.968D+00
   SIGMA_0 = 1.460D+00
   Y_A     = 2.557D+07
   Y_W     = 0.000D+00
   P       = 3.789D+00
   Q = 5.5D0+L-0.5D0*P

!!$!  Sodium (n=2,l=1)
!!$   L       = 1.000D+00
!!$   E_th    = 3.814D+01
!!$   E_0     = 3.655D+01
!!$   SIGMA_0 = 2.486D+02
!!$   Y_A     = 3.222D+02
!!$   Y_W     = 1.465D-01
!!$   P       = 3.570D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Sodium (n=2,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 7.084D+01
!!$   E_0     = 4.537D+01
!!$   SIGMA_0 = 1.142D+01
!!$   Y_A     = 2.395D+02
!!$   Y_W     = 0.000D+00
!!$   P       = 3.380D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Sodium (n=1,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 1.079D+03
!!$   E_0     = 4.216D+02
!!$   SIGMA_0 = 1.119D+01
!!$   Y_A     = 5.642D+07
!!$   Y_W     = 0.000D+00
!!$   P       = 7.736D-01
!!$   Q = 5.5D0+L-0.5D0*P

   DO I=1,NPOINT
      IF(ENERGY(I).GE.E_th) THEN
         Y = ENERGY(I)/E_0
         SIGMA(I) = (SIGMA_0*1.0D-18)*((Y-1.0D0)**2+Y_W**2)*Y**(-Q)*(1.0D0+SQRT(Y/Y_A))**(-P)
      END IF
   END DO

   RETURN
END FUNCTION SODIUM_PHOTOIONIZATION_CROSS_SECTION
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
FUNCTION MAGNESIUM_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT) RESULT(SIGMA)

   USE HEALPIX_TYPES

   IMPLICIT NONE

   INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
   REAL(KIND=DP),     INTENT(IN)  :: ENERGY(1:NPOINT)

   REAL(KIND=DP) :: SIGMA(1:NPOINT)

   INTEGER(KIND=I4B) :: I
   REAL(KIND=DP) :: E_th,E_0,SIGMA_0
   REAL(KIND=DP) :: L,P,Q,Y,Y_A,Y_W

!  Initialize the cross sections
   SIGMA = 0.0D0

!  Magnesium (n=3,l=0)
   L       = 0.000D+00
   E_th    = 7.646D+00
   E_0     = 9.393D+00
   SIGMA_0 = 3.034D+00
   Y_A     = 2.625D+07
   Y_W     = 0.000D+00
   P       = 3.923D+00
   Q = 5.5D0+L-0.5D0*P

!!$!  Magnesium (n=2,l=1)
!!$   L       = 1.000D+00
!!$   E_th    = 5.490D+01
!!$   E_0     = 4.937D+01
!!$   SIGMA_0 = 2.023D+02
!!$   Y_A     = 1.079D+04
!!$   Y_W     = 1.463D-02
!!$   P       = 2.960D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Magnesium (n=2,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 9.400D+01
!!$   E_0     = 4.587D+01
!!$   SIGMA_0 = 1.671D+01
!!$   Y_A     = 2.389D+01
!!$   Y_W     = 0.000D+00
!!$   P       = 4.742D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Magnesium (n=1,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 1.311D+03
!!$   E_0     = 2.711D+02
!!$   SIGMA_0 = 3.561D+01
!!$   Y_A     = 2.374D+01
!!$   Y_W     = 0.000D+00
!!$   P       = 1.952D+00
!!$   Q = 5.5D0+L-0.5D0*P

   DO I=1,NPOINT
      IF(ENERGY(I).GE.E_th) THEN
         Y = ENERGY(I)/E_0
         SIGMA(I) = (SIGMA_0*1.0D-18)*((Y-1.0D0)**2+Y_W**2)*Y**(-Q)*(1.0D0+SQRT(Y/Y_A))**(-P)
      END IF
   END DO

   RETURN
END FUNCTION MAGNESIUM_PHOTOIONIZATION_CROSS_SECTION
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
FUNCTION ALUMINIUM_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT) RESULT(SIGMA)

   USE HEALPIX_TYPES

   IMPLICIT NONE

   INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
   REAL(KIND=DP),     INTENT(IN)  :: ENERGY(1:NPOINT)

   REAL(KIND=DP) :: SIGMA(1:NPOINT)

   INTEGER(KIND=I4B) :: I
   REAL(KIND=DP) :: E_th,E_0,SIGMA_0
   REAL(KIND=DP) :: L,P,Q,Y,Y_A,Y_W

!  Initialize the cross sections
   SIGMA = 0.0D0

!  Aluminium (n=3,l=1)
   L       = 1.000D+00
   E_th    = 5.986D+00
   E_0     = 1.860D+01
   SIGMA_0 = 1.828D+02
   Y_A     = 2.797D+00
   Y_W     = 3.076D-01
   P       = 1.084D+01
   Q = 5.5D0+L-0.5D0*P

!!$!  Aluminium (n=3,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 1.133D+01
!!$   E_0     = 1.204D+01
!!$   SIGMA_0 = 5.384D+00
!!$   Y_A     = 4.341D+02
!!$   Y_W     = 0.000D+00
!!$   P       = 4.088D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Aluminium (n=2,l=1)
!!$   L       = 1.000D+00
!!$   E_th    = 8.040D+01
!!$   E_0     = 6.445D+01
!!$   SIGMA_0 = 1.735D+02
!!$   Y_A     = 1.131D+04
!!$   Y_W     = 2.337D-02
!!$   P       = 2.762D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Aluminium (n=2,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 1.256D+02
!!$   E_0     = 5.594D+01
!!$   SIGMA_0 = 1.425D+01
!!$   Y_A     = 3.094D+01
!!$   Y_W     = 0.000D+00
!!$   P       = 4.399D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Aluminium (n=1,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 1.567D+03
!!$   E_0     = 3.670D+02
!!$   SIGMA_0 = 2.206D+01
!!$   Y_A     = 4.405D+01
!!$   Y_W     = 0.000D+00
!!$   P       = 1.588D+00
!!$   Q = 5.5D0+L-0.5D0*P

   DO I=1,NPOINT
      IF(ENERGY(I).GE.E_th) THEN
         Y = ENERGY(I)/E_0
         SIGMA(I) = (SIGMA_0*1.0D-18)*((Y-1.0D0)**2+Y_W**2)*Y**(-Q)*(1.0D0+SQRT(Y/Y_A))**(-P)
      END IF
   END DO

   RETURN
END FUNCTION ALUMINIUM_PHOTOIONIZATION_CROSS_SECTION
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
FUNCTION SILICON_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT) RESULT(SIGMA)

   USE HEALPIX_TYPES

   IMPLICIT NONE

   INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
   REAL(KIND=DP),     INTENT(IN)  :: ENERGY(1:NPOINT)

   REAL(KIND=DP) :: SIGMA(1:NPOINT)

   INTEGER(KIND=I4B) :: I
   REAL(KIND=DP) :: E_th,E_0,SIGMA_0
   REAL(KIND=DP) :: L,P,Q,Y,Y_A,Y_W

!  Initialize the cross sections
   SIGMA = 0.0D0

!  Silicon (n=3,l=1)
   L       = 1.000D+00
   E_th    = 8.152D+00
   E_0     = 2.212D+01
   SIGMA_0 = 1.845D+02
   Y_A     = 3.849D+00
   Y_W     = 2.921D-01
   P       = 9.721D+00
   Q = 5.5D0+L-0.5D0*P

!!$!  Silicon (n=3,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 1.517D+01
!!$   E_0     = 1.413D+01
!!$   SIGMA_0 = 1.166D+01
!!$   Y_A     = 2.288D+01
!!$   Y_W     = 0.000D+00
!!$   P       = 5.334D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Silicon (n=2,l=1)
!!$   L       = 1.000D+00
!!$   E_th    = 1.060D+02
!!$   E_0     = 7.808D+01
!!$   SIGMA_0 = 1.532D+02
!!$   Y_A     = 5.765D+06
!!$   Y_W     = 2.774D-04
!!$   P       = 2.639D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Silicon (n=2,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 1.560D+02
!!$   E_0     = 7.017D+01
!!$   SIGMA_0 = 1.166D+01
!!$   Y_A     = 4.742D+01
!!$   Y_W     = 0.000D+00
!!$   P       = 3.933D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Silicon (n=1,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 1.846D+03
!!$   E_0     = 5.322D+02
!!$   SIGMA_0 = 1.184D+01
!!$   Y_A     = 2.580D+02
!!$   Y_W     = 0.000D+00
!!$   P       = 1.102D+00
!!$   Q = 5.5D0+L-0.5D0*P

   DO I=1,NPOINT
      IF(ENERGY(I).GE.E_th) THEN
         Y = ENERGY(I)/E_0
         SIGMA(I) = (SIGMA_0*1.0D-18)*((Y-1.0D0)**2+Y_W**2)*Y**(-Q)*(1.0D0+SQRT(Y/Y_A))**(-P)
      END IF
   END DO

   RETURN
END FUNCTION SILICON_PHOTOIONIZATION_CROSS_SECTION
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
FUNCTION PHOSPHORUS_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT) RESULT(SIGMA)

   USE HEALPIX_TYPES

   IMPLICIT NONE

   INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
   REAL(KIND=DP),     INTENT(IN)  :: ENERGY(1:NPOINT)

   REAL(KIND=DP) :: SIGMA(1:NPOINT)

   INTEGER(KIND=I4B) :: I
   REAL(KIND=DP) :: E_th,E_0,SIGMA_0
   REAL(KIND=DP) :: L,P,Q,Y,Y_A,Y_W

!  Initialize the cross sections
   SIGMA = 0.0D0

!  Phosphorus (n=3,l=1)
   L       = 1.000D+00
   E_th    = 1.049D+01
   E_0     = 2.580D+01
   SIGMA_0 = 9.925D+01
   Y_A     = 6.712D+00
   Y_W     = 2.765D-01
   P       = 8.516D+00
   Q = 5.5D0+L-0.5D0*P

!!$!  Phosphorus (n=3,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 2.017D+01
!!$   E_0     = 1.658D+01
!!$   SIGMA_0 = 1.125D+01
!!$   Y_A     = 2.613D+01
!!$   Y_W     = 0.000D+00
!!$   P       = 5.205D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Phosphorus (n=2,l=1)
!!$   L       = 1.000D+00
!!$   E_th    = 1.400D+02
!!$   E_0     = 8.812D+01
!!$   SIGMA_0 = 1.512D+02
!!$   Y_A     = 2.230D+03
!!$   Y_W     = 3.422D-04
!!$   P       = 2.795D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Phosphorus (n=2,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 1.940D+02
!!$   E_0     = 8.632D+01
!!$   SIGMA_0 = 9.931D+00
!!$   Y_A     = 6.594D+01
!!$   Y_W     = 0.000D+00
!!$   P       = 3.617D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Phosphorus (n=1,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 2.154D+03
!!$   E_0     = 6.472D+02
!!$   SIGMA_0 = 9.167D+00
!!$   Y_A     = 2.562D+02
!!$   Y_W     = 0.000D+00
!!$   P       = 1.063D+00
!!$   Q = 5.5D0+L-0.5D0*P

   DO I=1,NPOINT
      IF(ENERGY(I).GE.E_th) THEN
         Y = ENERGY(I)/E_0
         SIGMA(I) = (SIGMA_0*1.0D-18)*((Y-1.0D0)**2+Y_W**2)*Y**(-Q)*(1.0D0+SQRT(Y/Y_A))**(-P)
      END IF
   END DO

   RETURN
END FUNCTION PHOSPHORUS_PHOTOIONIZATION_CROSS_SECTION
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
FUNCTION SULPHUR_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT) RESULT(SIGMA)

   USE HEALPIX_TYPES

   IMPLICIT NONE

   INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
   REAL(KIND=DP),     INTENT(IN)  :: ENERGY(1:NPOINT)

   REAL(KIND=DP) :: SIGMA(1:NPOINT)

   INTEGER(KIND=I4B) :: I
   REAL(KIND=DP) :: E_th,E_0,SIGMA_0
   REAL(KIND=DP) :: L,P,Q,Y,Y_A,Y_W

!  Initialize the cross sections
   SIGMA = 0.0D0

!  Sulphur (n=3,l=1)
   L       = 1.000D+00
   E_th    = 1.036D+01
   E_0     = 2.975D+01
   SIGMA_0 = 5.644D+01
   Y_A     = 1.321D+01
   Y_W     = 2.621D-01
   P       = 7.513D+00
   Q = 5.5D0+L-0.5D0*P

!!$!  Sulphur (n=3,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 2.130D+01
!!$   E_0     = 1.916D+01
!!$   SIGMA_0 = 1.003D+01
!!$   Y_A     = 3.296D+01
!!$   Y_W     = 0.000D+00
!!$   P       = 5.038D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Sulphur (n=2,l=1)
!!$   L       = 1.000D+00
!!$   E_th    = 1.700D+02
!!$   E_0     = 9.152D+01
!!$   SIGMA_0 = 1.883D+02
!!$   Y_A     = 7.193D+01
!!$   Y_W     = 2.485D-01
!!$   P       = 3.633D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Sulphur (n=2,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 2.350D+02
!!$   E_0     = 1.047D+02
!!$   SIGMA_0 = 8.520D+00
!!$   Y_A     = 9.469D+01
!!$   Y_W     = 0.000D+00
!!$   P       = 3.346D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Sulphur (n=1,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 2.477D+03
!!$   E_0     = 8.114D+02
!!$   SIGMA_0 = 6.649D+00
!!$   Y_A     = 3.734D+03
!!$   Y_W     = 0.000D+00
!!$   P       = 8.646D-01
!!$   Q = 5.5D0+L-0.5D0*P

   DO I=1,NPOINT
      IF(ENERGY(I).GE.E_th) THEN
         Y = ENERGY(I)/E_0
         SIGMA(I) = (SIGMA_0*1.0D-18)*((Y-1.0D0)**2+Y_W**2)*Y**(-Q)*(1.0D0+SQRT(Y/Y_A))**(-P)
      END IF
   END DO

   RETURN
END FUNCTION SULPHUR_PHOTOIONIZATION_CROSS_SECTION
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
FUNCTION CHLORINE_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT) RESULT(SIGMA)

   USE HEALPIX_TYPES

   IMPLICIT NONE

   INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
   REAL(KIND=DP),     INTENT(IN)  :: ENERGY(1:NPOINT)

   REAL(KIND=DP) :: SIGMA(1:NPOINT)

   INTEGER(KIND=I4B) :: I
   REAL(KIND=DP) :: E_th,E_0,SIGMA_0
   REAL(KIND=DP) :: L,P,Q,Y,Y_A,Y_W

!  Initialize the cross sections
   SIGMA = 0.0D0

!  Chlorine (n=3,l=1)
   L       = 1.000D+00
   E_th    = 1.297D+01
   E_0     = 3.398D+01
   SIGMA_0 = 4.539D+01
   Y_A     = 2.232D+01
   Y_W     = 2.479D-01
   P       = 6.896D+00
   Q = 5.5D0+L-0.5D0*P

!!$!  Chlorine (n=3,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 2.531D+01
!!$   E_0     = 2.231D+01
!!$   SIGMA_0 = 6.628D+00
!!$   Y_A     = 1.843D+02
!!$   Y_W     = 0.000D+00
!!$   P       = 4.196D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Chlorine (n=2,l=1)
!!$   L       = 1.000D+00
!!$   E_th    = 2.090D+02
!!$   E_0     = 8.004D+01
!!$   SIGMA_0 = 3.053D+02
!!$   Y_A     = 3.498D+01
!!$   Y_W     = 2.017D-01
!!$   P       = 4.457D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Chlorine (n=2,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 2.780D+02
!!$   E_0     = 1.092D+02
!!$   SIGMA_0 = 1.059D+01
!!$   Y_A     = 2.491D+01
!!$   Y_W     = 0.000D+00
!!$   P       = 4.205D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Chlorine (n=1,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 2.830D+03
!!$   E_0     = 9.700D+02
!!$   SIGMA_0 = 5.255D+00
!!$   Y_A     = 1.856D+06
!!$   Y_W     = 0.000D+00
!!$   P       = 7.888D-01
!!$   Q = 5.5D0+L-0.5D0*P

   DO I=1,NPOINT
      IF(ENERGY(I).GE.E_th) THEN
         Y = ENERGY(I)/E_0
         SIGMA(I) = (SIGMA_0*1.0D-18)*((Y-1.0D0)**2+Y_W**2)*Y**(-Q)*(1.0D0+SQRT(Y/Y_A))**(-P)
      END IF
   END DO

   RETURN
END FUNCTION CHLORINE_PHOTOIONIZATION_CROSS_SECTION
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
FUNCTION ARGON_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT) RESULT(SIGMA)

   USE HEALPIX_TYPES

   IMPLICIT NONE

   INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
   REAL(KIND=DP),     INTENT(IN)  :: ENERGY(1:NPOINT)

   REAL(KIND=DP) :: SIGMA(1:NPOINT)

   INTEGER(KIND=I4B) :: I
   REAL(KIND=DP) :: E_th,E_0,SIGMA_0
   REAL(KIND=DP) :: L,P,Q,Y,Y_A,Y_W

!  Initialize the cross sections
   SIGMA = 0.0D0

!  Argon (n=3,l=1)
   L       = 1.000D+00
   E_th    = 1.576D+01
   E_0     = 3.854D+01
   SIGMA_0 = 4.872D+01
   Y_A     = 2.640D+01
   Y_W     = 2.355D-01
   P       = 6.662D+00
   Q = 5.5D0+L-0.5D0*P

!!$!  Argon (n=3,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 2.892D+01
!!$   E_0     = 2.525D+01
!!$   SIGMA_0 = 6.394D+00
!!$   Y_A     = 1.700D+02
!!$   Y_W     = 0.000D+00
!!$   P       = 4.223D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Argon (n=2,l=1)
!!$   L       = 1.000D+00
!!$   E_th    = 2.492D+02
!!$   E_0     = 1.647D+02
!!$   SIGMA_0 = 8.372D+01
!!$   Y_A     = 5.452D+01
!!$   Y_W     = 6.270D-01
!!$   P       = 3.328D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Argon (n=2,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 3.260D+02
!!$   E_0     = 1.302D+02
!!$   SIGMA_0 = 9.185D+00
!!$   Y_A     = 2.693D+01
!!$   Y_W     = 0.000D+00
!!$   P       = 4.021D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Argon (n=1,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 3.203D+03
!!$   E_0     = 1.135D+03
!!$   SIGMA_0 = 4.280D+00
!!$   Y_A     = 3.285D+07
!!$   Y_W     = 0.000D+00
!!$   P       = 7.631D-01
!!$   Q = 5.5D0+L-0.5D0*P

   DO I=1,NPOINT
      IF(ENERGY(I).GE.E_th) THEN
         Y = ENERGY(I)/E_0
         SIGMA(I) = (SIGMA_0*1.0D-18)*((Y-1.0D0)**2+Y_W**2)*Y**(-Q)*(1.0D0+SQRT(Y/Y_A))**(-P)
      END IF
   END DO

   RETURN
END FUNCTION ARGON_PHOTOIONIZATION_CROSS_SECTION
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
FUNCTION CALCIUM_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT) RESULT(SIGMA)

   USE HEALPIX_TYPES

   IMPLICIT NONE

   INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
   REAL(KIND=DP),     INTENT(IN)  :: ENERGY(1:NPOINT)

   REAL(KIND=DP) :: SIGMA(1:NPOINT)

   INTEGER(KIND=I4B) :: I
   REAL(KIND=DP) :: E_th,E_0,SIGMA_0
   REAL(KIND=DP) :: L,P,Q,Y,Y_A,Y_W

!  Initialize the cross sections
   SIGMA = 0.0D0

!  Calcium (n=4,l=0)
   L       = 0.000D+00
   E_th    = 6.113D+00
   E_0     = 7.366D+00
   SIGMA_0 = 2.373D+00
   Y_A     = 2.082D+02
   Y_W     = 5.841D-04
   P       = 4.841D+00
   Q = 5.5D0+L-0.5D0*P

!!$!  Calcium (n=3,l=1)
!!$   L       = 1.000D+00
!!$   E_th    = 3.443D+01
!!$   E_0     = 4.487D+01
!!$   SIGMA_0 = 9.017D+01
!!$   Y_A     = 1.465D+01
!!$   Y_W     = 2.754D-01
!!$   P       = 7.498D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Calcium (n=3,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 4.830D+01
!!$   E_0     = 3.012D+01
!!$   SIGMA_0 = 7.227D+00
!!$   Y_A     = 1.736D+02
!!$   Y_W     = 0.000D+00
!!$   P       = 4.165D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Calcium (n=2,l=1)
!!$   L       = 1.000D+00
!!$   E_th    = 3.523D+02
!!$   E_0     = 1.529D+02
!!$   SIGMA_0 = 1.282D+02
!!$   Y_A     = 2.217D+02
!!$   Y_W     = 3.343D-03
!!$   P       = 3.087D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Calcium (n=2,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 4.425D+02
!!$   E_0     = 1.201D+02
!!$   SIGMA_0 = 1.010D+01
!!$   Y_A     = 2.468D+01
!!$   Y_W     = 0.000D+00
!!$   P       = 4.592D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Calcium (n=1,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 4.043D+03
!!$   E_0     = 6.947D+02
!!$   SIGMA_0 = 1.586D+01
!!$   Y_A     = 2.563D+01
!!$   Y_W     = 0.000D+00
!!$   P       = 1.966D+00
!!$   Q = 5.5D0+L-0.5D0*P

   DO I=1,NPOINT
      IF(ENERGY(I).GE.E_th) THEN
         Y = ENERGY(I)/E_0
         SIGMA(I) = (SIGMA_0*1.0D-18)*((Y-1.0D0)**2+Y_W**2)*Y**(-Q)*(1.0D0+SQRT(Y/Y_A))**(-P)
      END IF
   END DO

   RETURN
END FUNCTION CALCIUM_PHOTOIONIZATION_CROSS_SECTION
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
FUNCTION CHROMIUM_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT) RESULT(SIGMA)

   USE HEALPIX_TYPES

   IMPLICIT NONE

   INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
   REAL(KIND=DP),     INTENT(IN)  :: ENERGY(1:NPOINT)

   REAL(KIND=DP) :: SIGMA(1:NPOINT)

   INTEGER(KIND=I4B) :: I
   REAL(KIND=DP) :: E_th,E_0,SIGMA_0
   REAL(KIND=DP) :: L,P,Q,Y,Y_A,Y_W

!  Initialize the cross sections
   SIGMA = 0.0D0

!  Chromium (n=4,l=0)
   L       = 0.000D+00
   E_th    = 6.767D+00
   E_0     = 9.636D+00
   SIGMA_0 = 6.532D-01
   Y_A     = 5.232D+02
   Y_W     = 9.332D-05
   P       = 4.641D+00
   Q = 5.5D0+L-0.5D0*P

!!$!  Chromium (n=3,l=2)
!!$   L       = 2.000D+00
!!$   E_th    = 8.660D+00
!!$   E_0     = 7.244D+00
!!$   SIGMA_0 = 1.485D+03
!!$   Y_A     = 9.671D+00
!!$   Y_W     = 7.760D-01
!!$   P       = 1.575D+01
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Chromium (n=3,l=1)
!!$   L       = 1.000D+00
!!$   E_th    = 4.900D+01
!!$   E_0     = 6.567D+01
!!$   SIGMA_0 = 5.313D+01
!!$   Y_A     = 1.981D+01
!!$   Y_W     = 2.601D-01
!!$   P       = 7.258D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Chromium (n=3,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 7.900D+01
!!$   E_0     = 4.234D+01
!!$   SIGMA_0 = 5.602D+00
!!$   Y_A     = 1.356D+02
!!$   Y_W     = 0.000D+00
!!$   P       = 4.374D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Chromium (n=2,l=1)
!!$   L       = 1.000D+00
!!$   E_th    = 5.850D+02
!!$   E_0     = 1.865D+02
!!$   SIGMA_0 = 1.540D+02
!!$   Y_A     = 5.560D+01
!!$   Y_W     = 5.785D-03
!!$   P       = 3.823D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Chromium (n=2,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 7.030D+02
!!$   E_0     = 1.588D+02
!!$   SIGMA_0 = 8.555D+00
!!$   Y_A     = 2.258D+01
!!$   Y_W     = 0.000D+00
!!$   P       = 4.789D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Chromium (n=1,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 5.996D+03
!!$   E_0     = 1.035D+03
!!$   SIGMA_0 = 1.025D+01
!!$   Y_A     = 3.343D+01
!!$   Y_W     = 0.000D+00
!!$   P       = 1.822D+00
!!$   Q = 5.5D0+L-0.5D0*P

   DO I=1,NPOINT
      IF(ENERGY(I).GE.E_th) THEN
         Y = ENERGY(I)/E_0
         SIGMA(I) = (SIGMA_0*1.0D-18)*((Y-1.0D0)**2+Y_W**2)*Y**(-Q)*(1.0D0+SQRT(Y/Y_A))**(-P)
      END IF
   END DO

   RETURN
END FUNCTION CHROMIUM_PHOTOIONIZATION_CROSS_SECTION
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
FUNCTION IRON_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT) RESULT(SIGMA)

   USE HEALPIX_TYPES

   IMPLICIT NONE

   INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
   REAL(KIND=DP),     INTENT(IN)  :: ENERGY(1:NPOINT)

   REAL(KIND=DP) :: SIGMA(1:NPOINT)

   INTEGER(KIND=I4B) :: I
   REAL(KIND=DP) :: E_th,E_0,SIGMA_0
   REAL(KIND=DP) :: L,P,Q,Y,Y_A,Y_W

!  Initialize the cross sections
   SIGMA = 0.0D0

!  Iron (n=4,l=0)
   L       = 0.000D+00
   E_th    = 7.902D+00
   E_0     = 1.277D+01
   SIGMA_0 = 1.468D+00
   Y_A     = 1.116D+05
   Y_W     = 3.238D-02
   P       = 4.112D+00
   Q = 5.5D0+L-0.5D0*P

!!$!  Iron (n=3,l=2)
!!$   L       = 2.000D+00
!!$   E_th    = 1.470D+01
!!$   E_0     = 1.407D+01
!!$   SIGMA_0 = 1.850D+04
!!$   Y_A     = 4.458D+00
!!$   Y_W     = 4.039D-01
!!$   P       = 1.691D+01
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Iron (n=3,l=1)
!!$   L       = 1.000D+00
!!$   E_th    = 6.600D+01
!!$   E_0     = 7.630D+01
!!$   SIGMA_0 = 6.298D+01
!!$   Y_A     = 1.479D+01
!!$   Y_W     = 2.646D-01
!!$   P       = 7.672D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Iron (n=3,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 1.040D+02
!!$   E_0     = 4.334D+01
!!$   SIGMA_0 = 5.921D+00
!!$   Y_A     = 5.293D+01
!!$   Y_W     = 0.000D+00
!!$   P       = 5.129D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Iron (n=2,l=1)
!!$   L       = 1.000D+00
!!$   E_th    = 7.240D+02
!!$   E_0     = 2.948D+02
!!$   SIGMA_0 = 7.191D+01
!!$   Y_A     = 3.219D+02
!!$   Y_W     = 6.314D-02
!!$   P       = 2.837D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Iron (n=2,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 8.570D+02
!!$   E_0     = 5.727D+01
!!$   SIGMA_0 = 1.076D+01
!!$   Y_A     = 2.785D+01
!!$   Y_W     = 0.000D+00
!!$   P       = 6.635D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Iron (n=1,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 7.124D+03
!!$   E_0     = 8.044D+02
!!$   SIGMA_0 = 2.055D+01
!!$   Y_A     = 3.633D+01
!!$   Y_W     = 0.000D+00
!!$   P       = 2.118D+00
!!$   Q = 5.5D0+L-0.5D0*P

   DO I=1,NPOINT
      IF(ENERGY(I).GE.E_th) THEN
         Y = ENERGY(I)/E_0
         SIGMA(I) = (SIGMA_0*1.0D-18)*((Y-1.0D0)**2+Y_W**2)*Y**(-Q)*(1.0D0+SQRT(Y/Y_A))**(-P)
      END IF
   END DO

   RETURN
END FUNCTION IRON_PHOTOIONIZATION_CROSS_SECTION
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
FUNCTION NICKEL_PHOTOIONIZATION_CROSS_SECTION(ENERGY,NPOINT) RESULT(SIGMA)

   USE HEALPIX_TYPES

   IMPLICIT NONE

   INTEGER(KIND=I4B), INTENT(IN)  :: NPOINT
   REAL(KIND=DP),     INTENT(IN)  :: ENERGY(1:NPOINT)

   REAL(KIND=DP) :: SIGMA(1:NPOINT)

   INTEGER(KIND=I4B) :: I
   REAL(KIND=DP) :: E_th,E_0,SIGMA_0
   REAL(KIND=DP) :: L,P,Q,Y,Y_A,Y_W

!  Initialize the cross sections
   SIGMA = 0.0D0

!  Nickel (n=4,l=0)
   L       = 0.000D+00
   E_th    = 7.637D+00
   E_0     = 1.468D+01
   SIGMA_0 = 1.437D+00
   Y_A     = 7.411D+02
   Y_W     = 3.908D-02
   P       = 4.342D+00
   Q = 5.5D0+L-0.5D0*P

!!$!  Nickel (n=3,l=2)
!!$   L       = 2.000D+00
!!$   E_th    = 1.700D+01
!!$   E_0     = 6.063D+00
!!$   SIGMA_0 = 1.186D+03
!!$   Y_A     = 6.823D+00
!!$   Y_W     = 6.227D-03
!!$   P       = 2.223D+01
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Nickel (n=3,l=1)
!!$   L       = 1.000D+00
!!$   E_th    = 8.200D+01
!!$   E_0     = 8.896D+01
!!$   SIGMA_0 = 4.351D+01
!!$   Y_A     = 1.942D+01
!!$   Y_W     = 2.566D-01
!!$   P       = 7.372D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Nickel (n=3,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 1.250D+02
!!$   E_0     = 5.448D+01
!!$   SIGMA_0 = 4.611D+00
!!$   Y_A     = 1.157D+02
!!$   Y_W     = 0.000D+00
!!$   P       = 4.548D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Nickel (n=2,l=1)
!!$   L       = 1.000D+00
!!$   E_th    = 8.760D+02
!!$   E_0     = 3.043D+02
!!$   SIGMA_0 = 8.611D+01
!!$   Y_A     = 7.868D+01
!!$   Y_W     = 1.680D-05
!!$   P       = 3.408D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Nickel (n=2,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 1.024D+03
!!$   E_0     = 1.132D+02
!!$   SIGMA_0 = 9.424D+00
!!$   Y_A     = 2.712D+01
!!$   Y_W     = 0.000D+00
!!$   P       = 5.643D+00
!!$   Q = 5.5D0+L-0.5D0*P

!!$!  Nickel (n=1,l=0)
!!$   L       = 0.000D+00
!!$   E_th    = 8.348D+03
!!$   E_0     = 7.366D+02
!!$   SIGMA_0 = 2.836D+01
!!$   Y_A     = 3.622D+01
!!$   Y_W     = 0.000D+00
!!$   P       = 2.316D+00
!!$   Q = 5.5D0+L-0.5D0*P

   DO I=1,NPOINT
      IF(ENERGY(I).GE.E_th) THEN
         Y = ENERGY(I)/E_0
         SIGMA(I) = (SIGMA_0*1.0D-18)*((Y-1.0D0)**2+Y_W**2)*Y**(-Q)*(1.0D0+SQRT(Y/Y_A))**(-P)
      END IF
   END DO

   RETURN
END FUNCTION NICKEL_PHOTOIONIZATION_CROSS_SECTION
!=======================================================================
