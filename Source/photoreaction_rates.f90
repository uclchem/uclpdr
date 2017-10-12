!=======================================================================
!
!  H2 photodissociation rate taking into
!  account shielding and grain extinction
!
!-----------------------------------------------------------------------
!
!  Input parameters:
!  K0  = Unattenuated photodissociation rate (in cm^3/s)
!  G0  = Incident FUV field (in Draine units)
!  AV  = visual extinction (in magnitudes)
!  NH2 = H2 column density (in cm^-2)
!
!  Program variables:
!  RATE   = H2 photodissociation rate taking into
!           account self-shielding and grain extinction
!  LAMBDA = wavelength (in Å) of a typical transition
!  DOPW   = Doppler linewidth (in Hz) of a typical transition
!           (assuming turbulent broadening with b=3 km/s)
!  RADW   = radiative linewidth (in Hz) of a typical transition
!
!  Functions called:
!  SCATTER  = attenuation due to scattering by dust
!  H2SHIELD = H2 shielding function (FGK or LEE)
!
!-----------------------------------------------------------------------
FUNCTION H2_PHOTODISSOCIATION_RATE(K0,G0,AV,NH2) RESULT(RATE)

   USE HEALPIX_TYPES
   USE GLOBAL_MODULE, ONLY : V_TURB
   IMPLICIT NONE

   REAL(KIND=DP) :: RATE
   REAL(KIND=DP), INTENT(IN) :: K0,G0,AV,NH2

   REAL(KIND=DP) :: LAMBDA,DOPW,RADW,SCATTER
   REAL(KIND=DP) :: H2SHIELD_FGK,H2SHIELD_LEE

   LAMBDA=1000.0D0
   DOPW=V_TURB/(LAMBDA*1.0D-8)
   RADW=8.0D7

!  Calculate the H2 photodissociation rate
   RATE=K0*G0*SCATTER(AV,LAMBDA)*H2SHIELD_FGK(NH2,DOPW,RADW)
!   RATE=K0*G0*SCATTER(AV,LAMBDA)*H2SHIELD_LEE(NH2)
!   RATE=K0*G0*H2SHIELD_LEE(NH2)

   RETURN
END FUNCTION H2_PHOTODISSOCIATION_RATE
!=======================================================================

!=======================================================================
!
!  CO photodissociation rate taking into
!  account shielding and grain extinction
!
!-----------------------------------------------------------------------
!
!  Input parameters:
!  K0  = Unattenuated photodissociation rate (in cm^3/s)
!  G0  = Incident FUV field (in Draine units)
!  AV  = visual extinction (in magnitudes)
!  NCO = CO column density (in cm^-2)
!  NH2 = H2 column density (in cm^-2)
!
!  Program variables:
!  RATE   = CO photodissociation rate taking into
!           account self-shielding and grain extinction
!  LAMBDA = wavelength (in Å) of a typical transition
!
!  Functions called:
!  LBAR     = function to determine the wavelength
!  SCATTER  = attenuation due to scattering by dust
!  COSHIELD = CO shielding function
!
!-----------------------------------------------------------------------
FUNCTION CO_PHOTODISSOCIATION_RATE(K0,G0,AV,NCO,NH2) RESULT(RATE)

   USE HEALPIX_TYPES
   IMPLICIT NONE

   REAL(KIND=DP) :: RATE
   REAL(KIND=DP), INTENT(IN) :: K0,G0,AV,NCO,NH2

   REAL(KIND=DP) :: LAMBDA,LBAR,SCATTER,COSHIELD

!  Calculate the mean wavelength of the CO dissociating bands
   LAMBDA=LBAR(NCO,NH2)

!  Calculate the CO photodissociation rate
   RATE=K0*G0*SCATTER(AV,LAMBDA)*COSHIELD(NCO,NH2)

   RETURN
END FUNCTION CO_PHOTODISSOCIATION_RATE
!=======================================================================

!=======================================================================
!
!  CI photoionization rate taking into account grain extinction
!  and shielding by CI and H2 lines, adopting the treatment of
!  Kamp & Bertoldi (2000, A&A, 353, 276, Equation 8)
!
!-----------------------------------------------------------------------
!
!  Input parameters:
!  K0   = Unattenuated photoionization rate (in cm^3/s)
!  G0   = Incident FUV field (in Draine units)
!  AV   = visual extinction (in magnitudes)
!  KAV  = tau(λ)/tau(V) correction factor
!  NCI  = CI column density (in cm^-2)
!  NH2  = H2 column density (in cm^-2)
!  TGAS = gas temperature (in K)
!
!  Program variables:
!  RATE = CI photoionization rate taking into
!         account shielding and grain extinction
!  TAUC = optical depth in the CI absorption band
!
!-----------------------------------------------------------------------
FUNCTION CI_PHOTOIONIZATION_RATE(K0,G0,AV,KAV,NCI,NH2,TGAS) RESULT(RATE)

   USE HEALPIX_TYPES
   IMPLICIT NONE

   REAL(KIND=DP) :: RATE
   REAL(KIND=DP), INTENT(IN) :: K0,G0,AV,KAV,NCI,NH2,TGAS

   REAL(KIND=DP) :: TAUC

!  Calculate the optical depth in the CI absorption band, accounting
!  for grain extinction and shielding by CI and overlapping H2 lines
   TAUC=KAV*AV+1.1D-17*NCI+(0.9D0*TGAS**0.27D0*(NH2/1.59D21)**0.45D0)

!  Calculate the CI photoionization rate
   RATE=K0*G0*EXP(-TAUC)

   RETURN
END FUNCTION CI_PHOTOIONIZATION_RATE
!=======================================================================

!=======================================================================
!
!  SI photoionization rate -- needs to be implemented!
!  For now, use the standard expression for photorates
!
!-----------------------------------------------------------------------
!
!  Input parameters:
!  K0   = Unattenuated photoionization rate (in cm^3/s)
!  G0   = Incident FUV field (in Draine units)
!  AV   = visual extinction (in magnitudes)
!  KAV  = tau(λ)/tau(V) correction factor
!  NSI  = SI column density (in cm^-2)
!
!  Program variables:
!  RATE = SI photoionization rate taking into
!         account shielding and grain extinction
!  TAUS = optical depth in the SI absorption band
!
!-----------------------------------------------------------------------
FUNCTION SI_PHOTOIONIZATION_RATE(K0,G0,AV,KAV,NSI) RESULT(RATE)

   USE HEALPIX_TYPES
   IMPLICIT NONE

   REAL(KIND=DP) :: RATE
   REAL(KIND=DP), INTENT(IN) :: K0,G0,AV,KAV,NSI

   REAL(KIND=DP) :: TAUS

!  Calculate the optical depth in the SI absorption band, accounting
!  for grain extinction and shielding by ???
   TAUS=KAV*AV

!  Calculate the SI photoionization rate
   RATE=K0*G0*EXP(-TAUS)

   RETURN
END FUNCTION SI_PHOTOIONIZATION_RATE
!=======================================================================

!=======================================================================
!
!  H2 line self-shielding, adopting the treatment of
!  Federman, Glassgold & Kwan (1979, ApJ, 227, 466)
!
!-----------------------------------------------------------------------
!
!  Input parameters:
!  NH2  = H2 column density (in cm^-2)
!  DOPW = Doppler linewidth (in Hz)
!  RADW = radiative linewidth (in Hz)
!
!  Program variables:
!  SHIELDING_FACTOR = H2 self-shielding factor containing
!              both Doppler and radiative contributions
!  FPARA     = fraction of H2 in para state: 1/(1+o/p ratio)
!  FOSC      = oscillator strength of a typical transition
!  TAUD      = parameter tauD (eq. A7) in Federman's paper
!              (optical depth at line centre)
!  R         = parameter r  (eq. A2) in Federman's paper
!  T         = parameter t1 (eq. A6) in Federman's paper
!  U         = parameter u1 (eq. A6) in Federman's paper
!  JD        = parameter JD (eq. A8) in Federman's paper
!              (Doppler contribution to self-shielding)
!  JR        = parameter JR (eq. A9) in Federman's paper
!              (radiative contribution to self-shielding)
!
!-----------------------------------------------------------------------
FUNCTION H2SHIELD_FGK(NH2,DOPW,RADW) RESULT(SHIELDING_FACTOR)

   USE HEALPIX_TYPES
   IMPLICIT NONE

   REAL(KIND=DP) :: SHIELDING_FACTOR
   REAL(KIND=DP), INTENT(IN) :: NH2,DOPW,RADW

   REAL(KIND=DP) :: FPARA,FOSC,TAUD
   REAL(KIND=DP) :: R,T,U,JD,JR

!  Calculate the optical depth at line centre = N(H2)*f_para*(πe^2/mc)*f/(√πß) ≈ N(H2)*f_para*(1.5E-2)*f/ß
   FPARA=0.5D0 ! (assume o/p ratio=1)
   FOSC=1.0D-2 ! Typical oscillator strength
   TAUD=NH2*FPARA*(1.497358985D-2)*FOSC/DOPW

!  Calculate the Doppler core contribution to the self-shielding (JD)
   IF(TAUD.EQ.0.0D0) THEN
      JD=1.0D0
   ELSE IF(TAUD.LT.2.0D0) THEN
      JD=EXP(-(0.666666667D0*TAUD))
   ELSE IF(TAUD.LT.10.0D0) THEN
      JD=0.638D0*TAUD**(-1.25D0)
   ELSE IF(TAUD.LT.100.0D0) THEN
      JD=0.505D0*TAUD**(-1.15D0)
   ELSE
      JD=0.344D0*TAUD**(-1.0667D0)
   ENDIF

!  Calculate the radiative wing contribution to self-shielding (JR)
   IF(RADW.EQ.0.0D0) THEN
      JR=0.0D0
   ELSE
      R=RADW/(1.772453851D0*DOPW)
      T=3.02D0*((R*1.0D3)**(-0.064D0))
      U=SQRT(TAUD*R)/T
      JR=R/(T*SQRT(0.785398163D0+U**2))
   ENDIF

!  Calculate the total self-shielding factor
   SHIELDING_FACTOR=JD+JR

   RETURN
END FUNCTION H2SHIELD_FGK
!=======================================================================

!=======================================================================
!
!  H2 line shielding, using the computed values listed in 
!  Lee et al. (1996, A&A, 311, 690, Table 10)
!
!  Appropriate shielding factors are determined by performing a
!  1-dimensional spline interpolation over the values listed in
!  Table 10 of Lee et al. which include contributions from line
!  overlap, self-shielding and continuum absorption
!
!-----------------------------------------------------------------------
!
!  Input parameters:
!  NH2 = H2 column density (in cm^-2)
!
!  Program variables:
!  SHIELDING_FACTOR = total H2 shielding factor with
!              contributions from both H2 and H lines
!              and dust continuum absorption by using
!              1D spline interpolation over the grid
!  SH2_GRID  = H2 shielding factors from Lee et al. (1996)
!              as a function of H2 column density
!  SH2_DERIV = 2nd derivative of SH2_GRID values from SPLINE
!  COL_GRID  = H2 column densities (in cm^-2)
!  NUM_COL   = number of H2 column densities
!  START     = .TRUE. when H2SHIELD_LEE is first called
!
!  Functions called:
!  SPLINE = second derivative of the supplied 1D function (in spline.f90)
!  SPLINT = 1-dimensional cubic spline interpolated value (in spline.f90)
!
!-----------------------------------------------------------------------
FUNCTION H2SHIELD_LEE(NH2) RESULT(SHIELDING_FACTOR)

   USE HEALPIX_TYPES
   USE SHIELDING_MODULE, ONLY : START,NUM_COL,COL_GRID,SH2_GRID,SH2_DERIV
   IMPLICIT NONE

   REAL(KIND=DP) :: SHIELDING_FACTOR
   REAL(KIND=DP), INTENT(IN) :: NH2

   REAL(KIND=DP) :: NH2_VALUE

   IF(START) THEN
      CALL SPLINE(COL_GRID,SH2_GRID,NUM_COL,1.0D30,1.0D30,SH2_DERIV)
   ENDIF

   NH2_VALUE=NH2
   IF(NH2.LT.COL_GRID(1)) NH2_VALUE=COL_GRID(1)
   IF(NH2.GT.COL_GRID(NUM_COL)) NH2_VALUE=COL_GRID(NUM_COL)

   CALL SPLINT(COL_GRID,SH2_GRID,SH2_DERIV,NUM_COL,NH2_VALUE,SHIELDING_FACTOR)
   IF(SHIELDING_FACTOR.LT.0.0D0) SHIELDING_FACTOR=0.0D0

   RETURN
END FUNCTION H2SHIELD_LEE
!=======================================================================

!=======================================================================
!
!  CO line shielding, using the computed values listed in 
!  van Dishoeck & Black (1988, ApJ, 334, 771, Table 5)
!
!  Appropriate shielding factors are determined by performing a
!  2-dimensional spline interpolation over the values listed in
!  Table 5 of van Dishoeck & Black, which include contributions
!  from self-shielding and H2 screening
!
!-----------------------------------------------------------------------
!
!  Input parameters:
!  NCO = CO column density (in cm^-2)
!  NH2 = H2 column density (in cm^-2)
!
!  Program variables:
!  SHIELDING_FACTOR = total 12CO shielding factor with
!              contributions from both H2 and CO lines
!              by 2D spline interpolation over the grid
!  SCO_GRID  = log10 values of the 12CO shielding factors 
!              from van Dishoeck & Black (1988) as a function
!              of CO column density (1st index) and H2 column
!              density (2nd index)
!  SCO_DERIV = 2nd derivative of SCO_GRID values from SPLIE2
!  NCO_GRID  = log10 values of CO column densities (in cm^-2)
!  NH2_GRID  = log10 values of H2 column densities (in cm^-2)
!  NUM_NCO   = number of CO column densities
!  NUM_NH2   = number of H2 column densities
!  START     = .TRUE. when COSHIELD is first called
!
!  Functions called:
!  SPLIE2 = second derivative of the supplied 2D function (in spline.f90)
!  SPLIN2 = 2-dimensional cubic spline interpolated value (in spline.f90)
!
!-----------------------------------------------------------------------
FUNCTION COSHIELD(NCO,NH2) RESULT(SHIELDING_FACTOR)

   USE HEALPIX_TYPES
   USE SHIELDING_MODULE, ONLY : START,NUM_NCO,NUM_NH2,NCO_GRID,NH2_GRID,SCO_GRID,SCO_DERIV
   IMPLICIT NONE

   REAL(KIND=DP) :: SHIELDING_FACTOR
   REAL(KIND=DP), INTENT(IN) :: NCO,NH2

   REAL(KIND=DP) :: LOG_NCO,LOG_NH2

   IF(START) THEN
      CALL SPLIE2(NCO_GRID,NH2_GRID,SCO_GRID,NUM_NCO,NUM_NH2,SCO_DERIV)
      START=.FALSE.
   ENDIF

   LOG_NCO=DLOG10(NCO+1.0D0)
   LOG_NH2=DLOG10(NH2+1.0D0)

   IF(LOG_NCO.LT.NCO_GRID(1)) LOG_NCO=NCO_GRID(1)
   IF(LOG_NH2.LT.NH2_GRID(1)) LOG_NH2=NH2_GRID(1)
   IF(LOG_NCO.GT.NCO_GRID(NUM_NCO)) LOG_NCO=NCO_GRID(NUM_NCO)
   IF(LOG_NH2.GT.NH2_GRID(NUM_NH2)) LOG_NH2=NH2_GRID(NUM_NH2)

   CALL SPLIN2(NCO_GRID,NH2_GRID,SCO_GRID,SCO_DERIV,NUM_NCO,NUM_NH2,LOG_NCO,LOG_NH2,SHIELDING_FACTOR)
   SHIELDING_FACTOR=10.0D0**SHIELDING_FACTOR

   RETURN
END FUNCTION COSHIELD
!=======================================================================

!=======================================================================
!
!  Scattering by dust grains, adopting the treatment of
!  Wagenblast & Hartquist (1989, MNRAS, 237, 1019) and
!  Flannery, Roberge & Rybicki (1980, ApJ, 236, 598)
!
!-----------------------------------------------------------------------
!
!  Input parameters:
!  AV     = visual extinction (in magnitudes)
!  LAMBDA = wavelength (in Å) of incident radiation
!
!  Program variables:
!  SCATTER = attenuation factor describing the influence of
!            grain scattering on the FUV flux, dependening
!            on the total column density and wavelength of
!            light (assuming albedo=0.3 gscat=0.8)
!  TAUV    = optical depth at visual wavelength (λ=5500Å)
!  TAUL    = optical depth at wavelength LAMBDA
!  A(0)    = a(0)*exp(-k(0)*tau)
!          = relative intensity decrease for 0 < tau < 1
!  A(I)    = ∑ a(i)*exp(-k(i)*tau) for i=1,5
!            relative intensity decrease for tau ≥ 1
!  K(0)    = see A0
!  K(I)    = see A(I)
!
!  Functions called:
!  XLAMBDA = function to determine tau(λ)/tau(V)
!
!-----------------------------------------------------------------------
FUNCTION SCATTER(AV,LAMBDA)

   USE HEALPIX_TYPES
   IMPLICIT NONE

   REAL(KIND=DP) :: SCATTER
   REAL(KIND=DP), INTENT(IN) :: AV,LAMBDA

   INTEGER(KIND=I4B) :: I
   REAL(KIND=DP)     :: TAUL,TAUV,XLAMBDA,EXPONENT
   REAL(KIND=DP), DIMENSION(0:5) :: A=(/1.000D0,2.006D0,-1.438D0,0.7364D0,-0.5076D0,-0.0592D0/)
   REAL(KIND=DP), DIMENSION(0:5) :: K=(/0.7514D0,0.8490D0,1.013D0,1.282D0,2.005D0,5.832D0/)

!  Calculate the optical depth at visual wavelength
   TAUV=AV/1.086D0

!  Convert the optical depth to that at the desired wavelength
   TAUL=TAUV*XLAMBDA(LAMBDA)

!  Calculate the attenuation due to scattering by dust (SCATTER)
   SCATTER=0.0D0
   IF(TAUL.LT.1.0D0) THEN
      EXPONENT=K(0)*TAUL
      IF(EXPONENT.LT.100.0D0) THEN
         SCATTER=A(0)*EXP(-EXPONENT)
      ENDIF
   ELSE
      DO I=1,5
         EXPONENT=K(I)*TAUL
         IF(EXPONENT.LT.100.0D0) THEN
            SCATTER=SCATTER+A(I)*EXP(-EXPONENT)
         ENDIF
      ENDDO
   ENDIF

   RETURN
END FUNCTION SCATTER
!=======================================================================

!=======================================================================
!
!  Determine the ratio of the optical depth at a given wavelength to
!  that at visual wavelength (λ=5500Å) using the extinction curve of
!  Savage & Mathis (1979, ARA&A, 17, 73, Table 2)
!
!-----------------------------------------------------------------------
!
!  Input parameters:
!  LAMBDA  = wavelength (in Å)
!
!  Program variables:
!  XLAMBDA       = Ratio of tau(λ)/tau(V) at the desired wavelength
!                  by 1D spline interpolation over the grid values
!  XLAMBDA_GRID  = tau(λ)/tau(V) values, determined by dividing the
!                  Aλ/E(B-V) values from Savage & Mathis (1979) by
!                  an assumed reddening of R=AV/E(B-V)=3.1
!  XLAMBDA_DERIV = 2nd derivative of XLAMBDA_GRID values from SPLINE
!  LAMBDA_GRID   = wavelengths (in Å) listed in Table 2
!  NUM_LAMBDA    = number of wavelengths
!  START         = .TRUE. when XLAMBDA is first called
!
!  Functions called:
!  SPLINE = second derivative of the supplied 1D function (in spline.f90)
!  SPLINT = 1-dimensional cubic spline interpolated value (in spline.f90)
!
!-----------------------------------------------------------------------
FUNCTION XLAMBDA(LAMBDA)

   USE HEALPIX_TYPES
   USE SHIELDING_MODULE, ONLY : START,NUM_LAMBDA,LAMBDA_GRID,XLAMBDA_GRID,XLAMBDA_DERIV
   IMPLICIT NONE

   REAL(KIND=DP) :: XLAMBDA
   REAL(KIND=DP), INTENT(IN) :: LAMBDA

   REAL(KIND=DP) :: LAMBDA_VALUE

   IF(START) THEN
      CALL SPLINE(LAMBDA_GRID,XLAMBDA_GRID,NUM_LAMBDA,1.0D30,1.0D30,XLAMBDA_DERIV)
   ENDIF

   LAMBDA_VALUE=LAMBDA
   IF(LAMBDA.LT.LAMBDA_GRID(1)) LAMBDA_VALUE=LAMBDA_GRID(1)
   IF(LAMBDA.GT.LAMBDA_GRID(NUM_LAMBDA)) LAMBDA_VALUE=LAMBDA_GRID(NUM_LAMBDA)

   CALL SPLINT(LAMBDA_GRID,XLAMBDA_GRID,XLAMBDA_DERIV,NUM_LAMBDA,LAMBDA_VALUE,XLAMBDA)
   IF(XLAMBDA.LT.0.0D0) XLAMBDA=0.0D0

   RETURN
END FUNCTION XLAMBDA
!=======================================================================

!=======================================================================
!
!  Calculate the mean wavelength (in Å) of the 33 dissociating bands of
!  CO, weighted by their fractional contribution to the total shielding
!  from van Dishoeck & Black (1988, ApJ, 334, 771, Equation 4)
!
!-----------------------------------------------------------------------
!
!  Input parameters:
!  NCO = CO column density (in cm^-2)
!  NH2 = H2 column density (in cm^-2)
!
!  Program variables:
!  LBAR = mean wavelength (in Å)
!  U    = log10(NCO)
!  W    = log10(NH2)
!
!-----------------------------------------------------------------------
FUNCTION LBAR(NCO,NH2)

   USE HEALPIX_TYPES
   IMPLICIT NONE

   REAL(KIND=DP) :: LBAR
   REAL(KIND=DP), INTENT(IN) :: NCO,NH2

   REAL(KIND=DP) :: U,W

   U=DLOG10(NCO+1.0D0)
   W=DLOG10(NH2+1.0D0)

   LBAR=(5675.0D0 - 200.6D0*W) &
     & - (571.6D0 - 24.09D0*W)*U &
     & + (18.22D0 - 0.7664D0*W)*U**2

!  LBAR cannot be smaller than the wavelength of band 1 (913.6Å)
!  and cannot be larger than the wavelength of band 33 (1076.1Å)
   IF(LBAR.LT.913.6D0)  LBAR=913.6D0
   IF(LBAR.GT.1076.1D0) LBAR=1076.1D0

   RETURN
END FUNCTION LBAR
!=======================================================================
