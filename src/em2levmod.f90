module em2levmod
  implicit none
  ! Routines to calculate emissivity of various lines as function of H
  ! ionization fraction (x), temperature (t), density (n).
  !
  ! For each line, one must specify parent ion, "excitation"
  ! wavelength (not necessarily the same as the line wavelength), and
  ! critical density for collisional deexcitation.
  !
  ! NOTE ON UNITS:
  ! Temperature must be fed in units of 10^4 K
  ! Density must be fed in units of cm^{-3}

  ! LIMITATIONS:
  ! 1. Ion fractions of metals are fixed functions of x. These are
  !    calibrated via Cloudy models using an ionizing spectrum from an
  !    O6 star. For harder or softer spectra, the results will be
  !    wrong!   
  ! 2. Electron fraction is assumed equal to x
  !
  ! AUTHOR:
  ! Will Henney (w.henney@astrosmo.unam.mx)
  !
  ! HISTORY: 
  ! 17 May 2007 - Major clean-up and rationalization
  !  
  ! Original version from the proplyd projects of 2001/2002. 
  !
  integer, parameter :: DP = kind(1.0d0)
  
  real(DP), private, parameter :: eV = 1.602176462e-12_dp, lambda_Ryd = 912.0_dp
  real(DP), private, parameter :: Ryd = 13.6_dp*eV, k_Bolt = 1.3806503e-16_dp
  real(DP), private, parameter :: T1 = 1.e4_dp

contains
  elemental function etaC(T,lambda)
    real(DP), intent(in) :: T,lambda
    real(DP) :: etaC
    etaC = exp( -lambda_Ryd*Ryd/k_Bolt/(T*T1)/lambda )*sqrt(1/T) / exp( -lambda_Ryd*Ryd/k_Bolt/T1/lambda )
  end function etaC
  elemental function deexc(x,T,n,ncrit) 
    real(DP), intent(in) :: x,T,n
    real(DP), intent(in) :: ncrit
    real(DP) :: deexc
    deexc = 1.0_dp/(1.0_dp+x*n/ncrit/sqrt(T))
  end function deexc
  elemental function eta_R(x,T,n) 
    real(DP), intent(in) :: x,T,n
    real(DP) :: eta_R
    eta_R = x**2 * n**2 / T
  end function eta_R
  elemental function eta_CI(x,T,n, lambda,ncrit)
    real(DP), intent(in) :: x,T,n
    real(DP), intent(in) :: lambda, ncrit
    real(DP) :: eta_CI
    eta_CI = x**2 * n**2 * etaC(T,lambda) * deexc(x,T,n,ncrit)
  end function eta_CI
  elemental function eta_Coll_Neutral(x,T,n, lambda,ncrit) result(eta)
    real(DP), intent(in) :: x,T,n
    real(DP), intent(in) :: lambda, ncrit
    real(DP) :: eta
    ! change 07 Feb 2004 to account for T-dependence of collision 
    ! strengths of a neutral line - assume linear in Temperature
    eta  = T * x*(1-x) * n**2 * etaC(T,lambda) * deexcN(x,T,n,ncrit)
  end function eta_Coll_Neutral

  ! function for a Sulfur+ line
  elemental function eta_Coll_Splus(x,T,n, lambda,ncrit) result(eta)
    real(DP), intent(in) :: x,T,n
    real(DP), intent(in) :: lambda, ncrit
    real(DP) :: eta
    eta = x*sulfur_frac(x) * n**2 * etaC(T,lambda) * deexc(x,T,n,ncrit)
  end function eta_Coll_Splus

  ! function for a Sulfur++ line
  elemental function eta_Coll_Splusplus(x,T,n, lambda,ncrit) result(eta)
    real(DP), intent(in) :: x,T,n
    real(DP), intent(in) :: lambda, ncrit
    real(DP) :: eta
    eta = x*(1.0_dp-sulfur_frac(x)) * n**2 * etaC(T,lambda) * deexc(x,T,n,ncrit)
  end function eta_Coll_Splusplus

  ! function for a Oxygen+ line
  elemental function eta_Coll_Oplus(x, T, n, lambda,ncrit) result(eta)
    real(DP), intent(in) :: x,T,n
    real(DP), intent(in) :: lambda, ncrit
    real(DP) :: eta
    ! We assume O+/O0 = H+/H0
    eta = x*x*(1.-o2plus_frac(x)) * n**2 * etaC(T,lambda) * deexc(x,T,n,ncrit)
  end function eta_Coll_Oplus

  ! function for a Oxygen++ line
  elemental function eta_Coll_Oplusplus(x, T, n, lambda,ncrit) result(eta)
    real(DP), intent(in) :: x,T,n
    real(DP), intent(in) :: lambda, ncrit
    real(DP) :: eta
    eta = x*o2plus_frac(x) * n**2 * etaC(T,lambda) * deexc(x,T,n,ncrit)
  end function eta_Coll_Oplusplus

  elemental function eta_Coll_Nplus(x, T, n, lambda, ncrit) result(eta)
    real(DP), intent(in) :: x,T,n
    real(DP), intent(in) :: lambda, ncrit
    real(DP) :: eta
    eta = x*nplus_frac(x) * n**2 * etaC(T,lambda) * deexc(x,T,n,ncrit)
  end function eta_Coll_Nplus

  elemental function sulfur_frac(x)
    real(DP), intent(in) :: x
    real(DP) :: sulfur_frac
    ! Fit to cloudy results for phi=1.e13, n=1.e4, see SiiVersusHzero in Wiki
    real(DP), parameter :: a = 6.15653_dp, b = 0.729506_dp, c = 0.629423_dp
    if (x<1.0) then
       sulfur_frac = a*(1.0_dp-x)**b/(1 + (a-1.0_dp)*(1.0_dp-x)**c)
    else
       sulfur_frac = 0.0_dp
    end if

  end function sulfur_frac

  elemental function nplus_frac(x)
    real(DP), intent(in) :: x
    real(DP) :: nplus_frac
    real(DP), parameter :: a = 1.01318_dp, chi1 = -0.0639264_dp
    real(DP), parameter :: chi2 = 2.63832_dp, w = 0.812014_dp
    nplus_frac = ion_frac(x, a, chi1, chi2, w)
  end function nplus_frac

  elemental function o2plus_frac(x)
    real(DP), intent(in) :: x
    real(DP) :: o2plus_frac
    real(DP), parameter :: a = 0.959463_dp, chi1 = 3.04291_dp
    real(DP), parameter :: chi2 = 6.3829_dp, w = 0.64178_dp
    o2plus_frac = ion_frac(x, a, chi1, chi2, w)
  end function o2plus_frac
  
  elemental function ion_frac(x, a, chi1, chi2, w)
    real(DP), intent(in) :: x, a, chi1, chi2, w
    real(DP) :: chi
    real(DP) :: ion_frac
    chi = log10(x/(1.0_dp-x))
    ion_frac = 0.25_dp*a*(1.0_dp+tanh((chi - chi1)/w))*(1.0_dp+tanh((chi2 - chi)/w))
  end function ion_frac
  
  elemental function deexcN(x,T,n,ncrit) 
    ! ditto
    real(DP), intent(in) :: x,T,n
    real(DP), intent(in) :: ncrit
    real(DP) :: deexcN
    deexcN = 1.0_dp/(1.0_dp+T*x*n/ncrit/sqrt(T))
  end function deexcN

end module em2levmod


