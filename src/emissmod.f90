module emissmod
  use em2levmod
  implicit none

  character(len=6) :: emtype

  ! wavelength corresponding to 1 Ghz
  real(kind=DP), parameter :: lamGhz = 29.9792458_dp
  
contains

  function extinct(emtype)
    ! Returns the frequency dependent dust extinction cross section
    !
    ! All is for R = 5.5, based on data in Osterbrock & Ferland, Table
    ! 7.1
    character(len=6), intent(in) :: emtype
    real :: extinct
    real, parameter :: Vband_sigma = 5.3e-22 ! cm^2 per H nucleon
    select case (emtype)
    case('neut00', 'neut01', 'FF06cm')
       extinct = 0.0
    case('Halpha', 'S26731', 'S26716', 'N26584', 'O16300', 'S36312')
       extinct = 0.858
    case('Hb4861') 
       extinct = 1.099
    case('N25755')
       extinct = 0.97
    case('O35007') 
       extinct = 1.076
    case('O34363', 'Hg4341') 
       extinct = 1.194
    case('O23729', 'O23726') 
       extinct = 1.28
    case default
       extinct = 1.0
    end select
    extinct = extinct*Vband_sigma
  end function extinct
  
  elemental function emissivity(d, x, t) result(e)
    ! Single precision wrapper of emissity_dp
    real, intent(in) :: d, x, t
    real :: e
    e = real(emissivity_dp(real(d, DP), real(x, DP), real(t, DP)))
  end function emissivity
  


  function pah_emissivity(d, x, av) result(e)
    ! cannot be elemental since depends on the grid indices
    real, dimension(:,:,:), intent(in) :: d, x, av
    real, dimension(size(d,1),size(d,2),size(d,3)) :: e
    integer :: i, j, k
    integer :: i0, j0, k0
    real, save :: rsq, fuvflux
    !$OMP THREADPRIVATE(rsq, fuvflux)
    real, parameter :: pi = 3.14159265358979
    ! Hardwire source position to center of grid for now
    ! FIXME - this will not work for globule sims
    i0 = size(d,1)/2
    j0 = size(d,2)/2
    k0 = size(d,3)/2
    !$OMP PARALLEL
    !$OMP DO
    do k=1, size(d,1); do j=1, size(d,2); do i=1, size(d,3)
       rsq = max(1.0, real((i-i0)**2 + (j-j0)**2 + (k-k0)**2))
       ! this is the attenuation and dilution of the FUV field, as in
       ! heat_fuv line of pdr_heat in pdr.f90
       fuvflux = exp(-1.9*av(i,j,k))/(4.0*pi*rsq)
       ! PAH abundance assumed proportional to neutral fraction
       e(i,j,k) = d(i,j,k)*(1.0-x(i,j,k))*fuvflux
    end do; end do; end do
    !$OMP END DO
    !$OMP END PARALLEL
    ! To put this in real units, it must be multiplied by a factor
    !
    ! L(fuv) eta sigma / (dx)^2 / (mu m_h)
    !
    ! where
    !
    ! L(fuv) : FUV luminosity of source (erg/s)
    ! eta    : Reprocessing efficiency of PAHs (about 0.1, see Bakes,
    !          Tielens, & Bauschlischer 2001, Figure 7)
    ! sigma  : PAH cross-section (cm^2)
    ! dx     : cell size (cm)
    ! mu m_h : mean mass per H atom
    ! 
    ! Units of e, as calculated above: g/cm^3
    ! Correct units, after multiplying by the factor: erg/s/cm^3
  end function pah_emissivity


  elemental function emissivity_dp(d, x, t) result(e)
    real(DP), intent(in) :: d, x, t
    real(DP) :: e
    real(DP) :: t4, lam

    t4 = t*1.e-4_dp             ! routines in em2levmod want T like this


    ! set radio wavelength if appropriate
    select case(emtype)
    case('FF07mm')
       lam = 0.7_dp
    case('FF06cm')
       lam = 6.0_dp
    case('FF02cm')
       lam = 2.0_dp
    case('FF20cm')
       lam = 20.0_dp
    case default
       continue
    end select
    
    select case (emtype)
    case('neut00')
       ! Just to pick out the emission from the neutral gas - will hopefully 
       ! be dominated by the cold, dense stuff
       e = d*(1.0_dp-x)
    case('neut01')
       ! This is a generic FIR collisional line from the PDR
       ! Excitation temperature of 100 K and excited by electrons or neutrals
       e = d*d*(1.0_dp-x)*(x + 1.0e-3_dp)*exp(-100.0_dp/t)*sqrt(100.0_dp/t)
    case('Halpha', 'Hb4861')
       e = eta_R(x, t4, d)
    case('FF06cm', 'FF02cm', 'FF20cm')
       ! from Osterbrock & Ferland eq 4.32 - this is actually kappa...
       e = 8.24e-2*d*d*x*x*(t**(-1.35_dp))*(lam/lamGHz)**2.1_dp
       ! ...in the sense that d I / d s = kappa (B - I)
    case('O16300')
       e = eta_Coll_Neutral(x, t4, d, 6300.0_dp, 1.e6_dp)
    case('N26584')
       e = eta_Coll_Nplus(x, t4, d, 6584.0_dp, 1.e4_dp)
    case('N25755')
       e = eta_Coll_Nplus(x, t4, d, 3063.0_dp, 1.e4_dp)
    case('O35007')
       e = eta_Coll_Oplusplus(x, t4, d, 5007.0_dp, 1.0e6_dp)
    case('O34363')
       e = eta_Coll_Oplusplus(x, t4, d, 2321.0_dp, 1.0e6_dp)
    case('S26731')
       e = eta_Coll_Splus(x, t4, d, 6731.0_dp, 2489.4_dp)
    case('S26716')
       e = eta_Coll_Splus(x, t4, d, 6716.0_dp, 742.3_dp)
    case('O23279')
       e = eta_Coll_Oplus(x, t4, d, 3729.0_dp, 3400.0_dp)
    case('O23276')
       e = eta_Coll_Oplus(x, t4, d, 3726.0_dp, 1.5e4_dp)
    case('S36312')
       e = eta_Coll_Splusplus(x, t4, d, 3722.0_dp, 1.e6_dp) 
    case default
       e = -1.0_dp
    end select
  end function emissivity_dp
  
end module emissmod
