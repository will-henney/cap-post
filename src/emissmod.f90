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
