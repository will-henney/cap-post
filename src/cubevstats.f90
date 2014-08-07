program cubevstats
  ! calculates mean velocities for the data cubes
  use wfitsutils, only: fitsread, fitscube
  use findcolumn, only: find_column_densities
  implicit none
  real, parameter :: sigma_dust_vband = 5.e-22
  real, parameter :: AV_per_tau = 1.086
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: boltzmann_k = 1.3806503e-16_dp
  real(dp), parameter :: mp = 1.67262158e-24_dp, mu = 1.3_dp
  real(dp), dimension(:,:,:), allocatable :: d, xi, vr, r, x, y, z, vx, vy, vz
  logical, dimension(:,:,:), allocatable :: m
  character(len=128) :: prefix, fitsfilename
  integer :: it1, it2, it, itstep
  integer :: nx, ny, nz, i, j, k
  real(dp) :: vr_vol_tot, vr_vol_n, vr_vol_i, vr_vol_m
  real(dp) :: vr_mass_tot, vr_mass_n, vr_mass_i, vr_mass_m
  real(dp) :: vr_em_i, vr_em_if
  character(len=1), parameter :: TAB = achar(9)
  character(len=15) :: itstring
  real(dp), parameter :: pi = 3.14159265358979_dp, cubesize = 4.0_dp*3.086e18_dp
!   integer, parameter :: itsmall = 70
  real(dp) :: rmax
  ! WJH 05 Jul 2010 - New distinction between neutral/molecular gas
  ! this is copied from the python code in mhd-pressures.py
  real(dp), parameter :: mol_AV0 = 3.0_dp                           ! position of molecular transition 
  real(dp), parameter :: mol_sharpness = 4.0_dp                     ! sharpness of molecular transition
  real(dp), allocatable, dimension(:,:,:) :: AV, xmol, wi, wn, wm, xm_arg

  print *, 'Run prefix (e.g., 30112005_c)?'
  read '(a)', prefix

  print *, 'First and last time index, and step?'
  read *, it1, it2, itstep

  write(itstring,'(3("-",i4.4))') it1, it2, itstep

  open(1, file=trim(prefix)//itstring//'.vstats', action='write')

  write(1, '("# ",11(a,"'//TAB//'"))') 'Time', &
       & 'Vr_vol_t', 'Vr_vol_m', 'Vr_vol_n', 'Vr_vol_i', &
       & 'Vr_mass_t', 'Vr_mass_m', 'Vr_mass_n', 'Vr_mass_i', &
       & 'Vr_em_i', 'Vr_em_if'

  do it = it1, it2, itstep
     write(fitsfilename, '(2a,i4.4,a)') trim(prefix), '_', it, 'd.fits'
     call fitsread(trim(fitsfilename))
     if (it ==it1) then
        ! first time setup
        nx = size(fitscube, 1)
        ny = size(fitscube, 2)
        nz = size(fitscube, 3)
        allocate( d(nx, ny, nz), xi(nx, ny, nz), vr(nx, ny, nz) )
        allocate( x(nx, ny, nz), y(nx, ny, nz), z(nx, ny, nz), r(nx, ny, nz) )
        allocate( vx(nx, ny, nz), vy(nx, ny, nz), vz(nx, ny, nz) )
        allocate( m(nx, ny, nz) )
        allocate( AV(nx, ny, nz) )
        allocate( xmol(nx, ny, nz), wi(nx, ny, nz), wn(nx, ny, nz), wm(nx, ny, nz), xm_arg(nx, ny, nz) )
     end if
     d = fitscube/mp/mu

     write(fitsfilename, '(2a,i4.4,a)') trim(prefix), '_', it, 'x.fits'
     call fitsread(trim(fitsfilename))
     xi = 1.0 - fitscube

     write(fitsfilename, '(2a,i4.4,a)') trim(prefix), '_', it, 'u.fits'
     call fitsread(trim(fitsfilename))
     vx = fitscube
     write(fitsfilename, '(2a,i4.4,a)') trim(prefix), '_', it, 'v.fits'
     call fitsread(trim(fitsfilename))
     vy = fitscube
     write(fitsfilename, '(2a,i4.4,a)') trim(prefix), '_', it, 'w.fits'
     call fitsread(trim(fitsfilename))
     vz = fitscube

     write(fitsfilename, '(2a,i4.4,a)') trim(prefix), '_', it, 'av.fits'
     call fitsread(trim(fitsfilename))
     AV = fitscube

     ! ! If necessary, Calculate AV instead of reading from file
     ! AV = AV_per_tau*sigma_dust_vband*find_column_densities(real(d))

     ! Position
     forall(i=1:nx, j=1:ny, k=1:nz)
        ! cartesian distances from the center
        x(i,j,k) = real(i, dp) - 0.5_dp*real(nx+1, dp)
        y(i,j,k) = real(j, dp) - 0.5_dp*real(ny+1, dp)
        z(i,j,k) = real(k, dp) - 0.5_dp*real(nz+1, dp)
     end forall
     r = sqrt(x**2 + y**2 + z**2)
     ! radial component of velocity
     vr = (vx*x + vy*y + vz*z)/r

     ! For early times, we have some partially ionized material
     ! around the edges of the grid, which is skewing our statistics.
     ! We cut it out by only considering a sphere of radius 1 pc. 
     rmax = 0.5*sqrt(3.)*real(nx)
     print*,'max min r',maxval(r),minval(r)
     m = r < rmax
     print*,'count m',count(m)


     ! WJH 05 Jul 2010 - New distinction between neutral/molecular
     ! this is copied from the python code in mhd-pressures.py
     xm_arg = mol_sharpness*(AV-mol_AV0)
     where (xm_arg > 10.0_dp) 
        xmol = 1.0_dp
     elsewhere
        xmol = 1.0_dp - 1.0_dp/(1.0_dp + exp(xm_arg))
     end where

     ! weights for ionized/neutral/molecular
     wi = max(xi, 0.0_dp)
     wn = max((1._dp-xi)*(1._dp-xmol), 0.0_dp)
     wm = max((1._dp-xi)*xmol, 0.0_dp)

     ! volume-weighted averages
     vr_vol_tot = sum(vr)/real(nx*ny*nz, dp)
     vr_vol_m = sum(vr*wm)/sum(wm)
     vr_vol_n = sum(vr*wn)/sum(wn)
     vr_vol_i = sum(vr*wi, mask=m)/sum(wi, mask=m)

     ! mass-weighted averages
     vr_mass_tot = sum(vr*d)/sum(d)
     vr_mass_m = sum(vr*wm*d)/sum(wm*d)
     vr_mass_n = sum(vr*wn*d)/sum(wn*d)
     vr_mass_i = sum(vr*wi*d, mask=m)/sum(wi*d, mask=m)

     ! emission-weighted average (assumed prop to n^2)
     vr_em_i = sum(vr*(xi*d)**2, mask=m)/sum((xi*d)**2, mask=m)
     ! and finally, with weighting for partially ionized gas at i-front
     vr_em_if = sum(vr*xi*(1._dp-xi)*d**2, mask=m)/sum(xi*(1._dp-xi)*d**2, mask=m)

     write(1, '(i4.4,"'//TAB//'",10(es11.3,"'//TAB//'")))') it, &
             & vr_vol_tot, vr_vol_m, vr_vol_n, vr_vol_i, &
             & vr_mass_tot, vr_mass_m, vr_mass_n, vr_mass_i, &
             & vr_em_i, vr_em_if

  end do
  
end program cubevstats
