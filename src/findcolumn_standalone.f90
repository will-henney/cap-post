module findcolumn
  !
  ! Finds column densities from a given point by method of short characteristics
  !
  ! WJH 14 Dec 2008 - modified from findtau.f90
  ! WJH 23 Oct 2010 - adapted for use with Capreole3D_V3
  ! WJH 02 Jul 2012 - adapted for use with FITS files from Capreole
  !                 - NOTE: various things are hard-wired that really shouldn't be
  
  !! WJH 23 Oct 2010
  ! coldens array is now local to moduleq
  !! End of modifications for Capreole version

  !! Original use statements
!!$  use fgjwcommon, only: dp, ieq_density, ni, nj, nk, i0, j0, k0
!!$  use fabgrids, only: coldens, state => primit
!!$  use runpars, only: dx, dy, dz, idebug => ilogging
  implicit none

  integer, parameter :: dp = kind(1.0)
  integer :: idebug = 0
  real(dp), pointer, dimension(:,:,:) :: coldens
  integer :: i0, j0, k0
  integer :: ni, nj, nk
  real :: dx, dy, dz

  !! WJH 24 Oct 2010
  !! It turns out that state array ordering is different between Phab-C2 and Capreole
  !! Deal with it with an extra layer of indirection
  real(dp), pointer, dimension(:,:,:) :: density ! associated in CapreoleCompat

  real, parameter :: BOX_SIZE = 4.0 ! This might need changin for different problems
  real, parameter :: PARSEC = 3.085677582e18_dp

contains

  subroutine CapreoleCompat
    ! Most of the translation between Phab and Capreole quantities is done
    ! via renaming in the import statments, but there is a small residue that
    ! must be done here
    ni = size(density, 1)
    nj = size(density, 2)
    nk = size(density, 3)

    i0 = ni/2
    j0 = nj/2
    k0 = nk/2

    ! Here we simply assume that we know the size of the simulation box (4 pc)
    dx = BOX_SIZE*PARSEC/ni
    dy = BOX_SIZE*PARSEC/nj
    dz = BOX_SIZE*PARSEC/nk

  end subroutine CapreoleCompat
  
  subroutine single_step(i, j, k, di, dj, dk)
    integer, intent(in) :: i, j, k, di, dj, dk

    real(DP), dimension(3) :: xyz1, xyz2, xyzold, xyzp, xyzpold, xyz0

    real(DP) :: lamold, xx_lamold, yy_lamold, total_path
    real(DP) :: a, b, coldens_lamold, dens_lamold

    real(DP) :: xx1, xx0, yy1, yy0, zz1, zz0, d_length

    integer :: ijkxx, ijkyy, ijkzz
    integer :: i1, j1, k1, i2, j2, k2, iold, jold, kold
    real(DP) :: coldens11, coldens12, coldens21, coldens22
    real(DP) :: dens11, dens12, dens21, dens22

    !     # finds the column density seen from point x0, y0, z0 using
    !     # the density of neutrals

    !!!
    !!! Find the points around the current point
    !!!

    if (di > 0) then
       i2 = i
       i1 = i - 1
       iold = i - 2
    else
       i2 = i - 1
       i1 = i
       iold = i + 1
    end if
    if (dj > 0) then
       j2 = j
       j1 = j - 1
       jold = j - 2
    else
       j2 = j - 1
       j1 = j
       jold = j + 1
    end if
    if (dk > 0) then
       k2 = k
       k1 = k - 1
       kold = k - 2
    else
       k2 = k - 1
       k1 = k
       kold = k + 1
    end if
    
    ! new version using 3-element arrays
    ! x2, x1, xold are at cell corners
    ! [xyz]2 is "forward" edge of current cell
    xyz2 = (/ real(i2,DP)*dx, real(j2,DP)*dy, real(k2,DP)*dz /)       
    ! x1 is "back" edge of current cell
    xyz1 = (/ real(i1,DP)*dx, real(j1,DP)*dy, real(k1,DP)*dz /)   
    ! xold is "back" edge of previous cell along path  
    xyzold = (/ real(iold,DP)*dx, real(jold,DP)*dy, real(kold,DP)*dz /)
    ! xp, xpold are at cell centers
    xyzp = 0.5*(xyz2 + xyz1)          ! current cell
    xyzpold = 0.5*(xyz1 + xyzold)     ! previous cell

    ! position of source - this is the top-right corner of (i0,j0
    ! ,k0)th cell
    xyz0 = (/ real(i0,DP)*dx, real(j0,DP)*dy, real(k0,DP)*dz /)       
    total_path = sqrt( sum((xyzp-xyz0)**2) )

    ! Treat separately the case that we are right next to the source
    if (total_path < dx) then
       if (idebug > 2) print *, 'Next to source: ', i, j, k, total_path 
       coldens(i,j,k) = total_path * density(i,j,k)
       return
    end if

    ! Find along which dimension we are furthest from the source
    ! We will call this the zz axis
    ijkzz = maxloc( abs(xyzp - xyz0), dim=1)
    ijkyy = 1 + modulo(ijkzz+1, 3)


    ! Use similar triangles to find point where path cuts the xx-yy plane
    lamold = (xyzp(ijkzz) - xyzpold(ijkzz))/(xyzp(ijkzz) - xyz0(ijkzz))
    xx_lamold = xyzp(ijkxx) - lamold*(xyzp(ijkxx) - xyz0(ijkxx))
    yy_lamold = xyzp(ijkyy) - lamold*(xyzp(ijkyy) - xyz0(ijkyy))
    ! now find the interpolation coefficients
    a = (xx_lamold - xyzpold(ijkxx))/(xyzp(ijkxx) - xyzpold(ijkxx))
    b = (yy_lamold - xyzpold(ijkyy))/(xyzp(ijkyy) - xyzpold(ijkyy))
    
    if (idebug > 2) print '(2(3i0,tr1),3a,tr1,5(es10.2,tr1))', i, j, k, &
         & ijkxx, ijkyy, ijkzz, &
         & merge('+', '-', (/di, dj, dk/) > 0), &
         & lamold, xx_lamold, yy_lamold, a, b
    if (idebug > 2) print '(a,es10.2)', 'Total path: ', total_path
    if (a < -1.0e-6_dp .or. a - 1.0_dp > 1.0e-6_dp) &
         & print *, '!!! a out of range: ', a, '!!!'
    if (b < -1.0e-6_dp .or. b - 1.0_dp > 1.0e-6_dp) &
         & print *, '!!! b out of range: ', b, '!!!'

    select case (ijkzz)
       ! These bits can't be generalized and must be treated
       !  separately for the 3 cases
    case(1)
       ! xx = y; yy = z; zz = x
       coldens22 = coldens(i-di,j,k)
       coldens21 = coldens(i-di,j-dj,k)
       coldens12 = coldens(i-di,j,k-dk)
       coldens11 = coldens(i-di,j-dj,k-dk)
       dens22 = density(i-di,j,k)
       dens21 = density(i-di,j-dj,k)
       dens12 = density(i-di,j,k-dk)
       dens11 = density(i-di,j-dj,k-dk)
    case(2)
       ! xx = z; yy = x; zz = y
       coldens22 = coldens(i,j-dj,k)
       coldens21 = coldens(i,j-dj,k-dk)
       coldens12 = coldens(i-di,j-dj,k)
       coldens11 = coldens(i-di,j-dj,k-dk)
       dens22 = density(i,j-dj,k)
       dens21 = density(i,j-dj,k-dk)
       dens12 = density(i-di,j-dj,k)
       dens11 = density(i-di,j-dj,k-dk)
    case(3)
       ! xx = x; yy = y; zz = z
       coldens22 = coldens(i,j,k-dk) 
       coldens21 = coldens(i-di,j,k-dk)
       coldens12 = coldens(i,j-dj,k-dk) 
       coldens11 = coldens(i-di,j-dj,k-dk)
       dens22 = density(i,j,k-dk) 
       dens21 = density(i-di,j,k-dk)
       dens12 = density(i,j-dj,k-dk) 
       dens11 = density(i-di,j-dj,k-dk)
    end select

    ! check that none of the coldens that we use are negative
    if (coldens22 < 0.0_dp) coldens22 = coldens11
    if (coldens11 < 0.0_dp) coldens11 = coldens22
    if (coldens12 < 0.0_dp) coldens12 = coldens22
    if (coldens21 < 0.0_dp) coldens21 = coldens22
    ! we pray that coldens22 and coldens11 are never negative at the same time

    ! do the linear interpolation in the xx-yy plane
    coldens_lamold = a*b*coldens22 + (1.0_dp-a)*b*coldens21 + &
         & a*(1.0_dp-b)*coldens12 + (1.0_dp-a)*(1.0_dp-b)*coldens11
    dens_lamold = a*b*dens22 + (1.0_dp-a)*b*dens21 + &
         & a*(1.0_dp-b)*dens12 + (1.0_dp-a)*(1.0_dp-b)*dens11
    
    xx1 = 1.0_dp
    xx0 = a
    yy1 = 1.0_dp
    yy0 = b
    zz1 = 1.0_dp
    zz0 = 0.0_dp


    d_length = sqrt(&
         & ((xx1 - xx0)*dx)**2 + &
         & ((yy1 - yy0)*dy)**2 + &
         & ((zz1 - zz0)*dz)**2 &
         & )

    
    ! simple version - just interpolate the density between two points
    coldens(i,j,k) = coldens_lamold + d_length*0.5_dp*(dens_lamold + density(i, j, k))

    if (idebug > 2) then
       print '("final: ",5(es10.2,tr1))', &
            & coldens_lamold, dens_lamold, &
            & density(i, j, k), d_length, coldens(i,j,k)
    end if


  end subroutine single_step
  
  subroutine sweep_octant(noctant)
    integer, intent(in) :: noctant
    integer :: i, j, k
    integer :: di, dj, dk
    integer :: i1, i2, j1, j2, k1, k2
    integer :: i1p, i2p, j1p, j2p, k1p, k2p
    integer :: i1m, i2m, j1m, j2m, k1m, k2m

    i1p = max(1, i0+1)
    i2p = ni
    i1m = min(ni, i0)
    i2m = 1 

    j1p = max(1, j0+1)
    j2p = nj
    j1m = min(nj, j0)
    j2m = 1 

    k1p = max(1, k0+1)
    k2p = nk
    k1m = min(nk, k0)
    k2m = 1 
   
    select case (noctant)
    case(1) ! +++
       di = 1; dj = 1; dk = 1
       i1 = i1p; i2 = i2p; j1 = j1p; j2 = j2p; k1 = k1p; k2 = k2p
    case(2) ! -++
       di = -1; dj = 1; dk = 1
       i1 = i1m; i2 = i2m; j1 = j1p; j2 = j2p; k1 = k1p; k2 = k2p
    case(3) ! +-+
       di = 1; dj = -1; dk = 1
       i1 = i1p; i2 = i2p; j1 = j1m; j2 = j2m; k1 = k1p; k2 = k2p
    case(4) ! --+
       di = -1; dj = -1; dk = 1
       i1 = i1m; i2 = i2m; j1 = j1m; j2 = j2m; k1 = k1p; k2 = k2p
    case(5) ! ++-
       di = 1; dj = 1; dk = -1
       i1 = i1p; i2 = i2p; j1 = j1p; j2 = j2p; k1 = k1m; k2 = k2m
    case(6) ! -+-
       di = -1; dj = 1; dk = -1
       i1 = i1m; i2 = i2m; j1 = j1p; j2 = j2p; k1 = k1m; k2 = k2m
    case(7) ! +--
       di = 1; dj = -1; dk = -1
       i1 = i1p; i2 = i2p; j1 = j1m; j2 = j2m; k1 = k1m; k2 = k2m
    case(8) ! ---
       di = -1; dj = -1; dk = -1
       i1 = i1m; i2 = i2m; j1 = j1m; j2 = j2m; k1 = k1m; k2 = k2m
    end select

    do i = i1, i2, di
       do j = j1, j2, dj
          do k = k1, k2, dk
             call single_step(i, j, k, di, dj, dk)
          end do
       end do
    end do
  end subroutine sweep_octant
  
  
  function find_column_densities(input_density)
    !! WJH 02 Jul 2012: this is now a function
    !! WJH 17 Jul 2007: rewritten for 3D Cartesian
    !!
    implicit none
    real, intent(in), target, dimension(:,:,:) :: input_density
    real, target, dimension(size(input_density,1), size(input_density,2), size(input_density,3)) :: find_column_densities

    integer :: noctant

    density => input_density
    coldens => find_column_densities
    call CapreoleCompat

    ! any negative value of coldens is known not to be set, so we can
    !  avoid using it in the interpolations
    coldens(:,:,:) = -1000.0_dp

    !$OMP PARALLEL
    !$OMP DO
    do noctant = 1, 8
       call sweep_octant(noctant)
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    

  end function find_column_densities
  

end module findcolumn
