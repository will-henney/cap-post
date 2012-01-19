module cuberotate
  ! Rotate a 3d array by arbitrary angles
  ! The "cube" need not actually be a cube - it is short for cuboid!
  !
  ! Two subroutines are provided: rotate and spin
  !
  ! For a scalar field, just call 'rotate' to transform the positions
  ! 
  ! For a vector field, call 'spin' to transform the vectors, and call
  ! 'rotate' on each component to transform the positions
  ! 
  ! Conventions:
  ! x, y, x refer to the axes of the original (input) cube
  ! xx, yy, zz refer to the axes of the transformed (output) cube
  ! theta is angle between x and xx
  ! phi is the angle between y and yy
  ! 
  ! 07 Jan 2008 - Will Henney
  implicit none
  integer :: interpolation = 2
  integer, parameter :: interpolation_nearest = 1
  integer, parameter :: interpolation_linear = 2

  real, private, parameter :: pi = 3.1415926535897932385e0, radians = pi/180.0
contains
  subroutine spin(ax, ay, az, theta, phi, axx, ayy, azz)
    ! Rotate a vector field to a new orientation (theta, phi). 
    ! Each component also needs to be rotated with 'rotate'
    implicit none
    real, intent(in), dimension(:,:,:) :: ax, ay, az
    real, intent(in) :: theta, phi
    real, intent(out), dimension(:,:,:) :: axx, ayy, azz
    real :: ct, st, cp, sp
    ! pre-calculate trig factors
    ct = cos(theta*radians)
    st = sin(theta*radians)
    cp = cos(phi*radians)
    sp = sin(phi*radians)
    axx = ax*ct + (ay*sp + az*cp)*st
    ayy = ay*cp - az*sp
    azz = -ax*st + (ay*sp + az*cp)*ct
  end subroutine spin
  
  subroutine rotate(cube, theta, phi, newcube)
    ! Rotate a cuboid grid to a new orientation (theta, phi), with
    ! optional volume clipping.  The input grid 'cube' may be a scalar
    ! field, or one component of a vector field. In the latter case,
    ! the components should also be transformed with the 'spin'
    ! subroutine.
    implicit none
    real, intent(in), dimension(:,:,:) :: cube
    real, intent(in) :: theta, phi
    real, intent(inout), allocatable, dimension(:,:,:) :: newcube
    ! Output shape is set to input shape by default (possibly clipping the data)
    ! Optionally specify different output shape by pre-allocating newcube

    integer :: nx, ny, nz, nxx, nyy, nzz, i, j, k, ii, jj, kk
    real, allocatable, dimension(:) :: xx, yy, zz
    real :: x, y, z
    ! trig factors
    real :: ct, st, cp, sp
    logical :: on_grid
    real :: a0, a1, b0, b1, c0, c1
    
    ! pre-calculate trig factors
    ct = cos(theta*radians)
    st = sin(theta*radians)
    cp = cos(phi*radians)
    sp = sin(phi*radians)

    ! size of input cube
    nx = size(cube, 1)
    ny = size(cube, 2)
    nz = size(cube, 3)

    print '(a,3(i0,tr1))', 'Input cube shape: ', shape(cube)

    if (allocated(newcube)) then
       ! Use whatever shape has been pre-allocated for newcube
       nxx = size(newcube, 1)
       nyy = size(newcube, 2)
       nzz = size(newcube, 3)
    else
       ! newcube not pre-allocated, so clip the output to the size of input cube
       nxx = nx
       nyy = ny
       nzz = nz
       allocate(newcube(nxx, nyy, nzz))
    end if
    print '(a,3(i0,tr1))', 'Output cube shape: ', nxx, nyy, nzz
    allocate(xx(nxx), yy(nyy), zz(nzz))
    ! coordinates wrt center of output cube
    xx = (/(real(i) - real(nxx)/2, i = 1, nxx)/)
    yy = (/(real(i) - real(nyy)/2, i = 1, nyy)/)
    zz = (/(real(i) - real(nzz)/2, i = 1, nzz)/)
       
    ! Loop over the output array 
    do kk = 1, nzz
       do jj = 1, nyy
          do ii = 1, nxx
             x = xx(ii)*ct - zz(kk)*st
             y = (xx(ii)*st + zz(kk)*ct)*sp + yy(jj)*cp
             z = (xx(ii)*st + zz(kk)*ct)*cp - yy(jj)*sp
             
             if (interpolation == interpolation_nearest) then
                i = nint(x + real(nx)/2)
                j = nint(y + real(ny)/2)
                k = nint(z + real(nz)/2)
                on_grid = all((/i>=1, i<=nx, j>=1, j<=ny, k>=1, k<=nz/))
                if (on_grid) then
                   newcube(ii,jj,kk) = cube(i,j,k)
                else
                   newcube(ii,jj,kk) = 0.0
                end if
             else if (interpolation == interpolation_linear) then
                i = int(x + real(nx)/2)
                j = int(y + real(ny)/2)
                k = int(z + real(nz)/2)
                on_grid = all((/i>=1, i<=nx-1, &
                     & j>=1, j<=ny-1, k>=1, k<=nz-1/))
                if (on_grid) then
                   ! Find linear interpolation coefficients
                   a1 = x - real(i) + real(nx)/2
                   a0 = 1.0 - a1
                   b1 = y - real(j) + real(ny)/2
                   b0 = 1.0 - b1
                   c1 = z - real(k) + real(nz)/2
                   c0 = 1.0 - c1
                   ! Interpolate from the 8 surrounding cells
                   newcube(ii,jj,kk) = &
                        & a0*b0*c0*cube(i,j,k) + &
                        & a0*b0*c1*cube(i,j,k+1) + &
                        & a0*b1*c0*cube(i,j+1,k) + &
                        & a0*b1*c1*cube(i,j+1,k+1) + &
                        & a1*b0*c0*cube(i+1,j,k) + &
                        & a1*b0*c1*cube(i+1,j,k+1) + &
                        & a1*b1*c0*cube(i+1,j+1,k) + &
                        & a1*b1*c1*cube(i+1,j+1,k+1)
                else
                   newcube(ii,jj,kk) = 0.0
                end if
             end if
          end do
       end do
    end do
  end subroutine rotate
  


end module cuberotate
