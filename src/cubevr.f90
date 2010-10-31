program cubevr
  ! calculates the radial component of velocity for a cube
  use wfitsutils
  implicit none
  real, dimension(:,:,:), allocatable :: u, v, w, vr, x, y, z, r
  character(len=128) :: prefix
  integer :: nx, ny, nz, i, j, k


  print *, 'File prefix?'
  read '(a)', prefix

  call fitsread(trim(prefix)//'u.fits')
  nx = size(fitscube, 1)
  ny = size(fitscube, 2)
  nz = size(fitscube, 3)

  allocate( u(nx, ny, nz), v(nx, ny, nz), w(nx, ny, nz), &
       &    vr(nx, ny, nz), x(nx, ny, nz), y(nx, ny, nz), z(nx, ny, nz), &
       &    r(nx, ny, nz) )

  u = fitscube                  
  call fitsread(trim(prefix)//'v.fits')
  v = fitscube            
  call fitsread(trim(prefix)//'p.fits')
  w = fitscube                  

  forall(i=1:nx, j=1:ny, k=1:nz)
     ! cartesian distances from the center
     x(i,j,k) = real(i) - 0.5*real(nx+1)
     y(i,j,k) = real(j) - 0.5*real(ny+1)
     z(i,j,k) = real(k) - 0.5*real(nz+1)
  end forall

  ! dot product of vel with rad
  r = sqrt(x**2 + y**2 + z**2)
  vr = (u*x + v*y + w*z)/r

  call fitswrite(vr, trim(prefix)//'vr.fits')

!   open(1, file=trim(prefix)//'.vrstats', action='write')
!   write(1, '(2es12.4)') (((r(i,j,k), vr(i,j,k), i = 1, nx), j = 1, ny), k = 1, nz)
!   close(1)
  

end program cubevr
