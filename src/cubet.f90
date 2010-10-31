program cubet
  ! calculates the temperature cube from the 
  ! cubes of density, ionization fraction, pressure
  use wfitsutils
  implicit none
  real, dimension(:,:,:), allocatable :: d, x, p, t
  character(len=128) :: prefix
  integer :: nx, ny, nz
  real, parameter :: boltzmann_k = 1.38e-16

  print *, 'File prefix?'
  read '(a)', prefix

  call fitsread(trim(prefix)//'d.fits')
  nx = size(fitscube, 1)
  ny = size(fitscube, 2)
  nz = size(fitscube, 3)

  allocate( d(nx, ny, nz), x(nx, ny, nz), p(nx, ny, nz), &
       &    t(nx, ny, nz) )

  d = fitscube                  ! number density
  call fitsread(trim(prefix)//'x.fits')
  x = 1.0 - fitscube            ! ion fraction 
  call fitsread(trim(prefix)//'p.fits')
  p = fitscube                  ! pressure

  ! p = n (1+x) k T
  t = p / d / (1.+x) / boltzmann_k 
  call fitswrite(t, trim(prefix)//'t.fits')

end program cubet
