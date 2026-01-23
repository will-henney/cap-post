program cubecdn
  ! Calculates neutral column density to each point on grid
  ! NOTE: This assumes that the source is at the center of the grid
  use wfitsutils
  use findcolumn, only: find_column_densities
  real, dimension(:,:,:), allocatable :: d, xn, cdn
  character(len=128) :: prefix
  integer :: nx, ny, nz

  print *, 'File prefix?'
  read '(a)', prefix

  call fitsread(trim(prefix)//'d.fits')
  nx = size(fitscube, 1)
  ny = size(fitscube, 2)
  nz = size(fitscube, 3)

  allocate( d(nx, ny, nz), cdn(nx, ny, nz), xn(nx, ny, nz) )

  d = fitscube                  ! number density

  ! neutral fraction
  call fitsread(trim(prefix)//'x.fits'); xn = fitscube
  
  cdn = find_column_densities(d * xn)

  call fitswrite(cdn, trim(prefix)//'cdn.fits')


end program cubecdn
