program cubeav
  ! Calculates extinction A_V for entire cube of Garrelt simulation
  ! based on the FITS cube of density
  ! NOTE: This assumes that the source is at the center of the grid
  use wfitsutils
  use findcolumn, only: find_column_densities
  real, parameter :: sigma_dust_vband = 5.e-22
  real, parameter :: AV_per_tau = 1.086
  real, dimension(:,:,:), allocatable :: d, av
  character(len=128) :: prefix
  integer :: nx, ny, nz

  print *, 'File prefix?'
  read '(a)', prefix

  call fitsread(trim(prefix)//'d.fits')
  nx = size(fitscube, 1)
  ny = size(fitscube, 2)
  nz = size(fitscube, 3)

  allocate( d(nx, ny, nz), av(nx, ny, nz) )

  d = fitscube                  ! number density

  av = AV_per_tau*sigma_dust_vband*find_column_densities(d)

  call fitswrite(av, trim(prefix)//'av.fits')


end program cubeav
