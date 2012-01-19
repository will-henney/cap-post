program cubeemiss
  ! calculates the emissivity cube in a given line from the
  ! datacubes of density, ionization fraction, temperature
  use wfitsutils, only: fitsread, fitswrite, fitscube
  use emissmod, only: emtype, emissivity
  implicit none
  real, dimension(:,:,:), allocatable :: d, x, t, e
  character(len=128) :: prefix
  integer :: nx, ny, nz

  print *, 'File prefix?'
  read '(a)', prefix

  print *, 'Emissivity type?'
  read '(a)', emtype

  call fitsread(trim(prefix)//'d.fits')
  nx = size(fitscube, 1)
  ny = size(fitscube, 2)
  nz = size(fitscube, 3)

  allocate( d(nx, ny, nz), x(nx, ny, nz), &
       &    t(nx, ny, nz), e(nx, ny, nz) )

  d = fitscube
  call fitsread(trim(prefix)//'x.fits')
  x = 1.0 - fitscube
  call fitsread(trim(prefix)//'t.fits')
  t = fitscube

  deallocate(fitscube)          ! Free up unneeded memory

  e = emissivity(d, x, t)
  call fitswrite(e, trim(prefix)//'e-'//emtype//'.fits')

end program cubeemiss

