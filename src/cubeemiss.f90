program cubeemiss
  ! calculates the emissivity cube in a given line from the
  ! datacubes of density, ionization fraction, temperature
  use wfitsutils, only: fitsread, fitswrite, fitscube
  use emissmod, only: emtype, emissivity, pah_emissivity
  implicit none
  real, dimension(:,:,:), allocatable :: d, x, t, e, av
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

  allocate( d(nx, ny, nz), x(nx, ny, nz), e(nx, ny, nz) )

  d = fitscube
  call fitsread(trim(prefix)//'x.fits')
  x = 1.0 - fitscube

  if (emtype(1:3)=='PAH') then
     ! PAH emission depends on UV flux, and hence on A_V array
     call fitsread(trim(prefix)//'av.fits')
     allocate ( av(nx, ny, nz) )
     av = fitscube
     e = pah_emissivity(d, x, av)
  else
     ! others depend on temperature
     call fitsread(trim(prefix)//'t.fits')
     allocate( t(nx, ny, nz) )
     t = fitscube
     e = emissivity(d, x, t)
  end if
  deallocate(fitscube)          ! Free up unneeded memory

  call fitswrite(e, trim(prefix)//'e-'//emtype//'.fits')

end program cubeemiss

