program makemap
  use wfitsutils
  use emissmod, only: extinct
  implicit none
  integer :: i, j, k, nx, ny, nz
  real, dimension(:), allocatable :: tau
  real, dimension(:,:), allocatable :: map
  real, dimension(:,:,:), allocatable :: e, d
  character(len=128) :: prefix
  character(len=6) :: emtype
  real :: extinct_sigma

  print *, 'File prefix?'
  read '(a)', prefix
  print *, 'Emissivity type?'
  read '(a)', emtype

  ! dust extinction cross-section
  extinct_sigma = extinct(emtype)

  call fitsread(trim(prefix)//'e-'//emtype//'.fits')
  nx = size(fitscube, 1)
  ny = size(fitscube, 2)
  nz = size(fitscube, 3)
  allocate( e(nx, ny, nz), d(nx, ny, nz), map(ny, nz), tau(nx) )
  e = fitscube
  call fitsread(trim(prefix)//'d.fits')
  d = fitscube
  deallocate(fitscube)

  call domap
  call fitswrite(map, trim(prefix)//'map-xp-'//emtype//'.fits')

  call flipx
  call domap
  call fitswrite(map, trim(prefix)//'map-xn-'//emtype//'.fits')

  call permute
  call domap
  call fitswrite(map, trim(prefix)//'map-zn-'//emtype//'.fits')

  call flipx
  call domap
  call fitswrite(map, trim(prefix)//'map-zp-'//emtype//'.fits')

  call permute
  call domap
  call fitswrite(map, trim(prefix)//'map-yp-'//emtype//'.fits')

  call flipx
  call domap
  call fitswrite(map, trim(prefix)//'map-yn-'//emtype//'.fits')
  
contains 

  subroutine flipx
    e(:,:,:) = e(nx:1:-1,:,:)
    d(:,:,:) = d(nx:1:-1,:,:)
  end subroutine flipx

  subroutine permute
    e = reshape(e, (/nx, ny, nz/), order = (/3, 1, 2/))
    d = reshape(d, (/nx, ny, nz/), order = (/3, 1, 2/))
  end subroutine permute

  subroutine domap
    do k = 1, nz
       do j = 1, ny
          tau(1) = 0.0
          do i = 2, nx
             tau(i) = tau(i-1) + d(i,j,k) 
          end do
          tau = tau*(4.0*3.1e18/real(nx))*extinct_sigma
          map(j, k) = sum(e(:,j,k)*exp(-tau))
       end do
    end do
  end subroutine domap

end program makemap
