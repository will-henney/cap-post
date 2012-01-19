program makerotmap
  use emissmod, only: extinct
  use wfitsutils, only: fitswrite
  use cuberotate, only: rotate, &
       & interpolation, interpolation_nearest, interpolation_linear
  implicit none
  integer :: i, j, k, nx=0, ny, nz, nxx, nyy, nzz
  real, dimension(:), allocatable :: tau
  real, dimension(:,:), allocatable :: map
  real, dimension(:,:,:), allocatable :: e, d ! original datacubes
  real :: theta, phi
  real, dimension(:,:,:), allocatable :: ee, dd ! rotated datacubes
  character(len=128) :: prefix
  character(len=6) :: emtype
  integer :: itime
  character(len=4) :: ctime
  character(len=13) :: rotid
  real :: extinct_sigma, zsize, dzz
  integer :: shrink = 4

  ! options are nearest/linear
  interpolation = interpolation_linear

  print *, 'File prefix?'
  read '(a)', prefix

  print *, 'Save time?'
  read *, itime
  write(ctime, '(i4.4)') itime
  
  print *, 'Emissivity type?'
  read '(a)', emtype

  ! dust extinction cross-section
  extinct_sigma = extinct(emtype)

  print *, 'Cube extent in z-direction (parsec)?'
  read *, zsize

  print *, 'Rotation theta, phi?'
  read *, theta, phi
  write(rotid, '(a,2(sp,i4.3),a)') '-rot', int(theta), int(phi), '-'

  call read_cube(d, 'd')
  print '("Original grid size: ",i0,"x",i0,"x",i0)', shape(d)
  ! Garrelt cubes are already in cm^-3 !!!
!   d = d/(1.3*1.67262158e-24)    ! convert to cm^-3
  dzz = zsize*3.085677582e18/real(nz) ! physical pixel scale
  
  nxx = size(d, 1)/shrink
  nyy = size(d, 1)/shrink
  nzz = size(d, 1)/shrink
  allocate(dd(nxx,nyy,nzz))
  call rotate(d, theta, phi, dd)
  deallocate(d)
  print '("Rotated grid size: ",i0,"x",i0,"x",i0)', nxx, nyy, nzz

  call read_cube(e, 'e-'//emtype)
  allocate(ee(nxx,nyy,nzz))
  call rotate(e, theta, phi, ee)
  deallocate(e)

  ! integration is along the new zz-axis
  allocate(map(nxx, nyy), tau(nzz) )

  call domap

  call fitswrite(map, trim(prefix)//'map'//rotid//emtype//ctime//'.fits')

contains

  subroutine read_cube(var, id)
    use wfitsutils, only: fitsread, fitscube
    real, intent(inout), dimension(:,:,:), allocatable :: var
    character(len=*), intent(in) :: id
    call fitsread(trim(prefix)//'_'//ctime//id//'.fits', verbose=.true.)
    if (nx==0) then
       nx = size(fitscube, 1)
       ny = size(fitscube, 2)
       nz = size(fitscube, 3)
    end if
    allocate(var(nx,ny,nz))
    var = fitscube
    deallocate(fitscube)
  end subroutine read_cube

  subroutine domap
    do i = 1, nxx
       do j = 1, nyy
          tau(1) = 0.0
          do k = 2, nzz
             tau(k) = tau(k-1) + dd(i,j,k) 
          end do
          tau = tau*dzz*extinct_sigma
          map(i, j) = sum(ee(i,j,:)*exp(-tau))
       end do
    end do
  end subroutine domap

end program makerotmap
