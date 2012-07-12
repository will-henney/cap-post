program cubewind
  ! calculates the ram pressure of the stellar wind and of the inward HII region flow
  use wfitsutils
  implicit none
  real, dimension(:,:,:), allocatable :: u, v, w, vr, x, y, z, r, p, d, pr, pw
  character(len=128) :: prefix
  integer :: nx, ny, nz, i, j, k
  real, parameter :: MU = 1.3, MH = 1.67262158e-24, KM = 1.0e5
  real, parameter :: BOX_SIZE = 4.0, PC = 3.085677582e18
  real :: dx, pw0                    ! cell size in centimeters, ram pressure at r=dx
  ! MDOT6VW3 is product of: 
  ! + Wind mass loss rate in 1.e-6 Msun/yr
  ! + Wind velocity in 1000 km/s
  ! which for Orion is well-constrained by observations (GAHA2001 sec 5.2)
  real, parameter :: MDOT6VW3 = 0.42
  real, parameter :: VW3 = 1.3 ! 
  real, parameter :: PI = 3.1415926535897932385e0
  real, parameter :: MSUN = 1.989e33, YR = 3.15576e7

  print *, 'File prefix?'
  read '(a)', prefix

  call fitsread(trim(prefix)//'u.fits')
  nx = size(fitscube, 1)
  ny = size(fitscube, 2)
  nz = size(fitscube, 3)


  ! normalization of wind ram pressure at radius of one cell width
  dx = BOX_SIZE*PC/real(nx)
  pw0 = MDOT6VW3 * (1.e-6 * MSUN / YR) * (1.e3*KM) / (4.0*PI*dx**2)
  print '(2(a,es10.3))', "Wind pressure ", pw0, " at r = ", dx

  allocate( u(nx, ny, nz), v(nx, ny, nz), w(nx, ny, nz), &
       &    vr(nx, ny, nz), x(nx, ny, nz), y(nx, ny, nz), z(nx, ny, nz), &
       &    r(nx, ny, nz) )

  u = fitscube                  
  call fitsread(trim(prefix)//'v.fits'); v = fitscube            
  call fitsread(trim(prefix)//'w.fits'); w = fitscube                  

  forall(i=1:nx, j=1:ny, k=1:nz)
     ! cartesian distances from the center
     x(i,j,k) = real(i) - 0.5*real(nx+1)
     y(i,j,k) = real(j) - 0.5*real(ny+1)
     z(i,j,k) = real(k) - 0.5*real(nz+1)
  end forall

  ! dot product of vel with rad
  r = sqrt(x**2 + y**2 + z**2)
  vr = (u*x + v*y + w*z)/r

  deallocate(x, y, z, u, v, w)  ! individual components no longer required


  allocate( pr(nx, ny, nz), pw(nx, ny, nz),  &
       &    p(nx, ny, nz), d(nx, ny, nz) )

  call fitsread(trim(prefix)//'d.fits'); d = fitscube            
  call fitsread(trim(prefix)//'p.fits'); p = fitscube            

  ! inward radial ram pressure of HII region gas
  pr = 0.0
  where (vr < 0.0)
     ! only set where radial velocity is inward
     pr = MU*MH*d*(vr*KM)**2
  end where
  
  ! wind ram pressure
  pw = pw0/r**2

  call fitswrite(vr, trim(prefix)//'vr.fits')
  call fitswrite(pr, trim(prefix)//'pr.fits')
  call fitswrite(pw, trim(prefix)//'pw.fits')

!   open(1, file=trim(prefix)//'.vrstats', action='write')
!   write(1, '(2es12.4)') (((r(i,j,k), vr(i,j,k), i = 1, nx), j = 1, ny), k = 1, nz)
!   close(1)
  

end program cubewind
