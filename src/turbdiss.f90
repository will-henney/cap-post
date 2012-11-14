program turbdiss
  ! calculates the turbulent dissipation rate: \div . ( \vec{v} 1/2 rho v^2)
  use wfitsutils
  implicit none
  real, dimension(:,:,:), allocatable :: u, v, w, d, x
  real, dimension(:,:,:), allocatable :: ke, tu, ptu, uke, vke, wke
  character(len=128) :: prefix
  integer :: nx, ny, nz
  real, parameter :: MU = 1.3, MH = 1.67262158e-24, KM = 1.0e5
  real, parameter :: BOX_SIZE = 4.0, PC = 3.085677582e18
  real :: dx                  ! cell size in centimeters
  real, parameter :: PI = 3.1415926535897932385e0
  real, parameter :: LSUN = 3.82e33, YR = 3.15576e7
  logical, parameter :: verbose = .True.

  print *, 'File prefix?'
  read '(a)', prefix

  call fitsread(trim(prefix)//'u.fits')
  nx = size(fitscube, 1)
  ny = size(fitscube, 2)
  nz = size(fitscube, 3)
  allocate( u(nx, ny, nz), v(nx, ny, nz), w(nx, ny, nz), d(nx, ny, nz))
  ! 2.5 GB allocated (512^3) 
  u = fitscube*KM
  call fitsread(trim(prefix)//'v.fits'); v = fitscube*KM            
  call fitsread(trim(prefix)//'w.fits'); w = fitscube*KM                  
  call fitsread(trim(prefix)//'d.fits'); d = fitscube*MU*MH
  deallocate(fitscube)
  ! 0.5 GB freed -> 2 GB

  if (verbose) print *, "Arrays u, v, w, d allocated and assigned"

  ! KE density and turbulent dissipation rate
  allocate(ke(nx, ny, nz))
  ! 0.5 GB allocated -> 2.5 GB
  ke = 0.5*d*(u**2 + v**2 + w**2) ! erg/cm^3

  if (verbose) print *, "KE assigned"
  if (verbose) print *, "RMS velocity: ", sqrt(sum(u**2 + v**2 + w**2)/(nx*ny*nz)) 
  if (verbose) print *, "Mean density: ", sum(d)/(nx*ny*nz) 
  if (verbose) print *, "Mean KE: ", sum(ke)/(nx*ny*nz) 

  deallocate(d)
  ! 0.5 GB freed -> 2 GB
  allocate( uke(nx, ny, nz), vke(nx, ny, nz), wke(nx, ny, nz))
  ! 1.5 GB allocated -> 3.5 GB
  ! find kinetic energy fluxes
  uke = ke*u                    ! erg/cm^2/s
  vke = ke*v
  wke = ke*w
  deallocate( u, v, w, ke )
  ! 2 GB freed -> 1.5 GB

  if (verbose) print *, "Vector KE flux assigned"

  allocate(tu(nx, ny, nz))
  ! 0.5 GB allocated -> 2 GB

  ! physical cell size
  dx = BOX_SIZE*PC/real(nx)

  tu = -divergence(uke, vke, wke) / dx ! erg/cm^3/s
  ! divergence allocates 4 new arrays: 2 GB allocated -> 4GB
  ! 2 GB freed -> 2 GB

  if (verbose) print *, "Divergence calculated"

  deallocate(uke, vke, wke)
  ! 1.5 GB freed -> 0.5 GB

  ! Now a strictly positive version of the same...
  ! ...which is necessary so we can sum it up
  allocate(ptu(nx, ny, nz))
  ! 0.5 GB allocated -> 1 GB
  where (tu > 0.0)
     ptu = tu
  elsewhere
     ptu = 0.0
  end where


  ! We save the full version to the FITS file, 
  ! which is negative where the flow is divergent
  call fitswrite(tu, trim(prefix)//'tu.fits')

  ! Also save separate versions for ionized and neutral gas
  ! TODO: separate the neutral into atomic and molecular
  allocate(x(nx, ny, nz))
  call fitsread(trim(prefix)//'x.fits'); x = 1.0 - fitscube
  ! Multiply by 1e18 to give numbers that ds9 can cope with 
  call fitswrite(1.e18*ptu*x, trim(prefix)//'e-Turb-i.fits')
  call fitswrite(1.e18*ptu*(1.0 - x), trim(prefix)//'e-Turb-n.fits')

  print "(a,es10.2,a)", "Dissipation rate in ionized gas: ", &
       & sum(ptu*x)*(dx/(LSUN**(1./3.)))**3, " Lsun"
  print "(a,es10.2,a)", "Dissipation rate in neutral gas: ", &
       & sum(ptu*(1.0 - x))*(dx/(LSUN**(1./3.)))**3, " Lsun"

  
contains
  function divergence(Ax, Ay, Az)
    ! Calculate divergence of the 3D vector field A = [Ax, Ay, Az]^T
    real, intent(in), dimension(:,:,:) :: Ax, Ay, Az
    real, dimension(size(Ax,1), size(Ax,2), size(Ax,3)) :: divergence
    real, dimension(size(Ax,1), size(Ax,2), size(Ax,3)) :: dAxdx, dAydy, dAzdz
    character(len=*), parameter :: div_method = "centered"
    integer :: nx, ny, nz

    nx = size(Ax,1); ny = size(Ax,2); nz = size(Ax,3)
    if (div_method == "centered") then
       ! Centered differences for the interior points
       dAxdx(2:nx-1,:,:) = 0.5*(Ax(3:nx,:,:) - Ax(1:nx-2,:,:))
       dAydy(:,2:ny-1,:) = 0.5*(Ay(:,3:ny,:) - Ay(:,1:ny-2,:))
       dAzdz(:,:,2:nz-1) = 0.5*(Az(:,:,3:nz) - Az(:,:,1:nz-2))
       ! Forward differences for the near face
       dAxdx(1,:,:) = Ax(2,:,:) - Ax(1,:,:)
       dAydy(:,1,:) = Ay(:,2,:) - Ay(:,1,:)
       dAzdz(:,:,1) = Az(:,:,2) - Az(:,:,1)
       ! Backward differences for the far face
       dAxdx(nx,:,:) = Ax(nx,:,:) - Ax(nx-1,:,:)
       dAydy(:,ny,:) = Ay(:,ny,:) - Ay(:,ny-1,:)
       dAzdz(:,:,nz) = Az(:,:,nz) - Az(:,:,nz-1)
    else
       print *, "ERROR: Method ", div_method, " not implemented yet!"
       stop
    endif

    divergence = dAxdx + dAydy + dAzdz

  end function divergence
  

end program turbdiss
