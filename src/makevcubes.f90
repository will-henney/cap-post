program makevcubes
  ! Created 18 Nov 2006 - Will Henney
  ! Based on makemoments, but creates the velocity cubes
  ! 
  ! The emissivity is binned according to the line-of-sight velocity
  ! to give a velocity cube. This can be sliced one way to give isovelocity images,  
  ! or sliced in the other two ways to give position velocity diagrams
  !
  ! Currently, there is no smoothing by the Doppler widths
  !
  use wfitsutils
  implicit none
  integer :: i, j, k, nx, ny, nz
  integer, parameter :: nm = 3  !Number of moments
  integer, parameter :: nd = 3  !Number of dimensions (components of vel vector)
  real, dimension(:,:,:), allocatable :: vcube
  real, dimension(:,:,:), allocatable :: e, d, p
  real, dimension(:,:,:,:), allocatable :: v
  character(len=128) :: prefix
  character(len=6) :: emtype
  real, parameter :: vmin = -40.0, vmax = +40.0
  integer, parameter :: nv = 81 ! implies dv = 1.0
  real :: sig, c0
  real, parameter :: mh = 1.67e-24


  print *, 'File prefix?'
  read '(a)', prefix
  print *, 'Emissivity type?'
  read '(a)', emtype

  ! read emissivity map
  call fitsread(trim(prefix)//'e-'//emtype//'.fits')
  ! find dimensions
  nx = size(fitscube, 1)
  ny = size(fitscube, 2)
  nz = size(fitscube, 3)
  ! allocate arrays to be read in 
  allocate( e(nx, ny, nz), d(nx, ny, nz), p(nx, ny, nz), v(nd, nx, ny, nz) )
  ! allocate new arrays to be calculated and written
  allocate(vcube(nv, ny, nz) )
  ! set emissivity
  e = fitscube

  ! read and set density (used for dust extinction)
  call fitsread(trim(prefix)//'d.fits')
  d = fitscube

  ! read and set velocities (used for dust extinction)
  call fitsread(trim(prefix)//'u.fits')
  v(1,:,:,:) = fitscube
  call fitsread(trim(prefix)//'v.fits')
  v(2,:,:,:) = fitscube
  call fitsread(trim(prefix)//'w.fits')
  v(3,:,:,:) = fitscube
  ! note, velocities are rank 4 array:
  ! - first dimension is over 3 vector components
  ! - other dimensions are over 3D position

  call fitsread(trim(prefix)//'p.fits')
  p = fitscube

  deallocate(fitscube)

  ! use a common smoothing for the entire line
  ! mean sound speed
  c0 = 1.e-5*sum(e*sqrt(p/(mh*d)))/sum(e) 
  ! sigma for the smoothing
  sig = c0/sqrt(atwt())
  print *, 'Mean sound speed, sigma are ', c0, sig


  ! Method: Rather than doing the integrations in different
  ! directions, what we do is to rotate the data cubes, so that we can
  ! always do the integration in the positive x direction. Elegant,
  ! but rather wasteful of computing resources.

  ! positive x
  call docubes
  call writecubes('xp')

  ! negative x
  call flipx
  call docubes
  call writecubes('xn')

  ! positive y
  call permute
  call docubes
  call writecubes('zn')

  ! negative y
  call flipx
  call docubes
  call writecubes('zp')

  ! positive z
  call permute
  call docubes
  call writecubes('yp')

  ! negative z
  call flipx
  call docubes
  call writecubes('yn')
  
contains 

  subroutine writecubes(view_id)
    character(len=2), intent(in) :: view_id
    call fitswrite(vcube, trim(prefix)//'vcube-'//view_id//'-'//emtype//'.fits')
  end subroutine writecubes
  

  subroutine flipx
    ! flips the current x-axis
    e(:,:,:) = e(nx:1:-1,:,:)
    d(:,:,:) = d(nx:1:-1,:,:)
    v(:,:,:,:) = v(:,nx:1:-1,:,:)
    ! also need to change sign of x-component of velocity
    v(1,:,:,:) = -v(1,:,:,:)
  end subroutine flipx

  subroutine permute
    ! shuffle the axes of the data cubes so that 
    ! new x = old y
    ! new y = old z
    ! new z = old x
    e = reshape(e, shape(e), order = (/3, 1, 2/))
    d = reshape(d, shape(d), order = (/3, 1, 2/))
    v = reshape(v, shape(v), order = (/1, 4, 2, 3/))
    ! also need to permute the components of the velocity
    v = cshift(v, 1, dim=1)     ! or should this be -1?
  end subroutine permute

  subroutine docubes
    use psfutils, only : smoothe
    integer :: mleft, mright, m1, m2, m, nm
    real :: tau, tauscale, ee
    real, dimension(:,:), allocatable :: vsmoothed, vraw
    ! This routine actually does the integration.
    ! Always in the positive x direction
    print *, 'In docubes'
    vcube = 0.0
    tauscale = (4.0*3.1e18/real(nx))*5e-22
    allocate(  vsmoothed(nv, ny), vraw(nv, ny) )
    do k = 1, nz
       do j = 1, ny
          tau = 0.0
          do i = 1, nx
             if (i > 1) tau = tau + tauscale*d(i,j,k) 
             ! First order version
             if (i > 1) then
                ! spread out the emissivity over a range of velocity cells
                mleft = 1 + nint( (v(1,i,j,k) - vmin)/(vmax - vmin)*real(nv-1) )
                mright = 1 + nint( (v(1,i-1,j,k) - vmin)/(vmax - vmin)*real(nv-1) )
                m1 = min(mleft, mright)
                m2 = max(mleft, mright)
                nm = m2 - m1 + 1 ! how many velocity cells
                ! mean emissivity per velocity cell
                ee = 0.5*(e(i,j,k) + e(i-1,j,k))/real(nm) 
                do m = m1, m2
                   if (m>1 .and. m<=nv) vcube(m,j,k) = vcube(m,j,k) +  ee * exp(-tau)
                end do
             end if
          end do
       end do
       ! Smooth this p-v plane in the v direction
       vraw = vcube(:,:,k)
       call smoothe(vraw, vsmoothed, sig, 0.0)
       vcube(:,:,k) = vsmoothed
    end do
    deallocate( vsmoothed, vraw )
  end subroutine docubes

  function atwt()
    real :: atwt
    if (emtype(1:1) == 'H') then
       atwt = 1.0
    else if (emtype(1:1) == 'N') then
       atwt = 14.0
    else if (emtype(1:1) == 'O') then
       atwt = 16.0
    else if (emtype(1:1) == 'S') then
       atwt = 32.0
    end if
  end function atwt
  

end program makevcubes
