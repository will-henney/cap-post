program cubestats
  ! calculates various statistics for the data cubes
  ! 17 Jul 2008 - minimal modifications to work with the Fabio datasets
  use wfitsutils, only: fitsread, fitscube
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: boltzmann_k = 1.3806503e-16, mp = 1.67262158e-24, mu = 1.3
  real(dp), dimension(:,:,:), allocatable :: d, x, u, v, w, r
  logical, dimension(:,:,:), allocatable :: ion_mask, m
  character(len=128) :: prefix, fitsfilename
  integer :: it1, it2, it, itstep
  integer :: nx, ny, nz, i, j, k
  real(dp) :: rmax
  real(dp) :: rif_max, rif_min
  real(dp) :: vrms_vol_t, vrms_vol_n, vrms_vol_i, sumx, sumy, sumz
  real(dp) :: vrms_mass_t, vrms_mass_n, vrms_mass_i
  real(dp) :: rx1, rx2, rmean_vol_i, rmean_mass_i
  real(dp) :: frac_ion_vol1, frac_ion_vol2, frac_ion_mass
  character(len=1), parameter :: TAB = achar(9)
  character(len=15) :: itstring
  real(dp), parameter :: pi = 3.14159265358979, cubesize = 4.0*3.086e18

  print *, 'Run prefix (e.g., 30112005_c)?'
  read '(a)', prefix


  print *, 'First and last time index, and step?'
!  read *, it1, it2, itstep
  read *, it

!  write(itstring,'(3("-",i4.4))') it1, it2, itstep
  write(itstring,'("_",i4.4)') it

  open(1, file=trim(prefix)//trim(itstring)//'.stats', action='write')
  open(2, file=trim(prefix)//trim(itstring)//'.rstats', action='write')

  write(1, '("# ",10(a,"'//TAB//'"))') 'Time', &
       & 'Ifrac_v', 'Ifrac_v2', 'Ifrac_m', &
       & 'Vrms_vol_t', 'Vrms_vol_n', 'Vrms_vol_i', &
       & 'Vrms_mass_t', 'Vrms_mass_n', 'Vrms_mass_i'
  write(2, '("# ",7(a,"'//TAB//'"))') 'Time', &
       & 'rx1', 'rx2', 'rmean_vol_i', 'rmean_mass_i', 'rif_min', 'rif_max'

!  do it = it1, it2, itstep
     write(fitsfilename, '(2a,i4.4,a)') trim(prefix), '_', it, 'd.fits'
     print*,'fits filename',trim(fitsfilename)
     call fitsread(trim(fitsfilename))
!     if (it ==it1) then
        ! first time setup
        nx = size(fitscube, 1)
        ny = size(fitscube, 2)
        nz = size(fitscube, 3)
        allocate( d(nx, ny, nz), x(nx, ny, nz), ion_mask(nx, ny, nz) )
        allocate( u(nx, ny, nz), v(nx, ny, nz), w(nx, ny, nz) )
        allocate( r(nx, ny, nz), m(nx, ny, nz) )
!     end if

     print*,'sizes ',nx,ny,nz
     d = fitscube/mp/mu

     print*,'done density'

     write(fitsfilename, '(2a,i4.4,a)') trim(prefix), '_', it, 'x.fits'
     call fitsread(trim(fitsfilename))
     x = fitscube
     x = 1.0 - x
     print*,'done ionized fraction', minval(x), maxval(x)

     write(fitsfilename, '(2a,i4.4,a)') trim(prefix), '_', it, 'u.fits'
     call fitsread(trim(fitsfilename))
     u = fitscube
     print*,'done u velocity'

     write(fitsfilename, '(2a,i4.4,a)') trim(prefix), '_', it, 'v.fits'
     call fitsread(trim(fitsfilename))
     v = fitscube
     print*,'done v velocity'

     write(fitsfilename, '(2a,i4.4,a)') trim(prefix), '_', it, 'w.fits'
     call fitsread(trim(fitsfilename))
     w = fitscube
     print*,'done w velocity'

     forall(i=1:nx, j=1:ny, k=1:nz)
        ! radial distance from center in grid units
        r(i, j, k) = &
             & sqrt( &
             &        (real(i) - 0.5*real(nx+1))**2 &
             &      + (real(j) - 0.5*real(ny+1))**2 &
             &      + (real(k) - 0.5*real(nz+1))**2 &
             &     )
     end forall
     print*,'done r'

     ! For early times, we have some partially ionized material
     ! around the edges of the grid, which is skewing our statistics.
     ! We cut it out by only considering a sphere of radius 1 pc at start,
     ! growing linearly with time. 
!     rmax = (1.0 + real(it)/50.0)*0.25*real(nx)
     rmax = 0.5*sqrt(3.)*real(nx)
     print*,'max min r',maxval(r),minval(r)
     m = r < rmax
     print*,'count m',count(m)

!     ion_mask = x>0.5
     ion_mask = x>0.9
     print*,'count ion_mask',count(ion_mask)

     frac_ion_vol1 = count(ion_mask)/real(nx*ny*nz)
     print*,'done frac_ion_vol1',frac_ion_vol1
     ! Note that we do not apply the radius cutoff mask to the denominator
     ! since we are considering that the stuff at the edges shoud "really"
     ! be neutral
     frac_ion_vol2 = sum(x, mask=m)/(real(nx)*real(ny)*real(nz))
     print*,'done frac_ion_vol2',frac_ion_vol2,sum(x)/real(nx*ny*nz), &
     & sum(x,mask=r<500.)/real(nx*ny*nz),sum(x,mask=x>=0.9)/real(nx*ny*nz)
     sumx = 0.0
     sumy = 0.0
     sumz = 0.0
     do k = 1,nz
        do j = 1,ny
           do i = 1, nx
              if(x(i,j,k)> 0.9)sumx = sumx + 1.0
              if(r(i,j,k)< rmax)sumy = sumy + x(i,j,k)
              sumz = sumz + 1.0
           enddo
        enddo
     enddo
     ! sumx = sumx/(real(nx)*real(ny)*real(nz))
     ! sumy = sumy/real(nx*ny*nz)
     ! sumz = sumz/real(nx*ny*nz)
     print*,'done sumx, sumy, sumz ', sumx, sumy, sumz
     print *, 'NX x NY x NZ = ', real(nx*ny*nz)

     frac_ion_mass = sum(x*d, mask=m)/sum(d)
     print*,'done frac_ion_mass'

     ! direction-averaged radius of ionized volume
     rx1 = cubesize*(frac_ion_vol1*3.0/4.0/pi)**(1./3.)
     rx2 = cubesize*(frac_ion_vol2*3.0/4.0/pi)**(1./3.)
     ! direction-averaged radius of ionized volume
     rx1 = cubesize*(frac_ion_vol1*3.0/4.0/pi)**(1./3.)
     rx2 = cubesize*(frac_ion_vol2*3.0/4.0/pi)**(1./3.)
     print*,'done rx1, rx2',rx1,rx2

     ! min/max i-front radius
     rif_max = maxval(r, mask=ion_mask)
     rif_min = minval(r, mask=.not.ion_mask)
     print*,'done rif_max, rif_min',rif_max,rif_min

     ! mean radius of ionized gas
     print*,'bits of rmean',sum(x,mask=m),sum(d*x,mask=m),sum(r*d*x,mask=m)
     rmean_mass_i = (cubesize/real(nx))/sum(d*x, mask=m)
     print*,'done rmean_mass_i',rmean_mass_i
     rmean_mass_i = rmean_mass_i*sum(r*d*x, mask=m)
     print*,'done rmean_mass_i -2',rmean_mass_i
     rmean_vol_i = (cubesize/real(nx))*sum(r*x, mask=m)/sum(x, mask=m)
     print*,'done rmean_vol_i',rmean_vol_i

     ! the 1D RMS velocity - this is the one we will use
     vrms_vol_t = sqrt(sum(u*u + v*v + w*w)/real(3*nx*ny*nz))
     ! same but mass-weighted
     vrms_mass_t = sqrt(sum(d*(u*u + v*v + w*w))/(3*sum(d)))
     print*,'done vrms_vol_t, vrms_mass_t',vrms_vol_t,vrms_mass_t

     vrms_vol_i = sqrt(sum((u*u + v*v + w*w)*x, mask=m)/(3*sum(x, mask=m)))
     vrms_mass_i = sqrt(sum((u*u + v*v + w*w)*d*x, mask=m)/(3*sum(d*x, mask=m)))
     print*,'done vrms_vol_i, vrms_mass_i',vrms_vol_i,vrms_mass_i

     vrms_vol_n = sqrt(sum((u*u + v*v + w*w)*(1.-x))/(3*sum((1.-x))))
     vrms_mass_n = sqrt(sum((u*u + v*v + w*w)*d*(1.-x))/(3*sum(d*(1.-x))))
     print*,'done vrms_vol_n, vrms_mass_n',vrms_vol_n,vrms_mass_n

!     if (mod(it,10)==0) print *, 'Done timestep: ', it
     print *, 'Done timestep: ', it
     write(1, '(i4.4,"'//TAB//'",9(es11.3,"'//TAB//'"))') it, &
          & frac_ion_vol1, frac_ion_vol2, frac_ion_mass, &
          & vrms_vol_t, vrms_vol_n, vrms_vol_i,  &
          & vrms_mass_t, vrms_mass_n, vrms_mass_i


     write(2, '(i4.4,"'//TAB//'",6(es11.3,"'//TAB//'"))') it, &
          & rx1, rx2, rmean_vol_i, rmean_mass_i, rif_min, rif_max
!  end do

end program cubestats
