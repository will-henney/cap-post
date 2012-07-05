module psfutils
  ! Version 1.2 14 Dec 1998 :
  ! Added  another opt arg to send a user-supplied PSF (e.g. from TinyTim)
  ! Restriction: array must be 63x63
  ! Version 1.1 03 Dec 1998 :
  ! Added opt arg to smoothe for anisotropic smoothing 
  implicit none
  real, private, parameter :: pi = 3.14159265358979

  interface smoothe
     module procedure smoothe1d, smoothe2d
  end interface
  

contains

  subroutine smoothe2d(inarray,outarray,sigma, sigma2,forcepsf,triplet,tripshift)
    ! smoothe a 2D array
    ! Note that the sigma must be in *pixels*
    implicit none
    real, intent(in), dimension(:,:) :: inarray
    real, intent(out), dimension(:,:) :: outarray
    real, intent(in) :: sigma

    real, intent(in), optional :: sigma2
    real, intent(in), dimension(:,:), optional :: forcepsf
    real, intent(in), dimension(3), optional :: triplet
    integer, intent(in), dimension(3), optional :: tripshift

    real, allocatable, dimension(:,:) :: psfarray
    integer :: nx, ny, i, j, np, nc, nfx, nfy, istart, jstart
    real xp, yp, xc, yc

    if ( sigma <= 0.0 .and. .not.present(forcepsf) ) then
       ! no smoothing: just return the same array
       outarray = inarray
       return
    end if

    nx = size(inarray,1) ; ny = size(inarray,2)
!!$    np = max(nx,ny)
    np = 63
    ! find the array center
    xc = 0.5*real(np+1) ; yc = xc
    
    allocate( psfarray(np,np) )
    if (present(forcepsf)) then
       ! if we were sent a PSF array, try to use it
       if (all(shape(forcepsf)==shape(psfarray))) then
          psfarray = forcepsf
       else
          nc = (np+1)/2
          nfx = size(forcepsf,dim=1)
          nfy = size(forcepsf,dim=2)
          istart = nc-nfx/2
          jstart = nc-nfy/2
          psfarray = 0.
          psfarray(istart:istart+nfx-1,jstart:jstart+nfx-1) = forcepsf
       end if
    else
       ! otherwise, we do a Gaussian

       ! fill the array
       do j = 1, np
          do i = 1, np
             xp = real(i)-xc
             yp = real(j)-yc
             if (present(sigma2)) then
                ! anisotropic smoothing if second sigma argument present
                if (sigma2 == 0.0) then
                   ! case of no smoothing in y
                   if (j==(np+1)/2) then
                      psfarray(i,j) = exp(-(xp*xp/(sigma*sigma)))
                   else
                      psfarray(i,j) = 0.0
                   end if
                else
                   ! case of smoothing in both
                   psfarray(i,j) = exp(-(xp*xp/(sigma*sigma) +&
                        & yp*yp/(sigma2*sigma2)))
                end if
             else
                psfarray(i,j) = exp(-(xp*xp+yp*yp)/(sigma*sigma))
             end if
             
          enddo
       end do

       ! If it is a triplet state, then use all 3 components
       if (present(triplet)) then
          if (.not.present(tripshift)) then
             print *, 'Both tripshift and triplet arguments must be present'
             stop
          end if
          
          psfarray = triplet(1)*eoshift(psfarray,-tripshift(1),dim=1) +&
               & triplet(2)*eoshift(psfarray,-tripshift(2),dim=1) +&
               & triplet(3)*eoshift(psfarray,-tripshift(3),dim=1) 
       end if
       
    end if
    
    ! normalise
    psfarray = psfarray/sum(psfarray)

    ! do the smoothing
    call psf(inarray,psfarray,outarray)

    deallocate( psfarray )
  end subroutine smoothe2d

  subroutine smoothe1d(inarray,outarray,sigma, forcepsf,triplet,tripshift)
    ! smoothe a 1D array
    ! Note that the sigma must be in *pixels*
    implicit none
    real, intent(in), dimension(:) :: inarray
    real, intent(out), dimension(:) :: outarray
    real, intent(in) :: sigma

    real, intent(in), dimension(:), optional :: forcepsf
    real, intent(in), dimension(3), optional :: triplet
    integer, intent(in), dimension(3), optional :: tripshift

    real, allocatable, dimension(:,:) :: psfarray, inarray2d, outarray2d, forcepsf2d
    integer :: nx, i, np, nc, nfx, istart
    real xp, xc

    if ( sigma <= 0.0 .and. .not.present(forcepsf) ) then
       ! no smoothing: just return the same array
       outarray = inarray
       return
    end if

    nx = size(inarray,1) 
    np = 63
    ! find the array center
    xc = 0.5*real(np+1) 
    
    allocate( psfarray(np,1), inarray2d(nx,1), outarray2d(nx,1) )
    if (present(forcepsf)) then
       allocate( forcepsf2d(size(forcepsf),1) )
       forcepsf2d(:,1) = forcepsf
       ! if we were sent a PSF array, try to use it
       if (all(shape(forcepsf2d)==shape(psfarray))) then
          psfarray = forcepsf2d
       else
          nc = (np+1)/2
          nfx = size(forcepsf2d,dim=1)
          istart = nc-nfx/2
          psfarray = 0.
          psfarray(istart:istart+nfx-1,1) = forcepsf2d(:,1)
       end if
    else
       ! otherwise, we do a Gaussian

       ! fill the array
       do i = 1, np
          xp = real(i)-xc
          psfarray(i,1) = exp(-(xp*xp)/(sigma*sigma))
       end do

       ! If it is a triplet state, then use all 3 components
       if (present(triplet)) then
          if (.not.present(tripshift)) then
             print *, 'Both tripshift and triplet arguments must be present'
             stop
          end if
          
          psfarray = triplet(1)*eoshift(psfarray,-tripshift(1),dim=1) +&
               & triplet(2)*eoshift(psfarray,-tripshift(2),dim=1) +&
               & triplet(3)*eoshift(psfarray,-tripshift(3),dim=1) 
       end if
       
    end if
    
    ! normalise
    psfarray = psfarray/sum(psfarray)

    ! copy to a 2d array 
    inarray2d(:,1) = inarray

    ! do the smoothing
    call psf(inarray2d,psfarray,outarray2d)

    ! copy back to a 1d array
    outarray = outarray2d(:,1)

    deallocate( psfarray )
  end subroutine smoothe1d



  !***********************************************************************
  !     Adapted from smooth2d_fft - 28 October 1997
  !     Now will use an input array as the smoothing function 
  !     Designed to use HST psf
  !     
  !     New version (4 Aug 1994) - Uses Fast Fourier Transform
  !     The array must be smaller than 32*32*32
  !***********************************************************************
  !***********************************************************************
  subroutine psf(ain,bin,aout)
    use fft, only: rlft2
    !     TAKES AN ARRAY AIN OF DIM aNX BY aNY and convolves it with
    !     another array bin of dim bnx by bny to give the array aout
    ! !!!IMPORTANT!!!! bnx, bny should be odd
    implicit none
    !     The size of the input arrays
    real, intent(in), dimension(:,:) :: ain, bin
    real, intent(out), dimension(:,:) :: aout
    !     
    !     The maximum size of the arrays we can deal with
    integer ipmax, n, ibiggest
    !     Use twice as many points for velocity

    real, allocatable ::  dat1(:,:), dat2(:,:)
    real sx2, sy2, xfac, yfac, x, y, norm
    integer :: anx, any, bnx, bny
    integer i, j, k

    complex, allocatable :: spec1(:,:), speq1(:), spec2(:,:), speq2(:)

    anx = size(ain,1) ; any = size(ain,2)
    bnx = size(bin,1) ; bny = size(bin,2)

    ! check that bnx, bny are odd 
    if (2*(bnx/2) == bnx .or. 2*(bny/2) == bny) then
       print *, 'bnx, bny must not be even in psf'
       stop
    end if

    !     Find smallest power of 2 that is bigger than the arrays
    ibiggest = max( max(anx,any), max(bnx,bny) )
    ipmax = ceiling(log10(float(ibiggest))/log10(2.))
    !     set the data arrays to this dimension
    n = 2**ipmax


    allocate ( dat1(N,N), dat2(N,N), spec1(N/2,N), spec2(N/2,N) )
    allocate ( speq1(N), speq2(N) )


    !     First find we need to pad out the input arrays to a power of 2
    dat1 = 0.0
    dat2 = 0.0
    dat1(1:anx,1:any) = ain
    dat2(1:bnx,1:bny) = bin
    !     now put the psf center in the top right corner of the  array
    dat2 = cshift(dat2,(bnx-1)/2,dim=1)
    dat2 = cshift(dat2,(bny-1)/2,dim=2)

    !     Now FFT the two of them
!!$    call rlft3(dat1,speq1,N,N,1,1)
!!$    call rlft3(dat2,speq2,N,N,1,1)
    call rlft2(dat1,spec1,speq1,1)
    call rlft2(dat2,spec2,speq2,1)

!!$    !     map real arrays dat1, dat2 onto complex arrays spec1, spec2
!!$    spec1 = cmplx ( dat1(::2,:), dat1(2::2,:) )
!!$    spec2 = cmplx ( dat2(::2,:), dat2(2::2,:) )

    norm = 2.0/(N*N)


    !     Multiply the transforms
    spec1 = norm*spec1*spec2
    speq1 = norm*speq1*speq2

    !     map back to real array
    dat1(::2,:) = real(spec1)
    dat1(2::2,:) = aimag(spec1)

    !     Now do the inverse transform
    call rlft2(dat1,spec1,speq1,-1)

    !     Now map onto the output array
    aout = dat1(:anx,:any)


    !     deallocate arrays so there is no problem next time we use the
    !     routine 
    deallocate ( dat1, dat2, spec1, spec2 )
    deallocate ( speq1, speq2 )

  end subroutine psf




end module psfutils

