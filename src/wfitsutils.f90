!!!****************************************************************
!!! Various utility routines for reading 2D Fits images 
!!! Makes use of fitsio library routines 
!!!****************************************************************
!!! Author: W. Henney
!!! v 1.2: 19 Jan 2003
!!!        Modified from kfitsutils
!!! v 1.1: 09 Oct 1998
!!!        Added optional argument to fitsread to allow the filename
!!!        to be explicitly specified, thus circumventing getfilename
!!! v 1.0: 29 Sep 1998
!!!        Added getfilename routine
!!! v 0.1: 25 Sep 1998 
!!!        Based on  `balmer/fitsextract.f90' and 'codes/orion/observ/extinct.f'
module wfitsutils
  implicit none

  ! variables accessible to calling routine
  real, public, allocatable, target, dimension(:,:) :: fitsimage
  real, public, allocatable, target, dimension(:,:,:) :: fitscube
  integer, public :: nx_fits, ny_fits, nz_fits, ndim_fits
  character(200), public :: fitsfilename

  type fits_keyword
    character(len=256) :: key, value, comment
  end type fits_keyword
  
  ! variables used by fitsio routines (all private to module) 
  integer, private ::  status, inunit, readwrite, blocksize, group, outunit
  integer, private ::  bitpix, naxis, fpixel
  real, private ::  nullval
  logical, private :: anynull, simple, extend
  integer, private :: naxes(3), nfound

  interface fitswrite
     module procedure fitswrite_2d, fitswrite_3d 
  end interface
  
contains

  subroutine fitswrite_2d(outarray,outfilename, keywords)
    real, intent(in), dimension(:,:) :: outarray
    character(len=*), intent(in) :: outfilename
    type(fits_keyword), optional, intent(in), dimension(:) :: keywords
    integer :: i

    !     initialise status flag
    status = 0

    ! wipe any previous versions of the file
    call deletefile(trim(outfilename),status)

    !     choose unit number
    call ftgiou(outunit,status)
    if (status.gt.0) call printerror(status)

    ! create new FITS file
    blocksize = 1
    call ftinit(outunit,trim(outfilename),blocksize,status)
    if (status.gt.0) call printerror(status)

    ! initialize parameters
    simple = .true.
    bitpix = -32
    naxis = 2
    naxes(1) = size(outarray,1)
    naxes(2) = size(outarray,2)
    extend = .true.

    ! write required header keywords
    call ftphpr(outunit,simple,bitpix,naxis,naxes,0,1,extend,status)
    if (status.gt.0) call printerror(status)

    ! write real array to FITS file
    group = 1
    fpixel = 1
    call ftppre(outunit,group,fpixel,naxes(1)*naxes(2),outarray,status)

    ! write some more headers
    call ftpkys(outunit,'PROGRAM','wfitsutils.f90','Program that generated this file',status)
    if (status.gt.0) call printerror(status)

    if (present(keywords)) then
       do i = 1, size(keywords)
          call ftpkys(outunit,&
               & trim(keywords(i)%key),&
               & trim(keywords(i)%value),&
               & trim(keywords(i)%comment),&
               & status)
       end do
    end if
    
    call ftclos(outunit,status)
    if (status.gt.0) call printerror(status)
    call ftfiou(outunit,status)
    if (status.gt.0) call printerror(status)

  end subroutine fitswrite_2d
  
  subroutine fitswrite_3d(outarray,outfilename, keywords)
    real, intent(in), dimension(:,:,:) :: outarray
    character(len=*), intent(in) :: outfilename
    type(fits_keyword), optional, intent(in), dimension(:) :: keywords
    integer :: i

    !     initialise status flag
    status = 0

    ! wipe any previous versions of the file
    call deletefile(trim(outfilename),status)

    !     choose unit number
    call ftgiou(outunit,status)
    if (status.gt.0) call printerror(status)

    ! create new FITS file
    blocksize = 1
    call ftinit(outunit,trim(outfilename),blocksize,status)
    if (status.gt.0) call printerror(status)

    ! initialize parameters
    simple = .true.
    bitpix = -32
    naxis = 3
    naxes(1) = size(outarray,1)
    naxes(2) = size(outarray,2)
    naxes(3) = size(outarray,3)
    extend = .true.

    ! write required header keywords
    call ftphpr(outunit,simple,bitpix,naxis,naxes,0,1,extend,status)
    if (status.gt.0) call printerror(status)

    ! write real array to FITS file
    group = 1
    fpixel = 1
    call ftppre(outunit,group,fpixel,naxes(1)*naxes(2)*naxes(3),outarray,status)

    ! write some more headers
    call ftpkys(outunit,'PROGRAM','wfitsutils.f90','Program that generated this file',status)
    if (status.gt.0) call printerror(status)

    if (present(keywords)) then
       do i = 1, size(keywords)
          call ftpkys(outunit,&
               & trim(keywords(i)%key),&
               & trim(keywords(i)%value),&
               & trim(keywords(i)%comment),&
               & status)
       end do
    end if
    
    call ftclos(outunit,status)
    if (status.gt.0) call printerror(status)
    call ftfiou(outunit,status)
    if (status.gt.0) call printerror(status)

  end subroutine fitswrite_3d
  

  subroutine fitsread(forcefilename, verbose)
    implicit none
    character(*), optional :: forcefilename
    logical, optional :: verbose

    logical :: printdebug

    if (present(verbose)) then
       printdebug = verbose
    else
       printdebug = .false.
    end if

    nullval = -1.
    group = 1
    !     initialise status flag
    status = 0
    call ftgiou(inunit,status)
    if (status.gt.0) call printerror(status)
    if (printdebug) print *, 'Unit #: ', inunit

    if (present(forcefilename)) then
       if (printdebug) print *, 'forcefilename: ', trim(forcefilename)
       fitsfilename = trim(forcefilename)
       if (printdebug) print *, 'fitsfilename: ', trim(fitsfilename)
    else
       if (printdebug) print *, 'No filename specified. Will ask for one'
!!$       call getfilename
       print *, 'Sorry, getfilename has been removed'
       stop
    end if
    
    !     readonly
    readwrite = 0
    
    !     open file
!!$    write(*,*) 'Opening file'
    call ftopen(inunit,trim(fitsfilename),readwrite,blocksize,status)
    if (status.gt.0) call printerror(status)
    status = 0
    
    !     determine size of image
    nfound = 0
    naxes = 0
    call ftgknj(inunit,'NAXIS',1,3,naxes,nfound,status)
    if (status.gt.0) call printerror(status)
    status = 0
    
    if (nfound == 2) then
       if (printdebug) print *, '2D image is ', naxes(1), ' by ', naxes(2)
       naxes(3) = 1
    else if (nfound == 3) then
       if (printdebug) print *, '3D image is ', naxes(1), ' by ', naxes(2), ' by ', naxes(3)
    else
       print *, '# of axes /= 2 or 3 in fitsread'
       stop
    endif

! Don't bother with the axis scaling keywords
!!$    call ftgkns(inunit,'CTYPE',1,3,axstring,nfound,status)
!!$    print *, axstring(1:nfound)
!!$    call ftgkne(inunit,'CRVAL',1,3,axoff,nfound,status)
!!$    print *, axoff(1:nfound)
!!$    call printerror(status)
!!$    call ftgkne(inunit,'CRPIX',1,3,axpix,nfound,status)
!!$    print *, axpix(1:nfound)
!!$    call printerror(status)
!!$    call ftgkne(inunit,'CDELT',1,3,axdelt,nfound,status)
!!$    print *, axdelt(1:nfound)
!!$    call printerror(status)
    
    
    ! allocate array for image
    ndim_fits = nfound
    nx_fits = naxes(1)
    ny_fits = naxes(2)
    nz_fits = naxes(3)

    if ( ndim_fits == 2 ) then 
       if (allocated(fitsimage)) deallocate(fitsimage)
       allocate( fitsimage(nx_fits,ny_fits))
       !     read image
       call ftgpve(inunit,group,1,nx_fits*ny_fits,nullval,fitsimage,anynull,status)
    else if ( ndim_fits == 3 ) then 
       if (allocated(fitscube)) deallocate(fitscube)
       allocate( fitscube(nx_fits,ny_fits,nz_fits))
       ! read cube
       call ftgpve(inunit,group,1,nx_fits*ny_fits*nz_fits,nullval,fitscube,anynull,status)
    end if
    if (status.gt.0) call printerror(status)
    status = 0

    call ftclos(inunit,status)
    if (status.gt.0) call printerror(status)
    status = 0
    call ftfiou(inunit,status)
    if (status.gt.0) call printerror(status)
    
  end subroutine fitsread


!!!
!!! Removed 03 Jun 2004 to eliminate dependency on newlinemod
!!!
!!$  subroutine getfilename
!!$    ! radically simplified 19 Jan 2003
!!$    use newlinemod
!!$    implicit none
!!$    ! all are routines from newlinemod
!!$    print *, 'choosing dataset'
!!$    call choose_which_dataset
!!$    print *, 'reading config file'
!!$    call read_config_file
!!$    print *, 'choosing single slit'
!!$    call choose_single_slit(fitsfilename)
!!$    print *, 'found fits filename'
!!$    return
!!$  end subroutine getfilename
  

  subroutine printerror(thisstatus)
    implicit none
    integer thisstatus
    character errtext*30, errmessage*80
    
    if (thisstatus.le.0) return
    call ftgerr(thisstatus,errtext)
    
    write(*,*) 'FITSIO Error Status =', thisstatus, ': ',errtext 
    
    call ftgmsg(errmessage)
    do while (errmessage.ne.' ')
       write(*,*) errmessage
       call ftgmsg(errmessage)
    enddo
!!$    print *, real(thisstatus)/0.0
    return
  end subroutine printerror

  subroutine deletefile(filename,status)
    ! A simple little routine to delete a FITS file
    character(len=*), intent(in) :: filename
    integer, intent(inout) ::  status
    integer :: unit, blocksize

    ! Simply return if status is greater than zero
    if (status .gt. 0)return

    ! Get an unused Logical Unit Number to use to open the FITS file
    call ftgiou(unit,status)

    ! Try to open the file, to see if it exists
    call ftopen(unit,filename,1,blocksize,status)

    if (status .eq. 0)then
       !        file was opened;  so now delete it 
       call ftdelt(unit,status)
    else if (status .eq. 103)then
       !        file doesn't exist, so just reset status to zero and clear errors
       status=0
       call ftcmsg
    else
       !        there was some other error opening the file; delete the file anyway
       status=0
       call ftcmsg
       call ftdelt(unit,status)
    end if

    ! Free the unit number for later reuse
    call ftfiou(unit, status)
  end subroutine deletefile

end module wfitsutils
