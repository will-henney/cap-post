module fft
  ! routines from NR for doing FFTs
  implicit none

  !! stuff from nrtype
  INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
  INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)
  INTEGER, PARAMETER :: SP = KIND(1.0)
  INTEGER, PARAMETER :: DP = KIND(1.0D0)
  INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
  INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0))
  INTEGER, PARAMETER :: LGT = KIND(.true.)
  REAL(SP), PARAMETER :: PI=3.141592653589793238462643383279502884197_sp
  REAL(SP), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_sp
  REAL(SP), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_sp
  REAL(SP), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_sp
  REAL(SP), PARAMETER :: EULER=0.5772156649015328606065120900824024310422_sp
  REAL(DP), PARAMETER :: PI_D=3.141592653589793238462643383279502884197_dp
  REAL(DP), PARAMETER :: PIO2_D=1.57079632679489661923132169163975144209858_dp
  REAL(DP), PARAMETER :: TWOPI_D=6.283185307179586476925286766559005768394_dp

  !! FFT routines
  interface fourrow
     module procedure fourrow_sp, fourrow_dp
  end interface

  !! utility routines
  INTEGER(I4B), PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
  INTERFACE assert
     MODULE PROCEDURE assert1,assert2,assert3,assert4,assert_v
  END INTERFACE
  INTERFACE assert_eq
     MODULE PROCEDURE assert_eq2,assert_eq3,assert_eq4,assert_eqn
  END INTERFACE
  INTERFACE arth
     MODULE PROCEDURE arth_r, arth_d, arth_i
  END INTERFACE
  INTERFACE swap
     MODULE PROCEDURE swap_i,swap_r,swap_rv,swap_c, &
          swap_cv,swap_cm,swap_z,swap_zv,swap_zm, &
          masked_swap_rs,masked_swap_rv,masked_swap_rm
  END INTERFACE
  INTERFACE reallocate
     MODULE PROCEDURE reallocate_rv,reallocate_rm,&
          reallocate_iv,reallocate_im,reallocate_hv
  END INTERFACE
  INTERFACE array_copy
     MODULE PROCEDURE array_copy_r, array_copy_d, array_copy_i
  END INTERFACE

contains

!!! basic routines
  SUBROUTINE fourrow_sp(data,isign)
    IMPLICIT NONE
    COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: data
    INTEGER(I4B), INTENT(IN) :: isign
    INTEGER(I4B) :: n,i,istep,j,m,mmax,n2
    REAL(DP) :: theta
    COMPLEX(SPC), DIMENSION(size(data,1)) :: temp
    COMPLEX(DPC) :: w,wp
    COMPLEX(SPC) :: ws
    n=size(data,2)
    call assert(iand(n,n-1)==0, 'n must be a power of 2 in fourrow_sp')
    n2=n/2
    j=n2
    do i=1,n-2
       if (j > i) call swap(data(:,j+1),data(:,i+1))
       m=n2
       do
          if (m < 2 .or. j < m) exit
          j=j-m
          m=m/2
       end do
       j=j+m
    end do
    mmax=1
    do
       if (n <= mmax) exit
       istep=2*mmax
       theta=PI_D/(isign*mmax)
       wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
       w=cmplx(1.0_dp,0.0_dp,kind=dpc)
       do m=1,mmax
          ws=w
          do i=m,n,istep
             j=i+mmax
             temp=ws*data(:,j)
             data(:,j)=data(:,i)-temp
             data(:,i)=data(:,i)+temp
          end do
          w=w*wp+w
       end do
       mmax=istep
    end do
  END SUBROUTINE fourrow_sp

  SUBROUTINE fourrow_dp(data,isign)
    IMPLICIT NONE
    COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: data
    INTEGER(I4B), INTENT(IN) :: isign
    INTEGER(I4B) :: n,i,istep,j,m,mmax,n2
    REAL(DP) :: theta
    COMPLEX(DPC), DIMENSION(size(data,1)) :: temp
    COMPLEX(DPC) :: w,wp
    COMPLEX(DPC) :: ws
    n=size(data,2)
    call assert(iand(n,n-1)==0, 'n must be a power of 2 in fourrow_dp')
    n2=n/2
    j=n2
    do i=1,n-2
       if (j > i) call swap(data(:,j+1),data(:,i+1))
       m=n2
       do
          if (m < 2 .or. j < m) exit
          j=j-m
          m=m/2
       end do
       j=j+m
    end do
    mmax=1
    do
       if (n <= mmax) exit
       istep=2*mmax
       theta=PI_D/(isign*mmax)
       wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
       w=cmplx(1.0_dp,0.0_dp,kind=dpc)
       do m=1,mmax
          ws=w
          do i=m,n,istep
             j=i+mmax
             temp=ws*data(:,j)
             data(:,j)=data(:,i)-temp
             data(:,i)=data(:,i)+temp
          end do
          w=w*wp+w
       end do
       mmax=istep
    end do
  END SUBROUTINE fourrow_dp


  SUBROUTINE four2(data,isign)
    IMPLICIT NONE
    COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: data
    INTEGER(I4B), INTENT(IN) :: isign
    COMPLEX(SPC), DIMENSION(size(data,2),size(data,1)) :: temp
    call fourrow(data,isign)
    temp=transpose(data)
    call fourrow(temp,isign)
    data=transpose(temp)
  END SUBROUTINE four2

  SUBROUTINE rlft2(data,spec,speq,isign)
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: data
    COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: spec
    COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: speq
    INTEGER(I4B), INTENT(IN) :: isign
    INTEGER :: i1,j1,nn1,nn2
    REAL(DP) :: theta
    COMPLEX(SPC) :: c1=(0.5_sp,0.0_sp),c2,h1,h2,w
    COMPLEX(SPC), DIMENSION(size(data,2)-1) :: h1a,h2a
    COMPLEX(DPC) :: ww,wp
    nn1=assert_eq(size(data,1),2*size(spec,1),'rlft2: nn1')
    nn2=assert_eq(size(data,2),size(spec,2),size(speq),'rlft2: nn2')
    call assert(iand((/nn1,nn2/),(/nn1,nn2/)-1)==0, &
         'dimensions must be powers of 2 in rlft2')
    c2=cmplx(0.0_sp,-0.5_sp*isign,kind=spc)
    theta=TWOPI_D/(isign*nn1)
    wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=spc)
    if (isign == 1) then
       spec(:,:)=cmplx(data(1:nn1:2,:),data(2:nn1:2,:),kind=spc)
       call four2(spec,isign)
       speq=spec(1,:)
    end if
    h1=c1*(spec(1,1)+conjg(speq(1)))
    h1a=c1*(spec(1,2:nn2)+conjg(speq(nn2:2:-1)))
    h2=c2*(spec(1,1)-conjg(speq(1)))
    h2a=c2*(spec(1,2:nn2)-conjg(speq(nn2:2:-1)))
    spec(1,1)=h1+h2
    spec(1,2:nn2)=h1a+h2a
    speq(1)=conjg(h1-h2)
    speq(nn2:2:-1)=conjg(h1a-h2a)
    ww=cmplx(1.0_dp,0.0_dp,kind=dpc)
    do i1=2,nn1/4+1
       j1=nn1/2-i1+2
       ww=ww*wp+ww
       w=ww
       h1=c1*(spec(i1,1)+conjg(spec(j1,1)))
       h1a=c1*(spec(i1,2:nn2)+conjg(spec(j1,nn2:2:-1)))
       h2=c2*(spec(i1,1)-conjg(spec(j1,1)))
       h2a=c2*(spec(i1,2:nn2)-conjg(spec(j1,nn2:2:-1)))
       spec(i1,1)=h1+w*h2
       spec(i1,2:nn2)=h1a+w*h2a
       spec(j1,1)=conjg(h1-w*h2)
       spec(j1,nn2:2:-1)=conjg(h1a-w*h2a)
    end do
    if (isign == -1) then
       call four2(spec,isign)
       data(1:nn1:2,:)=real(spec)
       data(2:nn1:2,:)=aimag(spec)
    end if
  END SUBROUTINE rlft2


!!! Util routines from nrutil

  !BL
  SUBROUTINE assert1(n1,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, INTENT(IN) :: n1
    if (.not. n1) then
       write (*,*) 'nrerror: an assertion failed with this tag:', &
            string
       STOP 'program terminated by assert1'
    end if
  END SUBROUTINE assert1
  !BL
  SUBROUTINE assert2(n1,n2,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, INTENT(IN) :: n1,n2
    if (.not. (n1 .and. n2)) then
       write (*,*) 'nrerror: an assertion failed with this tag:', &
            string
       STOP 'program terminated by assert2'
    end if
  END SUBROUTINE assert2
  !BL
  SUBROUTINE assert3(n1,n2,n3,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, INTENT(IN) :: n1,n2,n3
    if (.not. (n1 .and. n2 .and. n3)) then
       write (*,*) 'nrerror: an assertion failed with this tag:', &
            string
       STOP 'program terminated by assert3'
    end if
  END SUBROUTINE assert3
  !BL
  SUBROUTINE assert4(n1,n2,n3,n4,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, INTENT(IN) :: n1,n2,n3,n4
    if (.not. (n1 .and. n2 .and. n3 .and. n4)) then
       write (*,*) 'nrerror: an assertion failed with this tag:', &
            string
       STOP 'program terminated by assert4'
    end if
  END SUBROUTINE assert4
  !BL
  SUBROUTINE assert_v(n,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, DIMENSION(:), INTENT(IN) :: n
    if (.not. all(n)) then
       write (*,*) 'nrerror: an assertion failed with this tag:', &
            string
       STOP 'program terminated by assert_v'
    end if
  END SUBROUTINE assert_v
  !BL
  FUNCTION assert_eq2(n1,n2,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2
    INTEGER :: assert_eq2
    if (n1 == n2) then
       assert_eq2=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq2'
    end if
  END FUNCTION assert_eq2
  !BL
  FUNCTION assert_eq3(n1,n2,n3,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2,n3
    INTEGER :: assert_eq3
    if (n1 == n2 .and. n2 == n3) then
       assert_eq3=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq3'
    end if
  END FUNCTION assert_eq3
  !BL
  FUNCTION assert_eq4(n1,n2,n3,n4,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2,n3,n4
    INTEGER :: assert_eq4
    if (n1 == n2 .and. n2 == n3 .and. n3 == n4) then
       assert_eq4=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq4'
    end if
  END FUNCTION assert_eq4
  !BL
  FUNCTION assert_eqn(nn,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, DIMENSION(:), INTENT(IN) :: nn
    INTEGER :: assert_eqn
    if (all(nn(2:) == nn(1))) then
       assert_eqn=nn(1)
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eqn'
    end if
  END FUNCTION assert_eqn
  !BL
  SUBROUTINE nrerror(string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    write (*,*) 'nrerror: ',string
    STOP 'program terminated by nrerror'
  END SUBROUTINE nrerror
  !BL
  FUNCTION arth_r(first,increment,n)
    REAL(SP), INTENT(IN) :: first,increment
    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), DIMENSION(n) :: arth_r
    INTEGER(I4B) :: k,k2
    REAL(SP) :: temp
    if (n > 0) arth_r(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth_r(k)=arth_r(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth_r(k)=arth_r(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth_r(k+1:min(k2,n))=temp+arth_r(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  END FUNCTION arth_r
  !BL
  FUNCTION arth_d(first,increment,n)
    REAL(DP), INTENT(IN) :: first,increment
    INTEGER(I4B), INTENT(IN) :: n
    REAL(DP), DIMENSION(n) :: arth_d
    INTEGER(I4B) :: k,k2
    REAL(DP) :: temp
    if (n > 0) arth_d(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth_d(k)=arth_d(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth_d(k)=arth_d(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth_d(k+1:min(k2,n))=temp+arth_d(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  END FUNCTION arth_d
  !BL
  FUNCTION arth_i(first,increment,n)
    INTEGER(I4B), INTENT(IN) :: first,increment,n
    INTEGER(I4B), DIMENSION(n) :: arth_i
    INTEGER(I4B) :: k,k2,temp
    if (n > 0) arth_i(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth_i(k)=arth_i(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth_i(k)=arth_i(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  END FUNCTION arth_i
  !BL
  SUBROUTINE swap_i(a,b)
    INTEGER(I4B), INTENT(INOUT) :: a,b
    INTEGER(I4B) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_i
  !BL
  SUBROUTINE swap_r(a,b)
    REAL(SP), INTENT(INOUT) :: a,b
    REAL(SP) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_r
  !BL
  SUBROUTINE swap_rv(a,b)
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
    REAL(SP), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_rv
  !BL
  SUBROUTINE swap_c(a,b)
    COMPLEX(SPC), INTENT(INOUT) :: a,b
    COMPLEX(SPC) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_c
  !BL
  SUBROUTINE swap_cv(a,b)
    COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: a,b
    COMPLEX(SPC), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_cv
  !BL
  SUBROUTINE swap_cm(a,b)
    COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
    COMPLEX(SPC), DIMENSION(size(a,1),size(a,2)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_cm
  !BL
  SUBROUTINE swap_z(a,b)
    COMPLEX(DPC), INTENT(INOUT) :: a,b
    COMPLEX(DPC) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_z
  !BL
  SUBROUTINE swap_zv(a,b)
    COMPLEX(DPC), DIMENSION(:), INTENT(INOUT) :: a,b
    COMPLEX(DPC), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_zv
  !BL
  SUBROUTINE swap_zm(a,b)
    COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
    COMPLEX(DPC), DIMENSION(size(a,1),size(a,2)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_zm
  !BL
  SUBROUTINE masked_swap_rs(a,b,mask)
    REAL(SP), INTENT(INOUT) :: a,b
    LOGICAL(LGT), INTENT(IN) :: mask
    REAL(SP) :: swp
    if (mask) then
       swp=a
       a=b
       b=swp
    end if
  END SUBROUTINE masked_swap_rs
  !BL
  SUBROUTINE masked_swap_rv(a,b,mask)
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
    LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
    REAL(SP), DIMENSION(size(a)) :: swp
    where (mask)
       swp=a
       a=b
       b=swp
    end where
  END SUBROUTINE masked_swap_rv
  !BL
  SUBROUTINE masked_swap_rm(a,b,mask)
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a,b
    LOGICAL(LGT), DIMENSION(:,:), INTENT(IN) :: mask
    REAL(SP), DIMENSION(size(a,1),size(a,2)) :: swp
    where (mask)
       swp=a
       a=b
       b=swp
    end where
  END SUBROUTINE masked_swap_rm
  !BL
  FUNCTION zroots_unity(n,nn)
    INTEGER(I4B), INTENT(IN) :: n,nn
    COMPLEX(SPC), DIMENSION(nn) :: zroots_unity
    INTEGER(I4B) :: k
    REAL(SP) :: theta
    zroots_unity(1)=1.0
    theta=TWOPI/n
    k=1
    do
       if (k >= nn) exit
       zroots_unity(k+1)=cmplx(cos(k*theta),sin(k*theta),SPC)
       zroots_unity(k+2:min(2*k,nn))=zroots_unity(k+1)*&
            zroots_unity(2:min(k,nn-k))
       k=2*k
    end do
  END FUNCTION zroots_unity

  FUNCTION reallocate_rv(p,n)
    REAL(SP), DIMENSION(:), POINTER :: p, reallocate_rv
    INTEGER(I4B), INTENT(IN) :: n
    INTEGER(I4B) :: nold,ierr
    allocate(reallocate_rv(n),stat=ierr)
    if (ierr /= 0) call &
         nrerror('reallocate_rv: problem in attempt to allocate memory')
    if (.not. associated(p)) RETURN
    nold=size(p)
    reallocate_rv(1:min(nold,n))=p(1:min(nold,n))
    deallocate(p)
  END FUNCTION reallocate_rv
  !BL
  FUNCTION reallocate_iv(p,n)
    INTEGER(I4B), DIMENSION(:), POINTER :: p, reallocate_iv
    INTEGER(I4B), INTENT(IN) :: n
    INTEGER(I4B) :: nold,ierr
    allocate(reallocate_iv(n),stat=ierr)
    if (ierr /= 0) call &
         nrerror('reallocate_iv: problem in attempt to allocate memory')
    if (.not. associated(p)) RETURN
    nold=size(p)
    reallocate_iv(1:min(nold,n))=p(1:min(nold,n))
    deallocate(p)
  END FUNCTION reallocate_iv
  !BL
  FUNCTION reallocate_hv(p,n)
    CHARACTER(1), DIMENSION(:), POINTER :: p, reallocate_hv
    INTEGER(I4B), INTENT(IN) :: n
    INTEGER(I4B) :: nold,ierr
    allocate(reallocate_hv(n),stat=ierr)
    if (ierr /= 0) call &
         nrerror('reallocate_hv: problem in attempt to allocate memory')
    if (.not. associated(p)) RETURN
    nold=size(p)
    reallocate_hv(1:min(nold,n))=p(1:min(nold,n))
    deallocate(p)
  END FUNCTION reallocate_hv
  !BL
  FUNCTION reallocate_rm(p,n,m)
    REAL(SP), DIMENSION(:,:), POINTER :: p, reallocate_rm
    INTEGER(I4B), INTENT(IN) :: n,m
    INTEGER(I4B) :: nold,mold,ierr
    allocate(reallocate_rm(n,m),stat=ierr)
    if (ierr /= 0) call &
         nrerror('reallocate_rm: problem in attempt to allocate memory')
    if (.not. associated(p)) RETURN
    nold=size(p,1)
    mold=size(p,2)
    reallocate_rm(1:min(nold,n),1:min(mold,m))=&
         p(1:min(nold,n),1:min(mold,m))
    deallocate(p)
  END FUNCTION reallocate_rm
  !BL
  FUNCTION reallocate_im(p,n,m)
    INTEGER(I4B), DIMENSION(:,:), POINTER :: p, reallocate_im
    INTEGER(I4B), INTENT(IN) :: n,m
    INTEGER(I4B) :: nold,mold,ierr
    allocate(reallocate_im(n,m),stat=ierr)
    if (ierr /= 0) call &
         nrerror('reallocate_im: problem in attempt to allocate memory')
    if (.not. associated(p)) RETURN
    nold=size(p,1)
    mold=size(p,2)
    reallocate_im(1:min(nold,n),1:min(mold,m))=&
         p(1:min(nold,n),1:min(mold,m))
    deallocate(p)
  END FUNCTION reallocate_im
  !BL
!BL
  SUBROUTINE array_copy_r(src,dest,n_copied,n_not_copied)
    REAL(SP), DIMENSION(:), INTENT(IN) :: src
    REAL(SP), DIMENSION(:), INTENT(OUT) :: dest
    INTEGER(I4B), INTENT(OUT) :: n_copied, n_not_copied
    n_copied=min(size(src),size(dest))
    n_not_copied=size(src)-n_copied
    dest(1:n_copied)=src(1:n_copied)
  END SUBROUTINE array_copy_r
!BL
  SUBROUTINE array_copy_d(src,dest,n_copied,n_not_copied)
    REAL(DP), DIMENSION(:), INTENT(IN) :: src
    REAL(DP), DIMENSION(:), INTENT(OUT) :: dest
    INTEGER(I4B), INTENT(OUT) :: n_copied, n_not_copied
    n_copied=min(size(src),size(dest))
    n_not_copied=size(src)-n_copied
    dest(1:n_copied)=src(1:n_copied)
  END SUBROUTINE array_copy_d
!BL
  SUBROUTINE array_copy_i(src,dest,n_copied,n_not_copied)
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: src
    INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: dest
    INTEGER(I4B), INTENT(OUT) :: n_copied, n_not_copied
    n_copied=min(size(src),size(dest))
    n_not_copied=size(src)-n_copied
    dest(1:n_copied)=src(1:n_copied)
  END SUBROUTINE array_copy_i
!BL


end module fft


