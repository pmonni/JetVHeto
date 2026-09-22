! $Id: interpolation.f90 138 2003-01-07 19:17:32Z gsalam $
module interpolation
  use types
  implicit none
  private

  public :: interpolate, interpolate_uniform, polint
  public :: uniform_interpolation_weights, spline, ispline
contains

  !----------------------------------------------------------------------
  ! In an ordered table of xarr, find the interpolation which uses m
  ! points
  subroutine interpolate_uniform(xarr,yarr, x, y, dy, m)
    real(dp), intent(in)  :: xarr(:), yarr(:), x
    integer,  intent(in)  :: m
    real(dp), intent(out) :: y, dy
    !-------------------------------
    integer :: k,n
    
    n = size(xarr)
    k = ceiling((x - xarr(1))/(xarr(2) - xarr(1)))
    k = min(max(k-(m-1)/2,1),n+1-m)
    call polint(xarr(k:k+m-1),yarr(k:k+m-1), x, y, dy)
  end subroutine interpolate_uniform

  !----------------------------------------------------------------------
  ! In an ordered table of xarr, find the interpolation which uses m
  ! points
  subroutine interpolate(xarr,yarr, x, y, dy, m)
    real(dp), intent(in)  :: xarr(:), yarr(:), x
    integer,  intent(in)  :: m
    real(dp), intent(out) :: y, dy
    !-------------------------------
    integer :: k,n
    
    n = size(xarr)
    k = locate(xarr,x)
    k = min(max(k-(m-1)/2,1),n+1-m)
    call polint(xarr(k:k+m-1),yarr(k:k+m-1), x, y, dy)

  end subroutine interpolate
  
  !----------------------------------------------------------------------
  function locate(xx,x)
    implicit none
    real(dp), dimension(:), intent(in) :: xx
    real(dp), intent(in) :: x
    integer :: locate
    integer :: n,jl,jm,ju
    logical :: ascnd
    n=size(xx)
    ascnd = (xx(n) >= xx(1))
    jl=0
    ju=n+1
    do
       if (ju-jl <= 1) exit
       jm=(ju+jl)/2
       if (ascnd .eqv. (x >= xx(jm))) then
          jl=jm
       else
          ju=jm
       end if
    end do
    if (x == xx(1)) then
       locate=1
    else if (x == xx(n)) then
       locate=n-1
    else
       locate=jl
    end if
  end function locate


  !----------------------------------------------------------------------
  subroutine hunt(xx,x,jlo)
    implicit none
    integer, intent(inout) :: jlo
    real(dp), intent(in) :: x
    real(dp), dimension(:), intent(in) :: xx
    integer :: n,inc,jhi,jm
    logical :: ascnd
    n=size(xx)
    ascnd = (xx(n) >= xx(1))
    if (jlo <= 0 .or. jlo > n) then
       jlo=0
       jhi=n+1
    else
       inc=1
       if (x >= xx(jlo) .eqv. ascnd) then
          do
             jhi=jlo+inc
             if (jhi > n) then
                jhi=n+1
                exit
             else
                if (x < xx(jhi) .eqv. ascnd) exit
                jlo=jhi
                inc=inc+inc
             end if
          end do
       else
          jhi=jlo
          do
             jlo=jhi-inc
             if (jlo < 1) then
                jlo=0
                exit
             else
                if (x >= xx(jlo) .eqv. ascnd) exit
                jhi=jlo
                inc=inc+inc
             end if
          end do
       end if
    end if
    do
       if (jhi-jlo <= 1) then
          if (x == xx(n)) jlo=n-1
          if (x == xx(1)) jlo=1
          exit
       else
          jm=(jhi+jlo)/2
          if (x >= xx(jm) .eqv. ascnd) then
             jlo=jm
          else
             jhi=jm
          end if
       end if
    end do
  end subroutine hunt

  !----------------------------------------------------------------------
  subroutine polint(xa,ya,x,y,dy)
    use assertions, only: assert_eq
    implicit none
    real(dp), dimension(:), intent(in) :: xa,ya
    real(dp), intent(in) :: x
    real(dp), intent(out) :: y,dy
    integer :: m,n,ns
    real(dp), dimension(size(xa)) :: c,d,den,ho
    n=assert_eq(size(xa),size(ya),'polint')
    c=ya
    d=ya
    ho=xa-x
    ns=iminloc(abs(x-xa))
    y=ya(ns)
    ns=ns-1
    do m=1,n-1
       den(1:n-m)=ho(1:n-m)-ho(1+m:n)
       if (any(den(1:n-m) == 0.0)) stop 'polint: calculation failure'
       den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
       d(1:n-m)=ho(1+m:n)*den(1:n-m)
       c(1:n-m)=ho(1:n-m)*den(1:n-m)
       if (2*ns < n-m) then
          dy=c(ns+1)
       else
          dy=d(ns)
          ns=ns-1
       end if
       y=y+dy
    end do
  end subroutine polint


  function iminloc(arr)
    real(dp), dimension(:), intent(in) :: arr
    integer,  dimension(1) :: imin
    integer :: iminloc
    imin=minloc(arr(:))
    iminloc=imin(1)
  end function iminloc



  !--------------------------------------------------------------------
  ! returns the weights for the uniform interpolation,
  ! where the first entry of weights corresponds to x=1
  ! and the spacing is 1.
  !
  ! formula should be weight(i) = [Prod_{j/=i} (x-j)] (-1)^(n-1) / (i! (n-i)!)
  !
  ! Algorithm uses instead [Prod_{j} (x-j)] / (x-i) unless one
  ! of the x-i==0 in which case the result is just 0,...,0,,1,0...,0
  !
  ! For simplicity of caching, n == ubound(weights), is limited .le. nmax
  subroutine uniform_interpolation_weights(x,weights)
    use warnings_and_errors; use consts_dp
    real(dp), intent(in)  :: x
    real(dp), intent(out) :: weights(0:)
    !-----------------------------------------
    integer,  parameter :: nmax = 9
    !                                           order=n
    real(dp), save      :: normalisation(0:nmax,0:nmax) = zero
    real(dp)            :: dists(0:ubound(weights,dim=1)), prod
    integer             :: n, i
    
    n = ubound(weights,dim=1)
    if (n > nmax) call wae_error('uniform_interpolation_weights',&
         &'ubound of weights is too large:',intval=n)

    !-- intialise once for each n
    if (normalisation(0,n) == zero) then
       !-- calculate factorial
       normalisation(0,n) = one
       do i = 1, n
          normalisation(0,n) = normalisation(0,n) * (-i)
       end do
       !-- calculate inverse weight "normalisations"
       do i = 1, n
          normalisation(i,n) = (normalisation(i-1,n) * i)/(i-1-n)
       end do
       normalisation(:n,n) = one / normalisation(:n,n)
    end if
    
    do i = 0, n
       dists(i) = x - i
       if (dists(i) == zero) then
          weights(:) = zero
          weights(i) = one
          return
       end if
    end do
    prod = product(dists)
    do i = 0, n
       weights(i) = prod * normalisation(i,n) / dists(i)
    end do
    
  end subroutine uniform_interpolation_weights
  


  subroutine spline (x, y, b, c, d, n)
    !======================================================================
    !  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
    !  for cubic spline interpolation
    !  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
    !  for  x(i) <= x <= x(i+1)
    !  Alex G: January 2010
    !----------------------------------------------------------------------
    !  input..
    !  x = the arrays of data abscissas (in strictly increasing order)
    !  y = the arrays of data ordinates
    !  n = size of the arrays xi() and yi() (n>=2)
    !  output..
    !  b, c, d  = arrays of spline coefficients
    !  comments ...
    !  spline.f90 program is based on fortran version of program spline.f
    !  the accompanying function fspline can be used for interpolation
    !======================================================================
    implicit none
    integer n
    double precision x(n), y(n), b(n), c(n), d(n)
    integer i, j, gap
    double precision h

    gap = n-1
    ! check input
    if ( n < 2 ) return
    if ( n < 3 ) then
       b(1) = (y(2)-y(1))/(x(2)-x(1))   ! linear interpolation
       c(1) = 0.
       d(1) = 0.
       b(2) = b(1)
       c(2) = 0.
       d(2) = 0.
       return
    end if
    !
    ! step 1: preparation
    !
    d(1) = x(2) - x(1)
    c(2) = (y(2) - y(1))/d(1)
    do i = 2, gap
       d(i) = x(i+1) - x(i)
       b(i) = 2.0*(d(i-1) + d(i))
       c(i+1) = (y(i+1) - y(i))/d(i)
       c(i) = c(i+1) - c(i)
    end do
    !
    ! step 2: end conditions 
    !
    b(1) = -d(1)
    b(n) = -d(n-1)
    c(1) = 0.0
    c(n) = 0.0
    if(n /= 3) then
       c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
       c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
       c(1) = c(1)*d(1)**2/(x(4)-x(1))
       c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
    end if
    !
    ! step 3: forward elimination 
    !
    do i = 2, n
       h = d(i-1)/b(i-1)
       b(i) = b(i) - h*d(i-1)
       c(i) = c(i) - h*c(i-1)
    end do
    !
    ! step 4: back substitution
    !
    c(n) = c(n)/b(n)
    do j = 1, gap
       i = n-j
       c(i) = (c(i) - d(i)*c(i+1))/b(i)
    end do
    !
    ! step 5: compute spline coefficients
    !
    b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
    do i = 1, gap
       b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
       d(i) = (c(i+1) - c(i))/d(i)
       c(i) = 3.*c(i)
    end do
    c(n) = 3.0*c(n)
    d(n) = d(n-1)
  end subroutine spline

  function ispline(u, x, y, b, c, d, n)
    !======================================================================
    ! function ispline evaluates the cubic spline interpolation at point z
    ! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
    ! where  x(i) <= u <= x(i+1)
    !----------------------------------------------------------------------
    ! input..
    ! u       = the abscissa at which the spline is to be evaluated
    ! x, y    = the arrays of given data points
    ! b, c, d = arrays of spline coefficients computed by spline
    ! n       = the number of data points
    ! output:
    ! ispline = interpolated value at point u
    !=======================================================================
    implicit none
    double precision ispline
    integer n
    double precision  u, x(n), y(n), b(n), c(n), d(n)
    integer i, j, k
    double precision dx

    ! if u is ouside the x() interval take a boundary value (left or right)
    if(u <= x(1)) then
       ispline = y(1)
       return
    end if
    if(u >= x(n)) then
       ispline = y(n)
       return
    end if

    !*
    !  binary search for for i, such that x(i) <= u <= x(i+1)
    !*
    i = 1
    j = n+1
    do while (j > i+1)
       k = (i+j)/2
       if(u < x(k)) then
          j=k
       else
          i=k
       end if
    end do
    !*
    !  evaluate spline interpolation
    !*
    dx = u - x(i)
    ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
  end function ispline



end module interpolation

