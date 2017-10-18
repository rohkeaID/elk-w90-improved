
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: rfinterp
! !INTERFACE:
subroutine rfinterp(ni,xi,fi,no,xo,fo)
! !INPUT/OUTPUT PARAMETERS:
!   ni  : number of input points (in,integer)
!   xi  : input abscissa array (in,real(ni))
!   fi  : input data array (in,real(ni)
!   no  : number of output points (in,integer)
!   xo  : output abscissa array (in,real(ni))
!   fo  : output interpolated function (out,real(no))
! !DESCRIPTION:
!   Given a function defined on a set of input points, this routine uses a
!   clamped cubic spline to interpolate the function on a different set of
!   points. See routine {\tt spline}.
!
! !REVISION HISTORY:
!   Created January 2005 (JKD)
!   Arguments changed, April 2016 (JKD)
!EOP
!BOC
implicit none
integer, intent(in) :: ni
real(8), intent(in) :: xi(ni),fi(ni)
integer, intent(in) :: no
real(8), intent(in) :: xo(no)
real(8), intent(out) :: fo(no)
! local variables
integer i,j,k,l
real(8) x,dx
! automatic arrays
real(8) cf(3,ni)
if (ni.le.0) then
  write(*,*)
  write(*,'("Error(rfinterp): invalid number of input points : ",I8)') ni
  write(*,*)
  stop
end if
if (no.le.0) then
  write(*,*)
  write(*,'("Error(rfinterp): invalid number of output points : ",I8)') no
  write(*,*)
  stop
end if
if (ni.eq.1) then
  fo(:)=fi(1)
  return
end if
! compute the spline coefficients
call spline(ni,xi,fi,cf)
! evaluate spline at output points
i=1
do l=1,no
  x=xo(l)
  if (i.ge.ni) i=1
  if (x.lt.xi(i)) goto 10
  if (x.le.xi(i+1)) goto 30
! binary search
10 continue
  i=1
  j=ni+1
20 continue
  k=(i+j)/2
  if (x.lt.xi(k)) then
    j=k
  else
    i=k
  end if
  if (j.gt.i+1) goto 20
30 continue
  dx=x-xi(i)
  fo(l)=fi(i)+dx*(cf(1,i)+dx*(cf(2,i)+dx*cf(3,i)))
end do
return
end subroutine
!EOC

