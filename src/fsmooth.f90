
! Copyright (C) 2005 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: fsmooth
! !INTERFACE:
subroutine fsmooth(m,n,f)
! !INPUT/OUTPUT PARAMETERS:
!   m  : number of 3-point running averages to perform (in,integer)
!   n  : number of point (in,integer)
!   f  : function array (inout,real(n))
! !DESCRIPTION:
!   Removes numerical noise from a function by performing $m$ successive
!   3-point running averages on the data. The endpoints are kept fixed.
!
! !REVISION HISTORY:
!   Created December 2005 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: m,n
real(8), intent(inout) :: f(n)
! local variables
integer i,j
! automatic arrays
real(8) g(n)
do i=1,m
  do j=2,n-1
    g(j)=0.3333333333333333333d0*(f(j-1)+f(j)+f(j+1))
  end do
  f(2:n-1)=g(2:n-1)
end do
return
end subroutine
!EOC

