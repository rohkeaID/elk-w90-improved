
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine ssfext(iq,fext)
use modmain
implicit none
! arguments
integer, intent(in) :: iq
character(*), intent(out) :: fext
! local variables
integer i,j,m(3),n(3)
! external functions
integer gcd
external gcd
do i=1,3
  if (ivq(i,iq).ne.0) then
    j=gcd(ivq(i,iq),ngridq(i))
    m(i)=ivq(i,iq)/j
    n(i)=ngridq(i)/j
  else
    m(i)=0
    n(i)=0
  end if
end do
write(fext,'("_Q",2I2.2,"_",2I2.2,"_",2I2.2,".OUT")') m(1),n(1),m(2),n(2), &
 m(3),n(3)
return
end subroutine

