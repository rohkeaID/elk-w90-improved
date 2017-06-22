
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine achiinit(achi)
use modmain
use modscdft
use modrandom
implicit none
! arguments
complex(8), intent(out) :: achi(nbdg,nbdg)
! local variables
integer i,j
real(8) a,b
do i=1,nbdg
  do j=1,nbdg
    a=rndachi*(randomu()-0.5d0)
    b=rndachi*(randomu()-0.5d0)
    achi(i,j)=cmplx(a,b,8)
  end do
end do
return
end subroutine

