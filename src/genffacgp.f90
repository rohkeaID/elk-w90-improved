
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genffacgp(is,gpc,ffacgp)
use modmain
implicit none
! arguments
integer, intent(in) :: is
real(8), intent(in) :: gpc(ngtot)
real(8), intent(out) :: ffacgp(ngtot)
! local variables
integer ig
real(8) t1,t2
t1=fourpi/omega
do ig=1,ngtot
  if (gpc(ig).gt.epslat) then
    t2=gpc(ig)*rmt(is)
    ffacgp(ig)=t1*(sin(t2)-t2*cos(t2))/(gpc(ig)**3)
  else
    ffacgp(ig)=(t1/3.d0)*rmt(is)**3
  end if
end do
return
end subroutine

