
! Copyright (C) 2017 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gengclgq(treg,iq,ngq,gqc,gclgq)
use modmain
implicit none
! arguments
logical, intent(in) :: treg
integer, intent(in) :: iq,ngq
real(8), intent(in) :: gqc(ngq)
real(8), intent(out) :: gclgq(ngq)
! local variables
integer ig
real(8) t1,t2
if (treg) then
! regularise 1/(G+q)^2 for G+q in the first Brillouin zone
  t1=sqrt(vqc(1,iq)**2+vqc(2,iq)**2+vqc(3,iq)**2)
  do ig=1,ngq
    t2=gqc(ig)
    if (abs(t1-t2).lt.epslat) then
      gclgq(ig)=gclq(iq)
    else
      gclgq(ig)=fourpi/t2**2
    end if
  end do
else
! no regularisation
  do ig=1,ngq
    t1=gqc(ig)
    if (t1.gt.epslat) then
      gclgq(ig)=fourpi/t1**2
    else
      gclgq(ig)=0.d0
    end if
  end do
end if
return
end subroutine

