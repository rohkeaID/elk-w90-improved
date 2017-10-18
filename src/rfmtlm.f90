
! Copyright (C) 2016 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rfmtlm(lm,nr,nri,rfmt,fr)
use modmain
implicit none
! arguments
integer, intent(in) :: lm,nr,nri
real(8), intent(in) :: rfmt(npmtmax)
real(8), intent(out) :: fr(nrmtmax)
! local variables
integer iro,ir,npi,i
iro=nri+1
npi=lmmaxi*nri
if (lm.gt.lmmaxi) then
  fr(1:nri)=0.d0
else
  i=lm
  do ir=1,nri
    fr(ir)=rfmt(i)
    i=i+lmmaxi
  end do
end if
if (lm.gt.lmmaxo) then
  fr(iro:nr)=0.d0
else
  i=npi+lm
  do ir=iro,nr
    fr(ir)=rfmt(i)
    i=i+lmmaxo
  end do
end if
return
end subroutine

