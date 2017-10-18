
! Copyright (C) 2016 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine zfmtlm(lm,nr,nri,zfmt,fr1,fr2)
use modmain
implicit none
! arguments
integer, intent(in) :: lm,nr,nri
complex(8), intent(in) :: zfmt(npmtmax)
real(8), intent(out) :: fr1(nr),fr2(nr)
! local variables
integer iro,ir,npi,i
iro=nri+1
npi=lmmaxi*nri
if (lm.gt.lmmaxi) then
  fr1(1:nri)=0.d0
  fr2(1:nri)=0.d0
else
  i=lm
  do ir=1,nri
    fr1(ir)=dble(zfmt(i))
    fr2(ir)=aimag(zfmt(i))
    i=i+lmmaxi
  end do
end if
if (lm.gt.lmmaxo) then
  fr1(iro:nr)=0.d0
  fr2(iro:nr)=0.d0
else
  i=npi+lm
  do ir=iro,nr
    fr1(ir)=dble(zfmt(i))
    fr2(ir)=aimag(zfmt(i))
    i=i+lmmaxo
  end do
end if
return
end subroutine

