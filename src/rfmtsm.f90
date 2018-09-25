
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rfmtsm(m,nr,nri,rfmt)
use modmain
implicit none
! arguments
integer, intent(in) :: m,nr,nri
real(8), intent(inout) :: rfmt(*)
! local variables
integer nro,iro,ir
integer lm,npi,i
! automatic arrays
real(8) fr(nr)
if (m.le.0) return
nro=nr-nri
iro=nri+1
npi=lmmaxi*nri
do lm=1,lmmaxi
  i=lm
  do ir=1,nri
    fr(ir)=rfmt(i)
    i=i+lmmaxi
  end do
  do ir=iro,nr
    fr(ir)=rfmt(i)
    i=i+lmmaxo
  end do
  call fsmooth(m,nr,fr)
  i=lm
  do ir=1,nri
    rfmt(i)=fr(ir)
    i=i+lmmaxi
  end do
  do ir=iro,nr
    rfmt(i)=fr(ir)
    i=i+lmmaxo
  end do
end do
do lm=lmmaxi+1,lmmaxo
  i=npi+lm
  do ir=iro,nr
    fr(ir)=rfmt(i)
    i=i+lmmaxo
  end do
  call fsmooth(m,nro,fr(iro))
  i=npi+lm
  do ir=iro,nr
    rfmt(i)=fr(ir)
    i=i+lmmaxo
  end do
end do
return
end subroutine

