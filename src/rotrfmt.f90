
! Copyright (C) 2015 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rotrfmt(rot,nr,nri,lrstp,rfmt1,rfmt2)
use modmain
implicit none
! arguments
real(8), intent(in) :: rot(3,3)
integer, intent(in) :: nr,nri
integer, intent(in) :: lrstp
real(8), intent(in) :: rfmt1(lmmaxvr,nr)
real(8), intent(out) :: rfmt2(lmmaxvr,nr)
! local variables
integer nrci,nrco,iro,ld
ld=lmmaxvr*lrstp
! inner part of muffin-tin
nrci=(nri-1)/lrstp+1
call rotrflm(rot,lmaxinr,nrci,ld,rfmt1,rfmt2)
! outer part of muffin-tin
nrco=(nr-nri)/lrstp
if (nrco.eq.0) return
iro=nri+lrstp
call rotrflm(rot,lmaxvr,nrco,ld,rfmt1(:,iro),rfmt2(:,iro))
return
end subroutine

