
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

real(8) function rfint(rfmt,rfir)
use modmain
implicit none
! arguments
real(8), intent(in) :: rfmt(npmtmax,natmtot),rfir(ngtot)
! local variables
integer is,ias
integer nr,nri,ir,i
real(8) sum,t1
! automatic arrays
real(8) fr(nrmtmax)
! external functions
real(8) fintgt
external fintgt
! interstitial contribution
sum=dot_product(rfir(:),cfunir(:))
sum=sum*omega/dble(ngtot)
! muffin-tin contribution
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  i=1
  do ir=1,nri
    fr(ir)=rfmt(i,ias)*r2sp(ir,is)
    i=i+lmmaxi
  end do
  do ir=nri+1,nr
    fr(ir)=rfmt(i,ias)*r2sp(ir,is)
    i=i+lmmaxo
  end do
  t1=fintgt(-1,nr,rsp(:,is),fr)
  sum=sum+fourpi*y00*t1
end do
rfint=sum
return
end function

