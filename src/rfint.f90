
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

real(8) function rfint(rfmt,rfir)
use modmain
implicit none
! arguments
real(8), intent(in) :: rfmt(lmmaxvr,nrmtmax,natmtot),rfir(ngtot)
! local variables
integer is,ias,nr
real(8) sum,t1
! automatic arrays
real(8) fr(nrmtmax)
! external functions
real(8) ddot,fintgt
external ddot,fintgt
! interstitial contribution
sum=ddot(ngtot,rfir,1,cfunir,1)
sum=sum*omega/dble(ngtot)
! muffin-tin contribution
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  fr(1:nr)=rfmt(1,1:nr,ias)*r2sp(1:nr,is)
  t1=fintgt(-1,nr,rsp(:,is),fr)
  sum=sum+fourpi*y00*t1
end do
rfint=sum
return
end function

