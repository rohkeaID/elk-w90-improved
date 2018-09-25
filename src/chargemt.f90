
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine chargemt
use modmain
implicit none
! local variables
integer is,ias
integer nr,nri,ir,i
real(8) t1
! automatic arrays
real(8) fr(nrmtmax)
! external functions
real(8) fintgt
external fintgt
chgmttot=0.d0
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
! extract the l=0 component from the muffin-tin density
  i=1
  do ir=1,nri
    fr(ir)=rhomt(i,ias)*r2sp(ir,is)
    i=i+lmmaxi
  end do
  do ir=nri+1,nr
    fr(ir)=rhomt(i,ias)*r2sp(ir,is)
    i=i+lmmaxo
  end do
  t1=fintgt(-1,nr,rsp(:,is),fr)
  chgmt(ias)=fourpi*y00*t1
  chgmttot=chgmttot+chgmt(ias)
end do
return
end subroutine
