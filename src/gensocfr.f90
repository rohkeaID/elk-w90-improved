
! Copyright (C) 2009 J. K. Dewhurst, S. Sharma and E. K. U. Gross
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine gensocfr
use modmain
implicit none
! local variables
integer is,ias
integer nr,ir,irc
real(8) cso,rm
! allocatable arrays
real(8), allocatable :: vr(:),dvr(:)
if (.not.spinorb) return
! coefficient of spin-orbit coupling
cso=socscf/(4.d0*solsc**2)
allocate(vr(nrmtmax),dvr(nrmtmax))
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
! radial derivative of the spherical part of the Kohn-Sham potential
  vr(1:nr)=vsmt(1,1:nr,ias)*y00
  call fderiv(1,nr,rsp(:,is),vr,dvr)
  irc=0
  do ir=1,nr,lradstp
    irc=irc+1
    rm=1.d0-2.d0*cso*vr(ir)
    socfr(irc,ias)=cso*dvr(ir)/(rsp(ir,is)*rm**2)
  end do
end do
deallocate(vr,dvr)
return
end subroutine

