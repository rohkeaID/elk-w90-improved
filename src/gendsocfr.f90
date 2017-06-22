
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine gendsocfr
use modmain
use modphonon
implicit none
integer is,ias
integer nr,ir,irc
real(8) cso
complex(8) z1
! allocatable arrays
real(8), allocatable :: vr1(:),vr2(:)
real(8), allocatable :: dvr1(:),dvr2(:)
if (.not.spinorb) return
! coefficient of spin-orbit coupling
cso=1.d0/(4.d0*solsc**2)
allocate(vr1(nrmtmax),vr2(nrmtmax))
allocate(dvr1(nrmtmax),dvr2(nrmtmax))
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  vr1(1:nr)=dble(dvsmt(1,1:nr,ias))*y00
  vr2(1:nr)=aimag(dvsmt(1,1:nr,ias))*y00
  call fderiv(1,nr,rsp(:,is),vr1,dvr1)
  call fderiv(1,nr,rsp(:,is),vr2,dvr2)
  irc=0
  do ir=1,nr,lradstp
    irc=irc+1
    z1=cmplx(dvr1(ir),dvr2(ir),8)
    dsocfr(irc,ias)=(cso/rsp(ir,is))*z1
  end do
end do
deallocate(vr1,vr2,dvr1,dvr2)
return
end subroutine

