
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genrmesh
! !INTERFACE:
subroutine genrmesh
! !USES:
use modmain
use modvars
! !DESCRIPTION:
!   Generates the coarse and fine radial meshes for each atomic species in the
!   crystal. Also determines which points are in the inner part of the
!   muffin-tin using the value of {\tt fracinr}.
!
! !REVISION HISTORY:
!   Created September 2002 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ir,irc
real(8) t1,t2
! estimate the number of radial mesh points to infinity
nrspmax=1
do is=1,nspecies
! logarithmic mesh
  t1=log(rmaxsp(is)/rminsp(is))/log(rmt(is)/rminsp(is))
  t2=dble(nrmt(is)-1)*t1
  nrsp(is)=nint(t2)+1
  nrspmax=max(nrspmax,nrsp(is))
end do
! generate the radial meshes
if (allocated(rsp)) deallocate(rsp)
allocate(rsp(nrspmax,nspecies))
if (allocated(r2sp)) deallocate(r2sp)
allocate(r2sp(nrspmax,nspecies))
do is=1,nspecies
  t1=1.d0/dble(nrmt(is)-1)
! logarithmic mesh
  t2=log(rmt(is)/rminsp(is))
  do ir=1,nrsp(is)
    rsp(ir,is)=rminsp(is)*exp(dble(ir-1)*t1*t2)
    r2sp(ir,is)=rsp(ir,is)**2
  end do
end do
! set up the coarse radial meshes and find the inner part of the muffin-tin
! where rho is calculated with lmaxi
if (allocated(rcmt)) deallocate(rcmt)
allocate(rcmt(nrcmtmax,nspecies))
if (allocated(r2cmt)) deallocate(r2cmt)
allocate(r2cmt(nrcmtmax,nspecies))
do is=1,nspecies
  t1=fracinr*rmt(is)
  nrmti(is)=1
  nrcmti(is)=1
  irc=0
  do ir=1,nrmt(is),lradstp
    irc=irc+1
    rcmt(irc,is)=rsp(ir,is)
    r2cmt(irc,is)=r2sp(ir,is)
    if (rsp(ir,is).lt.t1) then
      nrmti(is)=ir
      nrcmti(is)=irc
    end if
  end do
end do
! write to VARIABLES.OUT
call writevars('nrsp',nv=nspecies,iva=nrsp)
call writevars('nrmt',nv=nspecies,iva=nrmt)
call writevars('nrmti',nv=nspecies,iva=nrmti)
call writevars('lradstp',iv=lradstp)
call writevars('nrcmt',nv=nspecies,iva=nrcmt)
call writevars('nrcmti',nv=nspecies,iva=nrcmti)
do is=1,nspecies
  call writevars('rsp',nv=nrmt(is),rva=rsp(:,is))
end do
return
end subroutine
!EOC

