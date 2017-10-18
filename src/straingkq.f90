
! Copyright (C) 2016 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine straingkq
use modmain
use modultra
use modstore
implicit none
integer is,ia,ig
integer nppt,ik,igk
integer jspn,iq
real(8) ta(3,3),tb(3,3),v(3)
if ((istrain.lt.1).or.(istrain.gt.nstrain)) return
! compute the strained lattice vectors
avec(:,:)=avec0(:,:)+deltast*strain(:,:,istrain)
! generate the strained reciprocal lattice vectors and unit cell volume
call reciplat(avec,bvec,omega,omegabz)
! determine the transformation matrix to the strained vectors
call r3mm(avec,ainv,ta)
call r3mm(bvec,binv,tb)
! recalculate all required variables which depend on avec
call r3minv(avec,ainv)
call r3minv(bvec,binv)
call r3mv(bvec,vqlss,vqcss)
do is=1,nspecies
  do ia=1,natoms(is)
    call r3mv(avec,atposl(:,ia,is),atposc(:,ia,is))
  end do
end do
if (ultracell) then
  call r3mm(ta,avecu0,avecu)
  call reciplat(avecu,bvecu,omegau,omegabzu)
end if
call r3mv(bvec,vecql,vecqc)
call r3mv(ainv,efieldc,efieldl)
! apply the transformation matrix to the G-vectors
do ig=1,ngtot
  v(:)=vgc(:,ig)
  call r3mv(tb,v,vgc(:,ig))
  gc(ig)=sqrt(vgc(1,ig)**2+vgc(2,ig)**2+vgc(3,ig)**2)
end do
! recalculate variables which depend on the G-vectors
call genjlgprmt(lnpsd,ngvec,gc,ngvec,jlgrmt)
call genylmg
call gensfacgp(ngvec,vgc,ngvec,sfacg)
do is=1,nspecies
  call genffacgp(is,gc,ffacg(:,is))
end do
call gencfun
call energynn
! apply the transformation to the k-vectors
do ik=1,nkptnr
  v(:)=vkc(:,ik)
  call r3mv(tb,v,vkc(:,ik))
end do
call genkpakq
! apply the transformation to G+k-vectors and recalculate dependent variables
if (xctype(1).lt.0) then
  nppt=nkptnr
else
  nppt=nkpt
end if
do ik=1,nppt
  do jspn=1,nspnfv
    do igk=1,ngk(jspn,ik)
      v(:)=vgkc(:,igk,jspn,ik)
      call r3mv(tb,v,vgkc(:,igk,jspn,ik))
      call sphcrd(vgkc(:,igk,jspn,ik),gkc(igk,jspn,ik),tpgkc(:,igk,jspn,ik))
    end do
    call gensfacgp(ngk(jspn,ik),vgkc(:,:,jspn,ik),ngkmax,sfacgk(:,:,jspn,ik))
  end do
end do
! apply the transformation to the q-vectors if required
if (xctype(1).lt.0) then
  do iq=1,nqptnr
    v(:)=vqc(:,iq)
    call r3mv(tb,v,vqc(:,iq))
  end do
  call genwiq2
end if
return
end subroutine

