
! Copyright (C) 2014 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine proj2d(rvfmt,rvfir)
use modmain
implicit none
! arguments
real(8), intent(inout) :: rvfmt(lmmaxvr,nrmtmax,natmtot,3)
real(8), intent(inout) :: rvfir(ngtot,3)
! local variables
integer is,ias,ir,lm
real(8) vl1(3),vl2(3),t1,t2,t3
real(8) vc1(3),vc2(3),vc3(3),vc4(3)
! determine the projection onto the plotting plane
vl1(:)=vclp2d(:,2)-vclp2d(:,1)
vl2(:)=vclp2d(:,3)-vclp2d(:,1)
call r3mv(avec,vl1,vc1)
call r3mv(avec,vl2,vc2)
call r3cross(vc1,vc2,vc3)
t1=sqrt(vc1(1)**2+vc1(2)**2+vc1(3)**2)
t2=sqrt(vc2(1)**2+vc2(2)**2+vc2(3)**2)
t3=sqrt(vc3(1)**2+vc3(2)**2+vc3(3)**2)
if ((t1.lt.epslat).or.(t2.lt.epslat).or.(t3.lt.epslat)) then
  write(*,*)
  write(*,'("Error(proj2d): degenerate 2D plotting directions")')
  write(*,*)
  stop
end if
vc1(:)=vc1(:)/t1
vc2(:)=vc2(:)/t2
vc3(:)=vc3(:)/t3
call r3cross(vc3,vc1,vc2)
t1=sqrt(vc2(1)**2+vc2(2)**2+vc2(3)**2)
vc2(:)=vc2(:)/t1
! muffin-tin part
do ias=1,natmtot
  is=idxis(ias)
  do ir=1,nrmt(is)
    do lm=1,lmmaxvr
      vc4(:)=rvfmt(lm,ir,ias,:)
      rvfmt(lm,ir,ias,1)=dot_product(vc4(:),vc1(:))
      rvfmt(lm,ir,ias,2)=dot_product(vc4(:),vc2(:))
      rvfmt(lm,ir,ias,3)=dot_product(vc4(:),vc3(:))
    end do
  end do
end do
! interstitial part
do ir=1,ngtot
  vc4(:)=rvfir(ir,:)
  rvfir(ir,1)=dot_product(vc4(:),vc1(:))
  rvfir(ir,2)=dot_product(vc4(:),vc2(:))
  rvfir(ir,3)=dot_product(vc4(:),vc3(:))
end do
return
end subroutine

