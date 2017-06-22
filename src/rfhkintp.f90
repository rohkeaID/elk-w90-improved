
! Copyright (C) 2015 D. Ernsting, S. Dugdale and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

real(8) function rfhkintp(vhpl,rfhk)
use modmain
use modpw
implicit none
! arguments
real(8), intent(in) :: vhpl(3)
real(4), intent(in) :: rfhk(nhkmax,nkpt)
! local variables
integer ivh0(3),ivk0(3),ihk
integer ivhb(3,0:1,0:1,0:1)
integer ivkb(3,0:1,0:1,0:1)
integer isym,lspl,ik,jk,i,j,k
real(8) vpl(3),fb(0:1,0:1,0:1)
real(8) f00,f01,f10,f11,f0,f1
real(8) v0(3),v1(3),v2(3),t1,t2
! find the H-vector and k-vector corresponding to the input H+p-vector
ivh0(:)=floor(vhpl(:))
vpl(:)=vhpl(:)-dble(ivh0(:))
v1(:)=vpl(:)*dble(ngridk(:))
ivk0(:)=floor(v1(:))
! determine the corners of the box containing the input point
do i=0,1; do j=0,1; do k=0,1
  ivhb(:,i,j,k)=ivh0(:)
  ivkb(:,i,j,k)=ivk0(:)
  ivkb(1,i,j,k)=ivkb(1,i,j,k)+i
  ivkb(2,i,j,k)=ivkb(2,i,j,k)+j
  ivkb(3,i,j,k)=ivkb(3,i,j,k)+k
  ivhb(:,i,j,k)=ivhb(:,i,j,k)+ivkb(:,i,j,k)/ngridk(:)
  ivkb(:,i,j,k)=modulo(ivkb(:,i,j,k),ngridk(:))
end do; end do; end do
! determine the function at each corner of the box
do i=0,1; do j=0,1; do k=0,1
  fb(i,j,k)=0.d0
! non-reduced k-point index
  jk=ivkiknr(ivkb(1,i,j,k),ivkb(2,i,j,k),ivkb(3,i,j,k))
! H+k-vector at corner of box
  v1(:)=dble(ivhb(:,i,j,k))+vkl(:,jk)
! store the origin of the box
  if ((i.eq.0).and.(j.eq.0).and.(k.eq.0)) v0(:)=v1(:)
! vector in Cartesian coordinates
  v2(:)=bvec(:,1)*v1(1)+bvec(:,2)*v1(2)+bvec(:,3)*v1(3)
! check length is within range
  t1=sqrt(v2(1)**2+v2(2)**2+v2(3)**2)
  if (t1.gt.hkmax) cycle
! find the lattice symmetry which maps the non-reduced to reduced k-point
  call findkpt(vkl(:,jk),isym,ik)
! index to spatial rotation in lattice point group
  lspl=lsplsymc(isym)
  v2(:)=symlat(1,:,lspl)*v1(1)+symlat(2,:,lspl)*v1(2)+symlat(3,:,lspl)*v1(3)
! find the H+k-vector for the reduced k-point
  do ihk=1,nhk(1,ik)
    t1=abs(v2(1)-vhkl(1,ihk,1,ik)) &
      +abs(v2(2)-vhkl(2,ihk,1,ik)) &
      +abs(v2(3)-vhkl(3,ihk,1,ik))
    if (t1.lt.epslat) then
      fb(i,j,k)=rfhk(ihk,ik)
      exit
    end if
  end do
end do; end do; end do
! interpolate function
t2=(vhpl(1)-v0(1))*dble(ngridk(1))
t1=1.d0-t2
f00=fb(0,0,0)*t1+fb(1,0,0)*t2
f01=fb(0,0,1)*t1+fb(1,0,1)*t2
f10=fb(0,1,0)*t1+fb(1,1,0)*t2
f11=fb(0,1,1)*t1+fb(1,1,1)*t2
t2=(vhpl(2)-v0(2))*dble(ngridk(2))
t1=1.d0-t2
f0=f00*t1+f10*t2
f1=f01*t1+f11*t2
t2=(vhpl(3)-v0(3))*dble(ngridk(3))
t1=1.d0-t2
rfhkintp=f0*t1+f1*t2
return
end function

