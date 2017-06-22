
! Copyright (C) 2015 D. Ernsting, S. Dugdale and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine emdplot1d(emds)
use modmain
use modpw
implicit none
! arguments
real(4), intent(in) :: emds(nhkmax,nkpt)
! local variables
integer nh(3),ip,n,i,j
real(8) vl1(3),vl2(3),vl3(3)
real(8) vc1(3),vc2(3),vc3(3),t1
! allocatable arrays
real(8), allocatable :: x(:),f1(:),f2(:),g(:)
! external functions
real(8) rfhkintp
external rfhkintp
! generate the 1D plotting points: use only the first segment
call plotpt1d(bvec,2,npp1d,vvlp1d,vplp1d,dvp1d,dpp1d)
! compute two vectors orthogonal to each other and the plotting vector; these
! are the directions to be used for integration
vl1(:)=vvlp1d(:,2)-vvlp1d(:,1)
call r3mv(bvec,vl1,vc1)
t1=sqrt(vc1(1)**2+vc1(2)**2+vc1(3)**2)
if (t1.lt.epslat) then
  write(*,*)
  write(*,'("Error(emdplot1d): zero length plotting vector")')
  write(*,*)
  stop
end if
vc1(:)=vc1(:)/t1
i=1
do j=2,3
  if (abs(vc1(j)).lt.abs(vc1(i))) i=j
end do
vc2(:)=0.d0
vc2(i)=1.d0
t1=dot_product(vc1,vc2)
vc2(:)=vc2(:)-t1*vc1(:)
t1=sqrt(vc2(1)**2+vc2(2)**2+vc2(3)**2)
vc2(:)=vc2(:)/t1
call r3cross(vc1,vc2,vc3)
! integration directions in lattice coordinates
call r3mv(binv,vc2,vl2)
call r3mv(binv,vc3,vl3)
! determine the number of integration points
nh(:)=int(hkmax*sqrt(avec(1,:)**2+avec(2,:)**2+avec(3,:)**2)/pi)+1
n=2*maxval(nh(:)*ngridk(:))
allocate(x(n))
do i=1,n
  t1=2.d0*dble(i-1)/dble(n-1)-1.d0
  x(i)=t1*hkmax
end do
open(50,file='EMD1D.OUT',action='WRITE',form='FORMATTED')
write(*,*)
! loop over plotting points along 1D line
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(f1,f2,g,i,j,vl1)
!$OMP DO ORDERED
do ip=1,npp1d
  allocate(f1(n),f2(n),g(n))
  do i=1,n
    do j=1,n
      vl1(:)=vplp1d(:,ip)+x(i)*vl2(:)+x(j)*vl3(:)
      f1(j)=rfhkintp(vl1,emds)
    end do
    call fderiv(-2,n,x,f1,g)
    f2(i)=g(n)
  end do
  call fderiv(-2,n,x,f2,g)
!$OMP ORDERED
  write(*,'("Info(emdplot1d): done ",I6," of ",I6," points")') ip,npp1d
  write(50,'(2G18.10)') dpp1d(ip),g(n)
  call flushifc(50)
!$OMP END ORDERED
  deallocate(f1,f2,g)
end do
!$OMP END DO
!$OMP END PARALLEL
close(50)
deallocate(x)
return
end subroutine

