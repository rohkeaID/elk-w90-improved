
! Copyright (C) 2014 D. Ernsting, S. Dugdale and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine emdplot2d(emds)
use modmain
use modpw
implicit none
! arguments
real(4), intent(in) :: emds(nhkmax,nkpt)
! local variables
integer nh(3),np,ip,n,i
real(8) vpnl(3),v1(3),t1
! allocatable arrays
real(8), allocatable :: vpl(:,:),vppc(:,:)
real(8), allocatable :: x(:),f(:),g(:)
! external functions
real(8) rfhkintp
external rfhkintp
! allocate local arrays
np=np2d(1)*np2d(2)
allocate(vpl(3,np),vppc(2,np))
! generate the 2D plotting points
call plotpt2d(bvec,binv,vpnl,vpl,vppc)
! determine the number of integration points
nh(:)=int(hkmax*sqrt(avec(1,:)**2+avec(2,:)**2+avec(3,:)**2)/pi)+1
n=2*maxval(nh(:)*ngridk(:))
allocate(x(n))
do i=1,n
  t1=2.d0*dble(i-1)/dble(n-1)-1.d0
  x(i)=t1*hkmax
end do
open(50,file='EMD2D.OUT',action='WRITE',form='FORMATTED')
write(50,'(2I6," : grid size")') np2d(:)
! loop over plotting points in the 2D plane
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(f,g,i,v1)
!$OMP DO ORDERED
do ip=1,np
  allocate(f(n),g(n))
! integrate along normal to plane
  do i=1,n
    v1(:)=vpl(:,ip)+x(i)*vpnl(:)
    f(i)=rfhkintp(v1,emds)
  end do
  call fderiv(-2,n,x,f,g)
!$OMP ORDERED
  write(50,'(3G18.10)') vppc(1,ip),vppc(2,ip),g(n)
!$OMP END ORDERED
  deallocate(f,g)
end do
!$OMP END DO
!$OMP END PARALLEL
close(50)
deallocate(vpl,vppc,x)
return
end subroutine

