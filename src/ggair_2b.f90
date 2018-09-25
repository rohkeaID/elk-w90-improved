
! Copyright (C) 2009 T. McQueen and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: ggair_2b
! !INTERFACE:
subroutine ggair_2b(g2rho,gvrho,vx,vc,dxdgr2,dcdgr2)
! !USES:
use modmain
! !DESCRIPTION:
!   Spin-unpolarised version of {\tt ggair\_sp\_2b}.
!
! !REVISION HISTORY:
!   Created November 2009 (JKD and TMcQ)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: g2rho(ngtot),gvrho(ngtot,3)
real(8), intent(inout) :: vx(ngtot),vc(ngtot)
real(8), intent(in) :: dxdgr2(ngtot),dcdgr2(ngtot)
! local variables
integer ig,ifg,i
! allocatable arrays
real(8), allocatable :: rfir(:)
complex(8), allocatable :: zfft1(:),zfft2(:)
allocate(rfir(ngtot))
allocate(zfft1(ngtot),zfft2(ngtot))
!------------------!
!     exchange     !
!------------------!
! compute grad dxdgr2
zfft1(:)=dxdgr2(:)
call zfftifc(3,ngridg,-1,zfft1)
! (grad dxdgr2).(grad rho)
rfir(:)=0.d0
do i=1,3
  zfft2(:)=0.d0
  do ig=1,ngvec
    ifg=igfft(ig)
    zfft2(ifg)=vgc(i,ig)*cmplx(-aimag(zfft1(ifg)),dble(zfft1(ifg)),8)
  end do
  call zfftifc(3,ngridg,1,zfft2)
  rfir(:)=rfir(:)+dble(zfft2(:))*gvrho(:,i)
end do
vx(:)=vx(:)-2.d0*(rfir(:)+dxdgr2(:)*g2rho(:))
!---------------------!
!     correlation     !
!---------------------!
! compute grad dcdgr2
zfft1(:)=dcdgr2(:)
call zfftifc(3,ngridg,-1,zfft1)
! (grad dcdgr2).(grad rho)
rfir(:)=0.d0
do i=1,3
  zfft2(:)=0.d0
  do ig=1,ngvec
    ifg=igfft(ig)
    zfft2(ifg)=vgc(i,ig)*cmplx(-aimag(zfft1(ifg)),dble(zfft1(ifg)),8)
  end do
  call zfftifc(3,ngridg,1,zfft2)
  rfir(:)=rfir(:)+dble(zfft2(:))*gvrho(:,i)
end do
vc(:)=vc(:)-2.d0*(rfir(:)+dcdgr2(:)*g2rho(:))
deallocate(rfir,zfft1,zfft2)
return
end subroutine
!EOC

