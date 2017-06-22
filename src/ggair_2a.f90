
! Copyright (C) 2009 T. McQueen and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: ggair_2a
! !INTERFACE:
subroutine ggair_2a(g2rho,gvrho,grho2)
! !USES:
use modmain
! !DESCRIPTION:
!   Spin-unpolarised version of {\tt ggair\_sp\_2a}.
!
! !REVISION HISTORY:
!   Created November 2009 (JKD and TMcQ)
!EOP
!BOC
implicit none
! arguments
real(8), intent(out) :: g2rho(ngtot)
real(8), intent(out) :: gvrho(ngtot,3)
real(8), intent(out) :: grho2(ngtot)
! local variables
integer i,ig,ifg
! allocatable arrays
complex(8), allocatable :: zfft1(:),zfft2(:)
allocate(zfft1(ngtot),zfft2(ngtot))
! Fourier transform density to G-space
zfft1(:)=rhoir(:)
call zfftifc(3,ngridg,-1,zfft1)
! grad^2 rho
zfft2(:)=0.d0
do ig=1,ngvec
  ifg=igfft(ig)
  zfft2(ifg)=-(gc(ig)**2)*zfft1(ifg)
end do
call zfftifc(3,ngridg,1,zfft2)
g2rho(:)=dble(zfft2(:))
! grad rho
do i=1,3
  zfft2(:)=0.d0
  do ig=1,ngvec
    ifg=igfft(ig)
    zfft2(ifg)=vgc(i,ig)*cmplx(-aimag(zfft1(ifg)),dble(zfft1(ifg)),8)
  end do
  call zfftifc(3,ngridg,1,zfft2)
  gvrho(:,i)=dble(zfft2(:))
end do
! (grad rho)^2
grho2(:)=gvrho(:,1)**2+gvrho(:,2)**2+gvrho(:,3)**2
deallocate(zfft1,zfft2)
return
end subroutine
!EOC

