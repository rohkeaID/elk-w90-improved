
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: ggair_sp_1
! !INTERFACE:
subroutine ggair_sp_1(rhoup,rhodn,grho,gup,gdn,g2up,g2dn,g3rho,g3up,g3dn)
! !INPUT/OUTPUT PARAMETERS:
!   rhoup : spin-up density (in,real(ngtot))
!   rhodn : spin-down density (in,real(ngtot))
!   grho  : |grad rho| (out,real(ngtot))
!   gup   : |grad rhoup| (out,real(ngtot))
!   gdn   : |grad rhodn| (out,real(ngtot))
!   g2up  : grad^2 rhoup (out,real(ngtot))
!   g2dn  : grad^2 rhodn (out,real(ngtot))
!   g3rho : (grad rho).(grad |grad rho|) (out,real(ngtot))
!   g3up  : (grad rhoup).(grad |grad rhoup|) (out,real(ngtot))
!   g3dn  : (grad rhodn).(grad |grad rhodn|) (out,real(ngtot))
! !DESCRIPTION:
!   Computes $|\nabla\rho|$, $|\nabla\rho^{\uparrow}|$,
!   $|\nabla\rho^{\downarrow}|$, $\nabla^2\rho^{\uparrow}$,
!   $\nabla^2\rho^{\downarrow}$, $\nabla\rho\cdot(\nabla|\nabla\rho|)$,
!   $\nabla\rho^{\uparrow}\cdot(\nabla|\nabla\rho^{\uparrow}|)$ and
!   $\nabla\rho^{\downarrow}\cdot(\nabla|\nabla\rho^{\downarrow}|)$ for the
!   interstitial charge density, as required by the generalised gradient
!   approximation functionals of type 1 for spin-polarised densities. See
!   routines {\tt potxc} and {\tt modxcifc}.
!
! !REVISION HISTORY:
!   Created October 2004 (JKD)
!   Simplified and improved, October 2009 (JKD)
!EOP
!BOC
use modmain
implicit none
! arguments
real(8), intent(in) :: rhoup(ngtot)
real(8), intent(in) :: rhodn(ngtot)
real(8), intent(out) :: grho(ngtot)
real(8), intent(out) :: gup(ngtot)
real(8), intent(out) :: gdn(ngtot)
real(8), intent(out) :: g2up(ngtot)
real(8), intent(out) :: g2dn(ngtot)
real(8), intent(out) :: g3rho(ngtot)
real(8), intent(out) :: g3up(ngtot)
real(8), intent(out) :: g3dn(ngtot)
! local variables
integer ig,ifg,i
! allocatable arrays
real(8), allocatable :: gvup(:,:),gvdn(:,:)
complex(8), allocatable :: zfft1(:),zfft2(:)
allocate(gvup(ngtot,3),gvdn(ngtot,3))
allocate(zfft1(ngtot),zfft2(ngtot))
!----------------!
!     rho up     !
!----------------!
zfft1(:)=rhoup(:)
call zfftifc(3,ngridg,-1,zfft1)
! |grad rhoup|
do i=1,3
  zfft2(:)=0.d0
  do ig=1,ngvec
    ifg=igfft(ig)
    zfft2(ifg)=vgc(i,ig)*cmplx(-aimag(zfft1(ifg)),dble(zfft1(ifg)),8)
  end do
  call zfftifc(3,ngridg,1,zfft2)
  gvup(:,i)=dble(zfft2(:))
end do
gup(:)=sqrt(gvup(:,1)**2+gvup(:,2)**2+gvup(:,3)**2)
! grad^2 rhoup
zfft2(:)=0.d0
do ig=1,ngvec
  ifg=igfft(ig)
  zfft2(ifg)=-(gc(ig)**2)*zfft1(ifg)
end do
call zfftifc(3,ngridg,1,zfft2)
g2up(:)=dble(zfft2(:))
! (grad rhoup).(grad |grad rhoup|)
zfft1(:)=gup(:)
call zfftifc(3,ngridg,-1,zfft1)
g3up(:)=0.d0
do i=1,3
  zfft2(:)=0.d0
  do ig=1,ngvec
    ifg=igfft(ig)
    zfft2(ifg)=vgc(i,ig)*cmplx(-aimag(zfft1(ifg)),dble(zfft1(ifg)),8)
  end do
  call zfftifc(3,ngridg,1,zfft2)
  g3up(:)=g3up(:)+gvup(:,i)*dble(zfft2(:))
end do
!------------------!
!     rho down     !
!------------------!
zfft1(:)=rhodn(:)
call zfftifc(3,ngridg,-1,zfft1)
! |grad rhodn|
do i=1,3
  zfft2(:)=0.d0
  do ig=1,ngvec
    ifg=igfft(ig)
    zfft2(ifg)=vgc(i,ig)*cmplx(-aimag(zfft1(ifg)),dble(zfft1(ifg)),8)
  end do
  call zfftifc(3,ngridg,1,zfft2)
  gvdn(:,i)=dble(zfft2(:))
end do
gdn(:)=sqrt(gvdn(:,1)**2+gvdn(:,2)**2+gvdn(:,3)**2)
! grad^2 rhodn
zfft2(:)=0.d0
do ig=1,ngvec
  ifg=igfft(ig)
  zfft2(ifg)=-(gc(ig)**2)*zfft1(ifg)
end do
call zfftifc(3,ngridg,1,zfft2)
g2dn(:)=dble(zfft2(:))
! (grad rhodn).(grad |grad rhodn|)
zfft1(:)=gdn(:)
call zfftifc(3,ngridg,-1,zfft1)
g3dn(:)=0.d0
do i=1,3
  zfft2(:)=0.d0
  do ig=1,ngvec
    ifg=igfft(ig)
    zfft2(ifg)=vgc(i,ig)*cmplx(-aimag(zfft1(ifg)),dble(zfft1(ifg)),8)
  end do
  call zfftifc(3,ngridg,1,zfft2)
  g3dn(:)=g3dn(:)+gvdn(:,i)*dble(zfft2(:))
end do
!-------------!
!     rho     !
!-------------!
! |grad rho|
grho(:)=sqrt((gvup(:,1)+gvdn(:,1))**2 &
            +(gvup(:,2)+gvdn(:,2))**2 &
            +(gvup(:,3)+gvdn(:,3))**2)
! (grad rho).(grad |grad rho|)
zfft1(:)=grho(:)
call zfftifc(3,ngridg,-1,zfft1)
g3rho(:)=0.d0
do i=1,3
  zfft2(:)=0.d0
  do ig=1,ngvec
    ifg=igfft(ig)
    zfft2(ifg)=vgc(i,ig)*cmplx(-aimag(zfft1(ifg)),dble(zfft1(ifg)),8)
  end do
  call zfftifc(3,ngridg,1,zfft2)
  g3rho(:)=g3rho(:)+(gvup(:,i)+gvdn(:,i))*dble(zfft2(:))
end do
deallocate(gvup,gvdn,zfft1,zfft2)
return
end subroutine
!EOC

