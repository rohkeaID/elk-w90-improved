
! Copyright (C) 2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gradrf(rfmt,rfir,grfmt,grfir)
use modmain
implicit none
! arguments
real(8), intent(in) :: rfmt(lmmaxvr,nrmtmax,natmtot),rfir(ngtot)
real(8), intent(out) :: grfmt(lmmaxvr,nrmtmax,natmtot,3),grfir(ngtot,3)
! local variables
integer is,ias,ig,ifg,i
! allocatable arrays
real(8), allocatable :: grfmt1(:,:,:)
complex(8), allocatable :: zfft1(:),zfft2(:)
! muffin-tin gradient
allocate(grfmt1(lmmaxvr,nrmtmax,3))
do ias=1,natmtot
  is=idxis(ias)
  call gradrfmt(nrmt(is),nrmtinr(is),rsp(:,is),rfmt(:,:,ias),nrmtmax,grfmt1)
  do i=1,3
    grfmt(:,1:nrmt(is),ias,i)=grfmt1(:,1:nrmt(is),i)
  end do
end do
deallocate(grfmt1)
! interstitial gradient
allocate(zfft1(ngtot),zfft2(ngtot))
zfft1(:)=rfir(:)
call zfftifc(3,ngridg,-1,zfft1)
do i=1,3
  zfft2(:)=0.d0
  do ig=1,ngvec
    ifg=igfft(ig)
    zfft2(ifg)=vgc(i,ig)*cmplx(-aimag(zfft1(ifg)),dble(zfft1(ifg)),8)
  end do
  call zfftifc(3,ngridg,1,zfft2)
  grfir(:,i)=dble(zfft2(:))
end do
deallocate(zfft1,zfft2)
return
end subroutine

