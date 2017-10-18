
! Copyright (C) 2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gradrf(rfmt,rfir,grfmt,grfir)
use modmain
implicit none
! arguments
real(8), intent(in) :: rfmt(npmtmax,natmtot),rfir(ngtot)
real(8), intent(out) :: grfmt(npmtmax,natmtot,3),grfir(ngtot,3)
! local variables
integer is,ias,np
integer ig,ifg,i
! allocatable arrays
real(8), allocatable :: grfmt1(:,:)
complex(8), allocatable :: zfft1(:),zfft2(:)
! muffin-tin gradient
allocate(grfmt1(npmtmax,3))
do ias=1,natmtot
  is=idxis(ias)
  np=npmt(is)
  call gradrfmt(nrmt(is),nrmti(is),rsp(:,is),rfmt(:,ias),npmtmax,grfmt1)
  do i=1,3
    grfmt(1:np,ias,i)=grfmt1(1:np,i)
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

