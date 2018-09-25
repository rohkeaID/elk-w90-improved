
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: symrfir
! !INTERFACE:
subroutine symrfir(ngv,rfir)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ngv  : number of G-vectors to be used for the Fourier space rotation
!          (in,integer)
!   rfir : real intersitial function (inout,real(ngtot))
! !DESCRIPTION:
!   Symmetrises a real scalar interstitial function. The function is first
!   Fourier transformed to $G$-space, and then averaged over each symmetry by
!   rotating the Fourier coefficients and multiplying them by a phase factor
!   corresponding to the symmetry translation.
!
! !REVISION HISTORY:
!   Created July 2007 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ngv
real(8), intent(inout) :: rfir(ngtot)
! local variables
logical tvz
integer isym,lspl,ilspl,sym(3,3)
integer iv(3),jv(3),ig,ifg,jfg
real(8) vtc(3),t1
complex(8) z1
! allocatable arrays
complex(8), allocatable :: zfft1(:),zfft2(:)
allocate(zfft1(ngtot),zfft2(ngtot))
! Fourier transform function to G-space
zfft1(:)=rfir(:)
call zfftifc(3,ngridg,-1,zfft1)
zfft2(:)=0.d0
! loop over crystal symmetries
do isym=1,nsymcrys
! translation in Cartesian coordinates
  call r3mv(avec,vtlsymc(:,isym),vtc)
! index to lattice symmetry of spatial rotation
  lspl=lsplsymc(isym)
! inverse rotation required for rotation of G-vectors
  ilspl=isymlat(lspl)
  sym(:,:)=symlat(:,:,ilspl)
! zero translation vector flag
  tvz=tv0symc(isym)
  do ig=1,ngv
    ifg=igfft(ig)
! multiply the transpose of the inverse symmetry matrix with the G-vector
    if (lspl.eq.1) then
      jfg=ifg
    else
      iv(:)=ivg(:,ig)
      jv(1)=sym(1,1)*iv(1)+sym(2,1)*iv(2)+sym(3,1)*iv(3)
      jv(2)=sym(1,2)*iv(1)+sym(2,2)*iv(2)+sym(3,2)*iv(3)
      jv(3)=sym(1,3)*iv(1)+sym(2,3)*iv(2)+sym(3,3)*iv(3)
      jfg=igfft(ivgig(jv(1),jv(2),jv(3)))
    end if
    if (tvz) then
! zero translation vector
      zfft2(jfg)=zfft2(jfg)+zfft1(ifg)
    else
! complex phase factor for translation
      t1=-(vgc(1,ig)*vtc(1)+vgc(2,ig)*vtc(2)+vgc(3,ig)*vtc(3))
      z1=cmplx(cos(t1),sin(t1),8)
      zfft2(jfg)=zfft2(jfg)+z1*zfft1(ifg)
    end if
  end do
end do
! Fourier transform to real-space and normalise
call zfftifc(3,ngridg,1,zfft2)
t1=1.d0/dble(nsymcrys)
rfir(:)=t1*dble(zfft2(:))
deallocate(zfft1,zfft2)
return
end subroutine
!EOC

