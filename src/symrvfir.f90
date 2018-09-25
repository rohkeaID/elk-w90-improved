
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: symrvfir
subroutine symrvfir(tspin,tnc,rvfir)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   tspin : .true. if spin rotations should be used (in,logical)
!   tnc   : .true. if the vector field is non-collinear, otherwise it is
!           collinear along the z-axis (in,logical)
!   rvfir : real interstitial vector function (inout,real(ngtot,*))
! !DESCRIPTION:
!   Symmetrises a real interstitial vector function. See routines {\tt symrvf}
!   and {\tt symrfir} for details.
!
! !REVISION HISTORY:
!   Created July 2007 (JKD)
!EOP
!BOC
implicit none
! arguments
logical, intent(in) :: tspin,tnc
real(8), intent(inout) :: rvfir(ngtot,*)
! local variables
integer isym,lspl,ilspl,lspn
integer nd,sym(3,3),iv(3),jv(3)
integer ig,ifg,jfg,i
real(8) sc(3,3),vtc(3),t1
complex(8) zv1(3),zv2(3),z1
! allocatable arrays
complex(8), allocatable :: zfft1(:,:),zfft2(:,:)
! dimension of the vector field
if (tnc) then
  nd=3
else
  nd=1
end if
allocate(zfft1(ngtot,nd),zfft2(ngtot,nd))
! Fourier transform vector function to G-space
do i=1,nd
  zfft1(:,i)=rvfir(:,i)
  call zfftifc(3,ngridg,-1,zfft1(:,i))
end do
zfft2(:,:)=0.d0
do isym=1,nsymcrys
! translation vector in Cartesian coordinates
  call r3mv(avec,vtlsymc(:,isym),vtc)
! index to spatial rotation lattice symmetry
  lspl=lsplsymc(isym)
! inverse rotation required for rotation of G-vectors
  ilspl=isymlat(lspl)
  sym(:,:)=symlat(:,:,ilspl)
  if (tspin) then
! global spin proper rotation in Cartesian coordinates
    lspn=lspnsymc(isym)
    sc(:,:)=symlatd(lspn)*symlatc(:,:,lspn)
  else
! set spin rotation equal to spatial rotation
    lspn=lspl
    sc(:,:)=symlatc(:,:,lspl)
  end if
  do ig=1,ngvec
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
! complex phase factor for translation
    t1=-(vgc(1,ig)*vtc(1)+vgc(2,ig)*vtc(2)+vgc(3,ig)*vtc(3))
    z1=cmplx(cos(t1),sin(t1),8)
! translation, spatial rotation and global spin rotation
    if (lspn.eq.1) then
! global spin symmetry is the identity
      zfft2(jfg,:)=zfft2(jfg,:)+z1*zfft1(ifg,:)
    else
      if (tnc) then
! non-collinear case
        zv1(:)=zfft1(ifg,:)
        zv2(1)=sc(1,1)*zv1(1)+sc(1,2)*zv1(2)+sc(1,3)*zv1(3)
        zv2(2)=sc(2,1)*zv1(1)+sc(2,2)*zv1(2)+sc(2,3)*zv1(3)
        zv2(3)=sc(3,1)*zv1(1)+sc(3,2)*zv1(2)+sc(3,3)*zv1(3)
        zfft2(jfg,:)=zfft2(jfg,:)+z1*zv2(:)
      else
! collinear case
        zfft2(jfg,1)=zfft2(jfg,1)+sc(3,3)*z1*zfft1(ifg,1)
      end if
    end if
  end do
end do
! Fourier transform to real-space and normalise
t1=1.d0/dble(nsymcrys)
do i=1,nd
  call zfftifc(3,ngridg,1,zfft2(:,i))
  rvfir(:,i)=t1*dble(zfft2(:,i))
end do
deallocate(zfft1,zfft2)
return
end subroutine
!EOC

