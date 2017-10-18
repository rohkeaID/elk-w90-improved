
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rvfcross
! !INTERFACE:
subroutine rvfcross(rvfmt1,rvfir1,rvfmt2,rvfir2,rvfmt3,rvfir3)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   rvfmt1 : first input muffin-tin field (in,real(npmtmax,natmtot,3))
!   rvfir1 : first input interstitial field (in,real(ngtot,3))
!   rvfmt2 : second input muffin-tin field (in,real(npmtmax,natmtot,3))
!   rvfir2 : second input interstitial field (in,real(ngtot,3))
!   rvfmt3 : output muffin-tin field (out,real(npmtmax,natmtot,3))
!   rvfir3 : output interstitial field (out,real(ngtot,3))
! !DESCRIPTION:
!   Given two real vector fields, ${\bf f}_1$ and ${\bf f}_2$, defined over the
!   entire unit cell, this routine computes the local cross product
!   $$ {\bf f}_3({\bf r})\equiv{\bf f}_1({\bf r})\times{\bf f}_2({\bf r}). $$
!
! !REVISION HISTORY:
!   Created February 2007 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: rvfmt1(npmtmax,natmtot,3),rvfir1(ngtot,3)
real(8), intent(in) :: rvfmt2(npmtmax,natmtot,3),rvfir2(ngtot,3)
real(8), intent(out) :: rvfmt3(npmtmax,natmtot,3),rvfir3(ngtot,3)
! local variables
integer is,ias,nr,nri,ir,i
real(8) v1(3),v2(3),v3(3)
! allocatable arrays
real(8), allocatable :: rvfmt4(:,:),rvfmt5(:,:)
!---------------------------!
!     muffin-tin region     !
!---------------------------!
allocate(rvfmt4(npmtmax,3),rvfmt5(npmtmax,3))
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  do i=1,3
    call rbsht(nr,nri,rvfmt1(:,ias,i),rvfmt4(:,i))
    call rbsht(nr,nri,rvfmt2(:,ias,i),rvfmt5(:,i))
  end do
  do i=1,npmt(is)
    v1(:)=rvfmt4(i,:)
    v2(:)=rvfmt5(i,:)
    call r3cross(v1,v2,v3)
    rvfmt4(i,:)=v3(:)
  end do
  do i=1,3
    call rfsht(nr,nri,rvfmt4(:,i),rvfmt3(:,ias,i))
  end do
end do
deallocate(rvfmt4,rvfmt5)
!-----------------------------!
!     interstitial region     !
!-----------------------------!
do ir=1,ngtot
  v1(:)=rvfir1(ir,:)
  v2(:)=rvfir2(ir,:)
  call r3cross(v1,v2,v3)
  rvfir3(ir,:)=v3(:)
end do
return
end subroutine
!EOC

