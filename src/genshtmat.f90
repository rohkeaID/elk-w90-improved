
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genshtmat
! !INTERFACE:
subroutine genshtmat
! !USES:
use modmain
! !DESCRIPTION:
!   Generates the forward and backward spherical harmonic transformation (SHT)
!   matrices using the spherical covering set produced by the routine
!   {\tt sphcover}. These matrices are used to transform a function between its
!   $(l,m)$-expansion coefficients and its values at the $(\theta,\phi)$ points
!   on the sphere.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer itp,info
real(8) v1(3),v2(3),r
! automatic arrays
integer ipiv(lmmaxo)
real(8) tp(2,lmmaxo),rlm(lmmaxo),work(lmmaxo)
complex(8) ylm(lmmaxo),zwork(lmmaxo)
!--------------------------------!
!     SHT matrices for lmaxo     !
!--------------------------------!
! allocate real SHT matrices
if (allocated(rbshto)) deallocate(rbshto)
allocate(rbshto(lmmaxo,lmmaxo))
if (allocated(rfshto)) deallocate(rfshto)
allocate(rfshto(lmmaxo,lmmaxo))
! allocate complex SHT matrices
if (allocated(zbshto)) deallocate(zbshto)
allocate(zbshto(lmmaxo,lmmaxo))
if (allocated(zfshto)) deallocate(zfshto)
allocate(zfshto(lmmaxo,lmmaxo))
! generate spherical covering set
call sphcover(lmmaxo,tp)
! rotate the spherical covering set if required
if (trotsht) then
  do itp=1,lmmaxo
    call sctovec(tp(:,itp),v1)
    call r3mv(rotsht,v1,v2)
    call sphcrd(v2,r,tp(:,itp))
  end do
end if
! generate real and complex spherical harmonics and set the backward SHT arrays
do itp=1,lmmaxo
  call genrlm(lmaxo,tp(:,itp),rlm)
  rbshto(itp,1:lmmaxo)=rlm(1:lmmaxo)
  call genylm(lmaxo,tp(:,itp),ylm)
  zbshto(itp,1:lmmaxo)=ylm(1:lmmaxo)
end do
! find the forward SHT arrays
! real
rfshto(:,:)=rbshto(:,:)
call dgetrf(lmmaxo,lmmaxo,rfshto,lmmaxo,ipiv,info)
if (info.ne.0) goto 10
call dgetri(lmmaxo,rfshto,lmmaxo,ipiv,work,lmmaxo,info)
if (info.ne.0) goto 10
! complex
zfshto(:,:)=zbshto(:,:)
call zgetrf(lmmaxo,lmmaxo,zfshto,lmmaxo,ipiv,info)
if (info.ne.0) goto 10
call zgetri(lmmaxo,zfshto,lmmaxo,ipiv,zwork,lmmaxo,info)
if (info.ne.0) goto 10
!--------------------------------!
!     SHT matrices for lmaxi     !
!--------------------------------!
! allocate real SHT matrices
if (allocated(rbshti)) deallocate(rbshti)
allocate(rbshti(lmmaxi,lmmaxi))
if (allocated(rfshti)) deallocate(rfshti)
allocate(rfshti(lmmaxi,lmmaxi))
! allocate complex SHT matrices
if (allocated(zbshti)) deallocate(zbshti)
allocate(zbshti(lmmaxi,lmmaxi))
if (allocated(zfshti)) deallocate(zfshti)
allocate(zfshti(lmmaxi,lmmaxi))
! generate spherical covering set for lmaxi
call sphcover(lmmaxi,tp)
! generate real and complex spherical harmonics and set the backward SHT arrays
do itp=1,lmmaxi
  call genrlm(lmaxi,tp(:,itp),rlm)
  rbshti(itp,1:lmmaxi)=rlm(1:lmmaxi)
  call genylm(lmaxi,tp(:,itp),ylm)
  zbshti(itp,1:lmmaxi)=ylm(1:lmmaxi)
end do
! find the forward SHT arrays
! real
rfshti(:,:)=rbshti(:,:)
call dgetrf(lmmaxi,lmmaxi,rfshti,lmmaxi,ipiv,info)
if (info.ne.0) goto 10
call dgetri(lmmaxi,rfshti,lmmaxi,ipiv,work,lmmaxi,info)
if (info.ne.0) goto 10
! complex
zfshti(:,:)=zbshti(:,:)
call zgetrf(lmmaxi,lmmaxi,zfshti,lmmaxi,ipiv,info)
if (info.ne.0) goto 10
call zgetri(lmmaxi,zfshti,lmmaxi,ipiv,zwork,lmmaxi,info)
if (info.ne.0) goto 10
return
10 continue
write(*,*)
write(*,'("Error(genshtmat): unable to find inverse spherical harmonic &
 &transform")')
write(*,'(" => improper spherical covering")')
write(*,*)
stop
end subroutine
!EOC

