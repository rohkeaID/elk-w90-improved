
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine projsbf
use modmain
implicit none
! local variables
integer is,ias,idm
real(8) t1
complex(8) zrho0
! allocatable arrays
real(8), allocatable :: rfmt(:,:,:),rfir(:)
real(8), allocatable :: grfmt(:,:,:,:),grfir(:,:)
complex(8), allocatable :: zrhomt(:,:,:),zrhoir(:)
complex(8), allocatable :: zvclmt(:,:,:),zvclir(:)
!allocate(rvfmt(lmmaxvr,nrmtmax,natmtot,3),rvfir(ngtot,3))
allocate(rfmt(lmmaxvr,nrmtmax,natmtot),rfir(ngtot))
allocate(grfmt(lmmaxvr,nrmtmax,natmtot,3),grfir(ngtot,3))
allocate(zrhomt(lmmaxvr,nrmtmax,natmtot),zrhoir(ngtot))
allocate(zvclmt(lmmaxvr,nrmtmax,natmtot),zvclir(ngtot))
if (.not.spinpol) then
  write(*,*)
  write(*,'("Error(projsbf): spin-unpolarised field is zero")')
  write(*,*)
  stop
end if
! compute the divergence of B_xc
rfmt(:,:,:)=0.d0
rfir(:)=0.d0
do idm=1,3
  call gradrf(bxcmt(:,:,:,idm),bxcir(:,idm),grfmt,grfir)
  do ias=1,natmtot
    is=idxis(ias)
    call rfmtadd(nrmt(is),nrmtinr(is),1,grfmt(:,:,ias,idm),rfmt(:,:,ias))
  end do
  rfir(:)=rfir(:)+grfir(:,idm)
end do
! convert real muffin-tin divergence to complex spherical harmonic expansion
do ias=1,natmtot
  is=idxis(ias)
  call rtozfmt(nrmt(is),nrmtinr(is),1,rfmt(:,:,ias),1,zrhomt(:,:,ias))
end do
! store real interstitial divergence in a complex array
zrhoir(:)=rfir(:)
! solve the complex Poisson's equation
call genzvclmt(nrmt,nrmtinr,nrspmax,rsp,nrmtmax,zrhomt,zvclmt)
call zpotcoul(nrmt,nrmtinr,nrspmax,rsp,1,gc,jlgr,ylmg,sfacg,zrhoir,nrmtmax, &
 zvclmt,zvclir,zrho0)
! convert complex muffin-tin potential to real spherical harmonic expansion
do ias=1,natmtot
  is=idxis(ias)
  call ztorfmt(nrmt(is),nrmtinr(is),1,zvclmt(:,:,ias),1,rfmt(:,:,ias))
end do
! store complex interstitial potential in real array
rfir(:)=dble(zvclir(:))
! compute the gradient
call gradrf(rfmt,rfir,grfmt,grfir)
! add gradient over 4*pi to existing B_xc
t1=1.d0/fourpi
bxcmt(:,:,:,:)=bxcmt(:,:,:,:)+t1*grfmt(:,:,:,:)
bxcir(:,:)=bxcir(:,:)+t1*grfir(:,:)
deallocate(rfmt,rfir,grfmt,grfir)
deallocate(zrhomt,zrhoir,zvclmt,zvclir)
return
end subroutine

