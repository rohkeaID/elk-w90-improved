
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine projsbf
use modmain
implicit none
! local variables
integer idm,is,ias,np
real(8) t1
! allocatable arrays
real(8), allocatable :: rfmt(:,:),rfir(:)
real(8), allocatable :: grfmt(:,:,:),grfir(:,:)
complex(8), allocatable :: zrhomt(:,:),zrhoir(:)
complex(8), allocatable :: zvclmt(:,:),zvclir(:)
allocate(rfmt(npmtmax,natmtot),rfir(ngtot))
allocate(grfmt(npmtmax,natmtot,3),grfir(ngtot,3))
allocate(zrhomt(npmtmax,natmtot),zrhoir(ngtot))
allocate(zvclmt(npmtmax,natmtot),zvclir(ngtot))
if (.not.spinpol) then
  write(*,*)
  write(*,'("Error(projsbf): spin-unpolarised field is zero")')
  write(*,*)
  stop
end if
! compute the divergence of B_xc
rfmt(:,:)=0.d0
rfir(:)=0.d0
do idm=1,3
  call gradrf(bxcmt(:,:,idm),bxcir(:,idm),grfmt,grfir)
  do ias=1,natmtot
    is=idxis(ias)
    np=npmt(is)
    rfmt(1:np,ias)=rfmt(1:np,ias)+grfmt(1:np,ias,idm)
  end do
  rfir(:)=rfir(:)+grfir(:,idm)
end do
! convert real muffin-tin divergence to complex spherical harmonic expansion
do ias=1,natmtot
  is=idxis(ias)
  call rtozfmt(nrmt(is),nrmti(is),rfmt(:,ias),zrhomt(:,ias))
end do
! store real interstitial divergence in a complex array
zrhoir(:)=rfir(:)
! solve the complex Poisson's equation
call genzvclmt(nrmt,nrmti,nrspmax,rsp,npmtmax,zrhomt,zvclmt)
call zpotcoul(nrmt,nrmti,npmt,npmti,nrspmax,rsp,ngridg,igfft,ngvec,gc,gclg, &
 ngvec,jlgrmt,ylmg,sfacg,zrhoir,npmtmax,zvclmt,zvclir)
! convert complex muffin-tin potential to real spherical harmonic expansion
do ias=1,natmtot
  is=idxis(ias)
  call ztorfmt(nrmt(is),nrmti(is),zvclmt(:,ias),rfmt(:,ias))
end do
! store complex interstitial potential in real array
rfir(:)=dble(zvclir(:))
! compute the gradient
call gradrf(rfmt,rfir,grfmt,grfir)
! add gradient over 4*pi to existing B_xc
t1=1.d0/fourpi
bxcmt(:,:,:)=bxcmt(:,:,:)+t1*grfmt(:,:,:)
bxcir(:,:)=bxcir(:,:)+t1*grfir(:,:)
deallocate(rfmt,rfir,grfmt,grfir)
deallocate(zrhomt,zrhoir,zvclmt,zvclir)
return
end subroutine

