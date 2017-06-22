
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: readstate
! !INTERFACE:
subroutine readstate
! !USES:
use modmain
use moddftu
! !DESCRIPTION:
!   Reads in the charge density and other relevant variables from the file
!   {\tt STATE.OUT}. Checks for version and parameter compatibility.
!
! !REVISION HISTORY:
!   Created May 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
logical spinpol_
integer iostat
integer is,ias,lmmax,lm,ir,jr
integer idm0,idm1,idm,jdm,mapidm(3)
integer i1,i2,i3,j1,j2,j3
integer version_(3)
integer nspecies_,natoms_,lmmaxvr_
integer nrmt_(maxspecies),nrmtmax_
integer nrcmt_(maxspecies),nrcmtmax_
integer ngridg_(3),ngtot_,ngvec_
integer ndmag_,nspinor_,fsmtype_,ftmtype_
integer dftu_,lmmaxdm_
real(8) t1
! allocatable arrays
integer, allocatable :: mapir(:)
real(8), allocatable :: rsp_(:,:),rcmt_(:,:)
real(8), allocatable :: rhomt_(:,:,:),rhoir_(:)
real(8), allocatable :: vclmt_(:,:,:),vclir_(:)
real(8), allocatable :: vxcmt_(:,:,:),vxcir_(:)
real(8), allocatable :: vsmt_(:,:,:),vsir_(:)
real(8), allocatable :: magmt_(:,:,:,:),magir_(:,:)
real(8), allocatable :: bxcmt_(:,:,:,:),bxcir_(:,:)
real(8), allocatable :: bsmt_(:,:,:,:),bsir_(:,:)
real(8), allocatable :: bfsmcmt_(:,:)
complex(8), allocatable :: vsig_(:)
complex(8), allocatable :: vmatmt_(:,:,:,:,:),vmftm_(:,:,:,:,:)
open(50,file='STATE'//trim(filext),action='READ',form='UNFORMATTED', &
 status='OLD',iostat=iostat)
if (iostat.ne.0) then
  write(*,*)
  write(*,'("Error(readstate): error opening ",A)') 'STATE'//trim(filext)
  write(*,*)
  stop
end if
read(50) version_
if (version_(1).lt.2) then
  write(*,*)
  write(*,'("Error(readstate): unable to read STATE.OUT from versions earlier &
   &than 2.0.0")')
  write(*,*)
  stop
end if
if ((version(1).ne.version_(1)).or.(version(2).ne.version_(2)).or. &
    (version(3).ne.version_(3))) then
  write(*,*)
  write(*,'("Warning(readstate): different versions")')
  write(*,'(" current   : ",I3.3,".",I3.3,".",I3.3)') version
  write(*,'(" STATE.OUT : ",I3.3,".",I3.3,".",I3.3)') version_
end if
read(50) spinpol_
read(50) nspecies_
if (nspecies.ne.nspecies_) then
  write(*,*)
  write(*,'("Error(readstate): differing nspecies")')
  write(*,'(" current   : ",I4)') nspecies
  write(*,'(" STATE.OUT : ",I4)') nspecies_
  write(*,*)
  stop
end if
read(50) lmmaxvr_
read(50) nrmtmax_
read(50) nrcmtmax_
allocate(rsp_(nrmtmax_,nspecies))
allocate(rcmt_(nrcmtmax_,nspecies))
do is=1,nspecies
  read(50) natoms_
  if (natoms(is).ne.natoms_) then
    write(*,*)
    write(*,'("Error(readstate): differing natoms for species ",I4)') is
    write(*,'(" current   : ",I4)') natoms(is)
    write(*,'(" STATE.OUT : ",I4)') natoms_
    write(*,*)
    stop
  end if
  read(50) nrmt_(is)
  read(50) rsp_(1:nrmt_(is),is)
  read(50) nrcmt_(is)
  read(50) rcmt_(1:nrcmt_(is),is)
end do
read(50) ngridg_
read(50) ngvec_
read(50) ndmag_
if ((spinpol_).and.(ndmag_.ne.1).and.(ndmag_.ne.3)) then
  write(*,*)
  write(*,'("Error(readstate): invalid ndmag in STATE.OUT : ",I8)') ndmag_
  write(*,*)
  stop
end if
read(50) nspinor_
read(50) fsmtype_
if ((version_(1).gt.2).or.(version_(2).ge.3)) then
  read(50) ftmtype_
else
  ftmtype_=0
end if
read(50) dftu_
read(50) lmmaxdm_
ngtot_=ngridg_(1)*ngridg_(2)*ngridg_(3)
allocate(mapir(ngtot))
allocate(rhomt_(lmmaxvr_,nrmtmax_,natmtot))
allocate(rhoir_(ngtot_))
allocate(vclmt_(lmmaxvr_,nrmtmax_,natmtot))
allocate(vclir_(ngtot_))
allocate(vxcmt_(lmmaxvr_,nrmtmax_,natmtot))
allocate(vxcir_(ngtot_))
allocate(vsmt_(lmmaxvr_,nrmtmax_,natmtot))
allocate(vsir_(ngtot_))
! read the muffin-tin density
read(50) rhomt_,rhoir_
! read the Coulomb potential (spin independent)
read(50) vclmt_,vclir_
! read the exchange-correlation potential
read(50) vxcmt_,vxcir_
! read the Kohn-Sham effective potential
if ((version_(1).gt.2).or.(version_(2).ge.2)) then
  read(50) vsmt_,vsir_
else
  allocate(vsig_(ngvec_))
  read(50) vsmt_,vsir_,vsig_
  deallocate(vsig_)
end if
! read the magnetisation, exchange-correlation and effective magnetic fields
if (spinpol_) then
  allocate(magmt_(lmmaxvr_,nrmtmax_,natmtot,ndmag_))
  allocate(magir_(ngtot_,ndmag_))
  allocate(bxcmt_(lmmaxvr_,nrmtmax_,natmtot,ndmag_))
  allocate(bxcir_(ngtot_,ndmag_))
  allocate(bsmt_(lmmaxvr_,nrcmtmax_,natmtot,ndmag_))
  allocate(bsir_(ngtot_,ndmag_))
  read(50) magmt_,magir_
  read(50) bxcmt_,bxcir_
  read(50) bsmt_,bsir_
! read fixed spin moment effective fields
  if (fsmtype_.ne.0) then
    allocate(bfsmcmt_(3,natmtot))
    read(50) bfsmc
    read(50) bfsmcmt_
    if (fsmtype.ne.0) bfsmcmt(:,:)=bfsmcmt_(:,:)
    deallocate(bfsmcmt_)
  end if
end if
! read DFT+U potential matrix in each muffin-tin
if (((dftu.ne.0).and.(dftu_.ne.0)).or. &
    ((ftmtype.ne.0).and.(ftmtype_.ne.0))) then
  allocate(vmatmt_(lmmaxdm_,nspinor_,lmmaxdm_,nspinor_,natmtot))
  read(50) vmatmt_
  lmmax=min(lmmaxdm,lmmaxdm_)
  vmatmt(:,:,:,:,:)=0.d0
  if (nspinor.eq.nspinor_) then
    vmatmt(1:lmmax,:,1:lmmax,:,:)=vmatmt_(1:lmmax,:,1:lmmax,:,:)
  else if ((nspinor.eq.1).and.(nspinor_.eq.2)) then
    vmatmt(1:lmmax,1,1:lmmax,1,:)=0.5d0*(vmatmt_(1:lmmax,1,1:lmmax,1,:) &
     +vmatmt_(1:lmmax,2,1:lmmax,2,:))
  else
    vmatmt(1:lmmax,1,1:lmmax,1,:)=vmatmt_(1:lmmax,1,1:lmmax,1,:)
    vmatmt(1:lmmax,2,1:lmmax,2,:)=vmatmt_(1:lmmax,1,1:lmmax,1,:)
  end if
  deallocate(vmatmt_)
end if
! read fixed tensor moment potential matrix elements
if ((ftmtype.ne.0).and.(ftmtype_.ne.0)) then
  allocate(vmftm_(lmmaxdm_,nspinor_,lmmaxdm_,nspinor_,natmtot))
  read(50) vmftm_
  lmmax=min(lmmaxdm,lmmaxdm_)
  vmftm_(:,:,:,:,:)=0.d0
  if (nspinor.eq.nspinor_) then
    vmftm(1:lmmax,:,1:lmmax,:,:)=vmftm_(1:lmmax,:,1:lmmax,:,:)
  else if ((nspinor.eq.1).and.(nspinor_.eq.2)) then
    vmftm(1:lmmax,1,1:lmmax,1,:)=0.5d0*(vmftm_(1:lmmax,1,1:lmmax,1,:) &
     +vmftm_(1:lmmax,2,1:lmmax,2,:))
  else
    vmftm(1:lmmax,1,1:lmmax,1,:)=vmftm_(1:lmmax,1,1:lmmax,1,:)
    vmftm(1:lmmax,2,1:lmmax,2,:)=vmftm_(1:lmmax,1,1:lmmax,1,:)
  end if
  deallocate(vmftm_)
end if
close(50)
! component map for spin-polarised case
if (spinpol) then
  if (ndmag.eq.ndmag_) then
    idm0=1; idm1=ndmag
    do idm=1,ndmag
      mapidm(idm)=idm
    end do
  else
    idm0=ndmag; idm1=ndmag
    mapidm(ndmag)=ndmag_
  end if
end if
!---------------------------!
!     muffin-tin arrays     !
!---------------------------!
rhomt(:,:,:)=0.d0
vclmt(:,:,:)=0.d0
vxcmt(:,:,:)=0.d0
vsmt(:,:,:)=0.d0
if (spinpol) then
  magmt(:,:,:,:)=0.d0
  bxcmt(:,:,:,:)=0.d0
  bsmt(:,:,:,:)=0.d0
end if
lmmax=min(lmmaxvr,lmmaxvr_)
! interpolate the old arrays on the new radial mesh
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(is,lm,idm,jdm)
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  do lm=1,lmmax
    call rfinterp(nrmt_(is),rsp_(:,is),lmmaxvr_,rhomt_(lm,1,ias),nrmt(is), &
     rsp(:,is),lmmaxvr,rhomt(lm,1,ias))
  end do
  do lm=1,lmmax
    call rfinterp(nrmt_(is),rsp_(:,is),lmmaxvr_,vclmt_(lm,1,ias),nrmt(is), &
     rsp(:,is),lmmaxvr,vclmt(lm,1,ias))
  end do
  do lm=1,lmmax
    call rfinterp(nrmt_(is),rsp_(:,is),lmmaxvr_,vxcmt_(lm,1,ias),nrmt(is), &
     rsp(:,is),lmmaxvr,vxcmt(lm,1,ias))
  end do
  do lm=1,lmmax
    call rfinterp(nrmt_(is),rsp_(:,is),lmmaxvr_,vsmt_(lm,1,ias),nrmt(is), &
     rsp(:,is),lmmaxvr,vsmt(lm,1,ias))
  end do
  if ((spinpol).and.(spinpol_)) then
    do idm=idm0,idm1
      jdm=mapidm(idm)
      do lm=1,lmmax
        call rfinterp(nrmt_(is),rsp_(:,is),lmmaxvr_,magmt_(lm,1,ias,jdm), &
         nrmt(is),rsp(:,is),lmmaxvr,magmt(lm,1,ias,idm))
      end do
      do lm=1,lmmax
        call rfinterp(nrmt_(is),rsp_(:,is),lmmaxvr_,bxcmt_(lm,1,ias,jdm), &
         nrmt(is),rsp(:,is),lmmaxvr,bxcmt(lm,1,ias,idm))
      end do
      do lm=1,lmmax
        call rfinterp(nrcmt_(is),rcmt_(:,is),lmmaxvr_,bsmt_(lm,1,ias,jdm), &
         nrcmt(is),rcmt(:,is),lmmaxvr,bsmt(lm,1,ias,idm))
      end do
    end do
  end if
end do
!$OMP END DO
!$OMP END PARALLEL
!-----------------------------!
!     interstitial arrays     !
!-----------------------------!
rhoir(:)=0.d0
vclir(:)=0.d0
vxcir(:)=0.d0
vsir(:)=0.d0
if (spinpol) then
  magir(:,:)=0.d0
  bxcir(:,:)=0.d0
  bsir(:,:)=0.d0
end if
! map from new grid to old
ir=0
do i3=0,ngridg(3)-1
  t1=dble(i3*ngridg_(3))/dble(ngridg(3))
  j3=modulo(nint(t1),ngridg_(3))
  do i2=0,ngridg(2)-1
    t1=dble(i2*ngridg_(2))/dble(ngridg(2))
    j2=modulo(nint(t1),ngridg_(2))
    do i1=0,ngridg(1)-1
      t1=dble(i1*ngridg_(1))/dble(ngridg(1))
      j1=modulo(nint(t1),ngridg_(1))
      ir=ir+1
      jr=j3*ngridg_(2)*ngridg_(1)+j2*ngridg_(1)+j1+1
      mapir(ir)=jr
    end do
  end do
end do
do ir=1,ngtot
  jr=mapir(ir)
  rhoir(ir)=rhoir_(jr)
  vclir(ir)=vclir_(jr)
  vxcir(ir)=vxcir_(jr)
  vsir(ir)=vsir_(jr)
end do
if ((spinpol).and.(spinpol_)) then
  do idm=idm0,idm1
    jdm=mapidm(idm)
    do ir=1,ngtot
      jr=mapir(ir)
      magir(ir,idm)=magir_(jr,jdm)
    end do
    do ir=1,ngtot
      jr=mapir(ir)
      bxcir(ir,idm)=bxcir_(jr,jdm)
    end do
    do ir=1,ngtot
      jr=mapir(ir)
      bsir(ir,idm)=bsir_(jr,jdm)
    end do
  end do
end if
deallocate(mapir,rsp_,rcmt_,rhomt_,rhoir_,vclmt_,vclir_)
deallocate(vxcmt_,vxcir_,vsmt_,vsir_)
if (spinpol_) deallocate(magmt_,magir_,bxcmt_,bxcir_,bsmt_,bsir_)
! Fourier transform Kohn-Sham potential to G-space
call genvsig
return
end subroutine
!EOC
