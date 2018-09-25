
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
integer is,ia,ias,lmmax,lm,ir,jr
integer idm,jdm,mapidm(3)
integer i1,i2,i3,j1,j2,j3,n
integer version_(3)
integer nspecies_,natoms_,lmmaxo_
integer nrmt_(maxspecies),nrmtmax_
integer nrcmt_(maxspecies),nrcmtmax_
integer ngridg_(3),ngtot_,ngvec_
integer ndmag_,nspinor_,fsmtype_,ftmtype_
integer dftu_,lmmaxdm_,xcgrad_
real(8) t1
! allocatable arrays
integer, allocatable :: mapir(:)
real(8), allocatable :: rsp_(:,:),rcmt_(:,:)
real(8), allocatable :: rfmt_(:,:,:),rfir_(:)
real(8), allocatable :: rvfmt_(:,:,:,:),rvfir_(:,:)
real(8), allocatable :: rvfcmt_(:,:,:,:),rfmt(:,:)
real(8), allocatable :: bfsmcmt_(:,:),fi(:),fo(:)
complex(8), allocatable :: vsig_(:)
complex(8), allocatable :: vmatmt_(:,:,:,:,:),vmftm_(:,:,:,:,:)
open(40,file='STATE'//trim(filext),form='UNFORMATTED',status='OLD', &
 iostat=iostat)
if (iostat.ne.0) then
  write(*,*)
  write(*,'("Error(readstate): error opening ",A)') 'STATE'//trim(filext)
  write(*,*)
  stop
end if
read(40) version_
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
read(40) spinpol_
read(40) nspecies_
if (nspecies.ne.nspecies_) then
  write(*,*)
  write(*,'("Error(readstate): differing nspecies")')
  write(*,'(" current   : ",I4)') nspecies
  write(*,'(" STATE.OUT : ",I4)') nspecies_
  write(*,*)
  stop
end if
read(40) lmmaxo_
lmmax=min(lmmaxo,lmmaxo_)
read(40) nrmtmax_
read(40) nrcmtmax_
allocate(rsp_(nrmtmax_,nspecies))
allocate(rcmt_(nrcmtmax_,nspecies))
do is=1,nspecies
  read(40) natoms_
  if (natoms(is).ne.natoms_) then
    write(*,*)
    write(*,'("Error(readstate): differing natoms for species ",I4)') is
    write(*,'(" current   : ",I4)') natoms(is)
    write(*,'(" STATE.OUT : ",I4)') natoms_
    write(*,*)
    stop
  end if
  read(40) nrmt_(is)
  read(40) rsp_(1:nrmt_(is),is)
  read(40) nrcmt_(is)
  read(40) rcmt_(1:nrcmt_(is),is)
end do
read(40) ngridg_
read(40) ngvec_
read(40) ndmag_
if ((spinpol_).and.(ndmag_.ne.1).and.(ndmag_.ne.3)) then
  write(*,*)
  write(*,'("Error(readstate): invalid ndmag in STATE.OUT : ",I8)') ndmag_
  write(*,*)
  stop
end if
read(40) nspinor_
read(40) fsmtype_
if ((version_(1).gt.2).or.(version_(2).ge.3)) then
  read(40) ftmtype_
else
  ftmtype_=0
end if
read(40) dftu_
read(40) lmmaxdm_
if ((version_(1).gt.5).or.((version_(1).eq.5).and.(version_(2).ge.1))) then
  read(40) xcgrad_
else
  xcgrad_=0
end if
ngtot_=ngridg_(1)*ngridg_(2)*ngridg_(3)
! map from old interstitial grid to new
allocate(mapir(ngtot))
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
allocate(rfmt_(lmmaxo_,nrmtmax_,natmtot),rfir_(ngtot_))
allocate(rfmt(lmmaxo,nrmtmax))
n=max(nrmtmax,nrmtmax_)
allocate(fi(n),fo(n))
! read the muffin-tin density
read(40) rfmt_,rfir_
! regrid and pack the muffin-tin function
call rgfmt(rhomt)
! regrid the interstitial function
rhoir(:)=rfir_(mapir(:))
! read the Coulomb potential, regrid and pack
read(40) rfmt_,rfir_
call rgfmt(vclmt)
vclir(:)=rfir_(mapir(:))
! read the exchange-correlation potential, regrid and pack
read(40) rfmt_,rfir_
call rgfmt(vxcmt)
vxcir(:)=rfir_(mapir(:))
! read the Kohn-Sham effective potential, regrid and pack
if ((version_(1).gt.2).or.(version_(2).ge.2)) then
  read(40) rfmt_,rfir_
else
  allocate(vsig_(ngvec_))
  read(40) rfmt_,rfir_,vsig_
  deallocate(vsig_)
end if
call rgfmt(vsmt)
vsir(:)=rfir_(mapir(:))
! read the magnetisation, exchange-correlation and effective magnetic fields
if (spinpol_) then
! component map for spin-polarised case
  mapidm(:)=0
  if (ndmag.eq.ndmag_) then
    do idm=1,ndmag
      mapidm(idm)=idm
    end do
  else
    mapidm(ndmag)=ndmag_
  end if
  allocate(rvfmt_(lmmaxo_,nrmtmax_,natmtot,ndmag_))
  allocate(rvfir_(ngtot_,ndmag_))
  allocate(rvfcmt_(lmmaxo_,nrcmtmax_,natmtot,ndmag_))
  read(40) rvfmt_,rvfir_
  call rgvfmt(magmt)
  call rgvir(magir)
  read(40) rvfmt_,rvfir_
  call rgvfmt(bxcmt)
  call rgvir(bxcir)
  read(40) rvfcmt_,rvfir_
  call rgvfcmt(bsmt)
  call rgvir(bsir)
  deallocate(rvfmt_,rvfir_,rvfcmt_)
! read fixed spin moment effective fields
  if (fsmtype_.ne.0) then
    allocate(bfsmcmt_(3,natmtot))
    read(40) bfsmc
    read(40) bfsmcmt_
    if (fsmtype.ne.0) bfsmcmt(:,:)=bfsmcmt_(:,:)
! make sure that the constraining fields are perpendicular to the fixed moments
! for fixed direction calculations (Y. Kvashnin and LN)
    if (fsmtype.lt.0) then
      if (ncmag) then
        call r3vo(momfix,bfsmc)
        do is=1,nspecies
          do ia=1,natoms(is)
            ias=idxas(ia,is)
            call r3vo(mommtfix(:,ia,is),bfsmcmt(:,ias))
          end do
        end do
      else
        bfsmc(:)=0.d0
        bfsmcmt(:,:)=0.d0
      end if
    end if
    deallocate(bfsmcmt_)
  end if
end if
if (xcgrad.eq.4) then
  if (xcgrad_.eq.4) then
    read(40) rfmt_,rfir_
    call rgfmt(wxcmt)
    wxcir(:)=rfir_(mapir(:))
  else
    wxcmt(:,:)=0.d0
    wxcir(:)=0.d0
  end if
  call genws
end if
deallocate(rfmt_,rfir_,rfmt,fi,fo)
! read DFT+U potential matrix in each muffin-tin
if (((dftu.ne.0).and.(dftu_.ne.0)).or. &
    ((ftmtype.ne.0).and.(ftmtype_.ne.0))) then
  allocate(vmatmt_(lmmaxdm_,nspinor_,lmmaxdm_,nspinor_,natmtot))
  read(40) vmatmt_
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
  read(40) vmftm_
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
close(40)
return

contains

subroutine rgfmt(rfmtp)
implicit none
! arguments
real(8), intent(out) :: rfmtp(npmtmax,natmtot)
do ias=1,natmtot
  is=idxis(ias)
! regrid the muffin-tin function
  do lm=1,lmmax
    fi(1:nrmt_(is))=rfmt_(lm,1:nrmt_(is),ias)
    call rfinterp(nrmt_(is),rsp_(:,is),fi,nrmt(is),rsp(:,is),fo)
    rfmt(lm,1:nrmt(is))=fo(1:nrmt(is))
  end do
  rfmt(lmmax+1:lmmaxo,1:nrmt(is))=0.d0
! pack the muffin-tin function
  call rfmtpack(.true.,nrmt(is),nrmti(is),rfmt,rfmtp(:,ias))
end do
return
end subroutine

subroutine rgvfmt(rvfmt)
implicit none
! arguments
real(8), intent(out) :: rvfmt(npmtmax,natmtot,ndmag)
do idm=1,ndmag
  jdm=mapidm(idm)
  if (jdm.eq.0) then
    rvfmt(:,:,idm)=0.d0
    cycle
  end if
  do ias=1,natmtot
    is=idxis(ias)
    do lm=1,lmmax
      fi(1:nrmt_(is))=rvfmt_(lm,1:nrmt_(is),ias,jdm)
      call rfinterp(nrmt_(is),rsp_(:,is),fi,nrmt(is),rsp(:,is),fo)
      rfmt(lm,1:nrmt(is))=fo(1:nrmt(is))
    end do
    rfmt(lmmax+1:lmmaxo,1:nrmt(is))=0.d0
    call rfmtpack(.true.,nrmt(is),nrmti(is),rfmt,rvfmt(:,ias,idm))
  end do
end do
return
end subroutine

subroutine rgvfcmt(rvfcmt)
implicit none
! arguments
real(8), intent(out) :: rvfcmt(npcmtmax,natmtot,ndmag)
do idm=1,ndmag
  jdm=mapidm(idm)
  if (jdm.eq.0) then
    rvfcmt(:,:,idm)=0.d0
    cycle
  end if
  do ias=1,natmtot
    is=idxis(ias)
    do lm=1,lmmax
      fi(1:nrcmt_(is))=rvfcmt_(lm,1:nrcmt_(is),ias,jdm)
      call rfinterp(nrcmt_(is),rcmt_(:,is),fi,nrcmt(is),rcmt(:,is),fo)
      rfmt(lm,1:nrcmt(is))=fo(1:nrcmt(is))
    end do
    rfmt(lmmax+1:lmmaxo,1:nrcmt(is))=0.d0
    call rfmtpack(.true.,nrcmt(is),nrcmti(is),rfmt,rvfcmt(:,ias,idm))
  end do
end do
return
end subroutine

subroutine rgvir(rvfir)
implicit none
! arguments
real(8), intent(out) :: rvfir(ngtot,ndmag)
do idm=1,ndmag
  jdm=mapidm(idm)
  if (jdm.eq.0) then
    rvfir(:,idm)=0.d0
    cycle
  end if
  rvfir(:,idm)=rvfir_(mapir(:),jdm)
end do
return
end subroutine

end subroutine
!EOC

