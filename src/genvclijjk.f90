
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genvclijjk
! !INTERFACE:
subroutine genvclijjk(ikp,vclijjk)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ikp     : k-point from non-reduced set (in,integer)
!   vclijjk : Coulomb matrix elements (out,complex(nstsv,nstsv,nstsv,nkpt))
! !DESCRIPTION:
!   Calculates Coulomb matrix elements of the type $(i-jj-k)$.
!
! !REVISION HISTORY:
!   Created 2008 (Sharma)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ikp
complex(8), intent(out) :: vclijjk(nstsv,nstsv,nstsv,nkpt)
! local variables
integer ik,iv(3)
integer ig,iq,igq0
integer ist1,ist2,ist3
real(8) cfq,v(3),t1
complex(8) zrho01,zrho02,z1,z2
complex(8) sfacgq0(natmtot)
! automatic arrays
integer idx(nstsv)
! allocatable arrays
real(8), allocatable :: vgqc(:,:),tpgqc(:,:),gqc(:)
real(8), allocatable :: jlgqr(:,:,:),jlgq0r(:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: evecfv(:,:),evecsv(:,:)
complex(8), allocatable :: ylmgq(:,:),sfacgq(:,:)
complex(8), allocatable :: wfmt1(:,:,:,:,:),wfmt2(:,:,:,:,:)
complex(8), allocatable :: wfir1(:,:,:),wfir2(:,:,:)
complex(8), allocatable :: zrhomt(:,:,:,:),zrhoir(:,:)
complex(8), allocatable :: zvclmt(:,:,:),zvclir(:)
! external functions
complex(8) zfinp
external zfinp
! allocate local arrays
allocate(vgqc(3,ngvec),tpgqc(2,ngvec),gqc(ngvec))
allocate(jlgqr(0:lnpsd,ngvec,nspecies),jlgq0r(0:lmaxvr,nrcmtmax,nspecies))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv(nmatmax,nstfv),evecsv(nstsv,nstsv))
allocate(ylmgq(lmmaxvr,ngvec),sfacgq(ngvec,natmtot))
allocate(wfmt1(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfmt2(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfir1(ngtot,nspinor,nstsv),wfir2(ngtot,nspinor,nstsv))
allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot,nstsv),zrhoir(ngtot,nstsv))
allocate(zvclmt(lmmaxvr,nrcmtmax,natmtot),zvclir(ngtot))
! factor for long-range term
cfq=0.5d0*(omega/pi)**2
! get the eigenvectors from file for non-reduced k-point ikp
call getevecfv(filext,vkl(:,ikp),vgkl(:,:,1,ikp),evecfv)
call getevecsv(filext,vkl(:,ikp),evecsv)
! find the matching coefficients
call match(ngk(1,ikp),gkc(:,1,ikp),tpgkc(:,:,1,ikp),sfacgk(:,:,1,ikp),apwalm)
! index to all states
do ist1=1,nstsv
  idx(ist1)=ist1
end do
! calculate the wavefunctions for all states of passed non-reduced k-point ikp
call genwfsv(.false.,.false.,nstsv,idx,ngk(1,ikp),igkig(:,1,ikp),apwalm, &
 evecfv,evecsv,wfmt2,ngtot,wfir2)
! start loop over reduced k-point set
do ik=1,nkpt
! get the eigenvectors from file
  call getevecfv(filext,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
  call getevecsv(filext,vkl(:,ik),evecsv)
! find the matching coefficients
  call match(ngk(1,ik),gkc(:,1,ik),tpgkc(:,:,1,ik),sfacgk(:,:,1,ik),apwalm)
! calculate the wavefunctions for all states of the reduced k-point
  call genwfsv(.false.,.false.,nstsv,idx,ngk(1,ik),igkig(:,1,ik),apwalm, &
   evecfv,evecsv,wfmt1,ngtot,wfir1)
! determine q-vector
  iv(:)=ivk(:,ik)-ivk(:,ikp)
  iv(:)=modulo(iv(:),ngridk(:))
  iq=iqmap(iv(1),iv(2),iv(3))
  v(:)=vkc(:,ik)-vkc(:,ikp)
  do ig=1,ngvec
! determine G+q-vectors
    vgqc(:,ig)=vgc(:,ig)+v(:)
! G+q-vector length and (theta, phi) coordinates
    call sphcrd(vgqc(:,ig),gqc(ig),tpgqc(:,ig))
! spherical harmonics for G+q-vectors
    call genylm(lmaxvr,tpgqc(:,ig),ylmgq(:,ig))
  end do
! structure factors for G+q
  call gensfacgp(ngvec,vgqc,ngvec,sfacgq)
! find the shortest G+q-vector
  call findigp0(ngvec,gqc,igq0)
  sfacgq0(:)=sfacgq(igq0,:)
! compute the required spherical Bessel functions
  call genjlgpr(lnpsd,gqc,jlgqr)
  call genjlgq0r(gqc(igq0),jlgq0r)
!----------------------------------------------!
!     valence-valence-valence contribution     !
!----------------------------------------------!
  do ist2=1,nstsv
    do ist1=1,nstsv
! calculate the complex overlap density for all states
      call genzrho(.true.,.true.,wfmt2(:,:,:,:,ist2),wfir2(:,:,ist2), &
       wfmt1(:,:,:,:,ist1),wfir1(:,:,ist1),zrhomt(:,:,:,ist1),zrhoir(:,ist1))
    end do
    do ist1=1,nstsv
! compute the potential and G=0 coefficient of the density
      call genzvclmt(nrcmt,nrcmtinr,nrcmtmax,rcmt,nrcmtmax,zrhomt(:,:,:,ist1), &
       zvclmt)
      call zpotcoul(nrcmt,nrcmtinr,nrcmtmax,rcmt,igq0,gqc,jlgqr,ylmgq,sfacgq, &
       zrhoir(:,ist1),nrcmtmax,zvclmt,zvclir,zrho02)
      z1=zfinp(zrhomt(:,:,:,ist1),zrhoir(:,ist1),zvclmt,zvclir)
      t1=cfq*wiq2(iq)*(dble(zrho02)**2+aimag(zrho02)**2)
      vclijjk(ist1,ist1,ist2,ik)=wkptnr*dble(z1)+t1
      do ist3=ist1+1,nstsv
        z1=zfinp(zrhomt(:,:,:,ist3),zrhoir(:,ist3),zvclmt,zvclir)
! compute the density coefficient of the smallest G+q-vector
        call zrhogp(jlgq0r,ylmgq(:,igq0),sfacgq0,zrhomt(:,:,:,ist3), &
         zrhoir(:,ist3),zrho01)
        z2=cfq*wiq2(iq)*conjg(zrho01)*zrho02
        vclijjk(ist3,ist1,ist2,ik)=wkptnr*z1+z2
      end do
    end do
  end do
! calculate the lower diagonal
  do ist1=1,nstsv
    do ist3=1,ist1-1
      vclijjk(ist3,ist1,:,ik)=conjg(vclijjk(ist1,ist3,:,ik))
    end do
  end do
! end loop over reduced k-point set
end do
deallocate(vgqc,tpgqc,gqc,jlgqr,jlgq0r)
deallocate(apwalm,evecfv,evecsv,ylmgq,sfacgq)
deallocate(wfmt1,wfmt2,wfir1,wfir2)
deallocate(zrhomt,zrhoir,zvclmt,zvclir)
return
end subroutine
!EOC

