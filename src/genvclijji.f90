
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genvclijji
! !INTERFACE:
subroutine genvclijji(ikp,vclijji)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ikp     : k-point from non-reduced set (in,integer)
!   vclijji : Coulomb matrix elements (out,real(nstsv,nstsv,nkpt))
! !DESCRIPTION:
!   Calculates the Coulomb matrix elements of the type $(i-jj-i)$.
!
! !REVISION HISTORY:
!   Created June 2008 (Sharma)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ikp
real(8), intent(out) :: vclijji(nstsv,nstsv,nkpt)
! local variables
integer ik,iv(3)
integer ig,iq,igq0
integer ist1,ist2
real(8) cfq,v(3),tp(2),t1
complex(8) zgq0,z1
! automatic arrays
integer idx(nstsv)
! allocatable arrays
real(8), allocatable :: vgqc(:,:),gqc(:),jlgqrmt(:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: evecfv(:,:),evecsv(:,:)
complex(8), allocatable :: ylmgq(:,:),sfacgq(:,:)
complex(8), allocatable :: wfmt2(:,:,:,:),wfir2(:,:,:)
complex(8), allocatable :: wfmt1(:,:,:,:),wfir1(:,:,:)
complex(8), allocatable :: zrhomt(:,:),zrhoir(:)
complex(8), allocatable :: zvclmt(:,:),zvclir(:)
! external functions
complex(8) zfinp
external zfinp
! allocate local arrays
allocate(vgqc(3,ng2gk),gqc(ng2gk))
allocate(jlgqrmt(0:lnpsd,ng2gk,nspecies))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv(nmatmax,nstfv),evecsv(nstsv,nstsv))
allocate(ylmgq(lmmaxo,ng2gk),sfacgq(ng2gk,natmtot))
allocate(wfmt1(npcmtmax,natmtot,nspinor,nstsv),wfir1(ngtot,nspinor,nstsv))
allocate(wfmt2(npcmtmax,natmtot,nspinor,nstsv),wfir2(ngtot,nspinor,nstsv))
allocate(zrhomt(npcmtmax,natmtot),zrhoir(ngtot))
allocate(zvclmt(npcmtmax,natmtot),zvclir(ngtot))
! factor for long-range term
cfq=0.5d0*(omega/pi)**2
! get the eigenvectors from file for non-reduced k-point ikp
call getevecfv(filext,0,vkl(:,ikp),vgkl(:,:,1,ikp),evecfv)
call getevecsv(filext,0,vkl(:,ikp),evecsv)
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
  call getevecfv(filext,ik,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
  call getevecsv(filext,ik,vkl(:,ik),evecsv)
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
  do ig=1,ng2gk
! determine G+q-vectors
    vgqc(:,ig)=vgc(:,ig)+v(:)
! G+q-vector length and (theta, phi) coordinates
    call sphcrd(vgqc(:,ig),gqc(ig),tp)
! spherical harmonics for G+q-vectors
    call genylm(lmaxo,tp,ylmgq(:,ig))
  end do
! structure factors for G+q
  call gensfacgp(ng2gk,vgqc,ng2gk,sfacgq)
! find the shortest G+q-vector
  call findigp0(ng2gk,gqc,igq0)
! compute the required spherical Bessel functions
  call genjlgprmt(lnpsd,ng2gk,gqc,ng2gk,jlgqrmt)
!----------------------------------------------!
!     valence-valence-valence contribution     !
!----------------------------------------------!
  do ist1=1,nstsv
    do ist2=1,nstsv
! calculate the complex overlap density
      call genzrho(.true.,.true.,wfmt2(:,:,:,ist2),wfir2(:,:,ist2), &
       wfmt1(:,:,:,ist1),wfir1(:,:,ist1),zrhomt,zrhoir)
! compute the potential and G=0 coefficient of the density
      call genzvclmt(nrcmt,nrcmti,nrcmtmax,rcmt,npcmtmax,zrhomt,zvclmt)
      call zpotcoul(nrcmt,nrcmti,npcmt,npcmti,nrcmtmax,rcmt,ng2gk,igq0,gqc, &
       ng2gk,jlgqrmt,ylmgq,sfacgq,zrhoir,npcmtmax,zvclmt,zvclir,zgq0)
      z1=zfinp(zrhomt,zrhoir,zvclmt,zvclir)
      t1=cfq*wiq2(iq)*(dble(zgq0)**2+aimag(zgq0)**2)
      vclijji(ist1,ist2,ik)=wkptnr*dble(z1)+t1
! end loop over ist2
    end do
! end loop over ist1
  end do
! end loop over reduced k-point set
end do
deallocate(vgqc,gqc,jlgqrmt)
deallocate(apwalm,evecfv,evecsv,ylmgq,sfacgq)
deallocate(wfmt1,wfmt2,wfir1,wfir2)
deallocate(zrhomt,zrhoir,zvclmt,zvclir)
return
end subroutine
!EOC

