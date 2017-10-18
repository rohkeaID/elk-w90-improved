
! Copyright (C) 2006 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine eveqnhf(ikp,vmt,vir,bmt,bir,evecsvp)
use modmain
implicit none
! arguments
integer, intent(in) :: ikp
real(8), intent(in) :: vmt(npcmtmax,natmtot),vir(ngtot)
real(8), intent(in) :: bmt(npcmtmax,natmtot,ndmag),bir(ngtot,ndmag)
complex(8), intent(inout) :: evecsvp(nstsv,nstsv)
! local variables
integer ik,jk,ist1,ist2,ist3
integer iv(3),ig,iq,igq0
real(8) cfq,vc(3),tp(2),t1
complex(8) zgq01,zgq02,z1,z2
! automatic arrays
integer idx(nstsv)
complex(8) sfacgq0(natmtot)
! allocatable arrays
real(8), allocatable :: vgqc(:,:),gqc(:)
real(8), allocatable :: jlgqrmt(:,:,:),jlgq0r(:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:),evecfv(:,:),evecsv(:,:)
complex(8), allocatable :: h(:,:),v(:,:),kmat(:,:)
complex(8), allocatable :: ylmgq(:,:),sfacgq(:,:)
complex(8), allocatable :: wfmt1(:,:,:,:),wfir1(:,:,:)
complex(8), allocatable :: wfmt2(:,:,:,:),wfir2(:,:,:)
complex(8), allocatable :: zrhomt(:,:,:),zrhoir(:,:)
complex(8), allocatable :: zvclmt(:,:),zvclir(:)
! external functions
complex(8) zfinp
external zfinp
!$OMP CRITICAL
write(*,'("Info(eveqnhf): ",I6," of ",I6," k-points")') ikp,nkpt
!$OMP END CRITICAL
! allocate local arrays
allocate(vgqc(3,ng2gk),gqc(ng2gk))
allocate(jlgqrmt(0:lnpsd,ng2gk,nspecies),jlgq0r(0:lmaxo,nrcmtmax,nspecies))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv(nmatmax,nstfv),evecsv(nstsv,nstsv))
allocate(h(nstsv,nstsv),v(nstsv,nstsv))
allocate(ylmgq(lmmaxo,ng2gk),sfacgq(ng2gk,natmtot))
allocate(wfmt1(npcmtmax,natmtot,nspinor,nstsv),wfir1(ngtot,nspinor,nstsv))
allocate(wfmt2(npcmtmax,natmtot,nspinor,nstsv),wfir2(ngtot,nspinor,nstsv))
allocate(zrhomt(npcmtmax,natmtot,nstsv),zrhoir(ngtot,nstsv))
! coefficient of long-range term
cfq=0.5d0*(omega/pi)**2
! get the eigenvectors from file for input reduced k-point
call getevecfv(filext,ikp,vkl(:,ikp),vgkl(:,:,:,ikp),evecfv)
! find the matching coefficients
call match(ngk(1,ikp),gkc(:,1,ikp),tpgkc(:,:,1,ikp),sfacgk(:,:,1,ikp),apwalm)
! index to all states
do ist1=1,nstsv
  idx(ist1)=ist1
end do
! calculate the wavefunctions for all states of the input k-point
call genwfsv(.false.,.true.,nstsv,idx,ngk(1,ikp),igkig(:,1,ikp),apwalm, &
 evecfv,evecsvp,wfmt1,ngtot,wfir1)
!-----------------------------------------!
!     local potential matrix elements     !
!-----------------------------------------!
if (hybrid.and.spinpol) then
! magnetic field matrix elements in hybrid case
  call genvbmatk(vmt,vir,bmt,bir,ngk(1,ikp),igkig(:,1,ikp),wfmt1,ngtot,wfir1,h)
else
  call genvmatk(vmt,vir,ngk(1,ikp),igkig(:,1,ikp),wfmt1,ngtot,wfir1,h)
end if
! Fourier transform wavefunctions to real-space
call zftwfir(ngk(1,ikp),igkig(:,1,ikp),wfir1)
!---------------------------------!
!     kinetic matrix elements     !
!---------------------------------!
allocate(kmat(nstsv,nstsv))
call getkmat(ikp,kmat)
call zgemm('N','N',nstsv,nstsv,nstsv,zone,kmat,nstsv,evecsvp,nstsv,zzero,v, &
 nstsv)
call zgemm('C','N',nstsv,nstsv,nstsv,zone,evecsvp,nstsv,v,nstsv,zone,h,nstsv)
deallocate(kmat)
!------------------------------!
!     Fock matrix elements     !
!------------------------------!
v(:,:)=0.d0
! loop over non-reduced k-point set
do ik=1,nkptnr
! find the equivalent reduced k-point
  jk=ivkik(ivk(1,ik),ivk(2,ik),ivk(3,ik))
! determine the q-vector
  iv(:)=ivk(:,ikp)-ivk(:,ik)
  iv(:)=modulo(iv(:),ngridk(:))
  iq=iqmap(iv(1),iv(2),iv(3))
  vc(:)=vkc(:,ikp)-vkc(:,ik)
  do ig=1,ng2gk
! determine G+q-vectors
    vgqc(:,ig)=vgc(:,ig)+vc(:)
! G+q-vector length and (theta, phi) coordinates
    call sphcrd(vgqc(:,ig),gqc(ig),tp)
! spherical harmonics for G+q-vector
    call genylm(lmaxo,tp,ylmgq(:,ig))
  end do
! structure factors for G+q
  call gensfacgp(ng2gk,vgqc,ng2gk,sfacgq)
! find the shortest G+q-vector
  call findigp0(ng2gk,gqc,igq0)
  sfacgq0(:)=sfacgq(igq0,:)
! compute the required spherical Bessel functions
  call genjlgprmt(lnpsd,ng2gk,gqc,ng2gk,jlgqrmt)
  call genjlgq0r(gqc(igq0),jlgq0r)
! find the matching coefficients
  call match(ngk(1,ik),gkc(:,1,ik),tpgkc(:,:,1,ik),sfacgk(:,:,1,ik),apwalm)
! get the eigenvectors from file for non-reduced k-point
  call getevecfv(filext,0,vkl(:,ik),vgkl(:,:,1,ik),evecfv)
  call getevecsv(filext,0,vkl(:,ik),evecsv)
! calculate the wavefunctions for all states
  call genwfsv(.false.,.false.,nstsv,idx,ngk(1,ik),igkig(:,1,ik),apwalm, &
   evecfv,evecsv,wfmt2,ngtot,wfir2)
  do ist3=1,nstsv
    if (abs(occsv(ist3,jk)).lt.epsocc) cycle
! calculate the complex overlap densities for all states (T. McQueen)
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
    do ist1=1,nstsv
      call genzrho(.true.,.true.,wfmt2(:,:,:,ist3),wfir2(:,:,ist3), &
       wfmt1(:,:,:,ist1),wfir1(:,:,ist1),zrhomt(:,:,ist1),zrhoir(:,ist1))
    end do
!$OMP END DO
!$OMP END PARALLEL
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zvclmt,zvclir,zgq01,zgq02) &
!$OMP PRIVATE(ist1,z1,z2,t1)
!$OMP DO
    do ist2=1,nstsv
      allocate(zvclmt(npcmtmax,natmtot),zvclir(ngtot))
! calculate the Coulomb potential
      call genzvclmt(nrcmt,nrcmti,nrcmtmax,rcmt,npcmtmax,zrhomt(:,:,ist2), &
       zvclmt)
      call zpotcoul(nrcmt,nrcmti,npcmt,npcmti,nrcmtmax,rcmt,ng2gk,igq0,gqc, &
       ng2gk,jlgqrmt,ylmgq,sfacgq,zrhoir(:,ist2),npcmtmax,zvclmt,zvclir,zgq02)
      do ist1=1,ist2
        z1=zfinp(zrhomt(:,:,ist1),zrhoir(:,ist1),zvclmt,zvclir)
! compute the density coefficient of the smallest G+q-vector
        call zftgq0(jlgq0r,ylmgq(:,igq0),sfacgq0,zrhomt(:,:,ist1), &
         zrhoir(:,ist1),zgq01)
        z2=cfq*wiq2(iq)*conjg(zgq01)*zgq02
        t1=occsv(ist3,jk)/occmax
        v(ist1,ist2)=v(ist1,ist2)-t1*(wkptnr*z1+z2)
      end do
      deallocate(zvclmt,zvclir)
    end do
!$OMP END DO
!$OMP END PARALLEL
  end do
! end loop over non-reduced k-point set
end do
deallocate(vgqc,gqc,jlgqrmt,jlgq0r)
deallocate(ylmgq,sfacgq,apwalm,evecfv)
deallocate(wfmt1,wfmt2,wfir1,wfir2)
deallocate(zrhomt,zrhoir)
! scale the Coulomb matrix elements in the case of a hybrid functional
if (hybrid) v(:,:)=hybridc*v(:,:)
! add the Coulomb matrix elements to Hamiltonian
h(:,:)=h(:,:)+v(:,:)
!----------------------------------------------!
!     diagonalise Hartree-Fock Hamiltonian     !
!----------------------------------------------!
call eveqnz(nstsv,nstsv,h,evalsv(:,ikp))
! apply unitary transformation to second-variational states
evecsv(:,:)=evecsvp(:,:)
call zgemm('N','N',nstsv,nstsv,nstsv,zone,evecsv,nstsv,h,nstsv,zzero,evecsvp, &
 nstsv)
deallocate(evecsv,h,v)
return
end subroutine

