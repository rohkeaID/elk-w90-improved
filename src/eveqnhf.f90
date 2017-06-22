
! Copyright (C) 2006 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine eveqnhf(ikp,vmt,vir,bmt,bir,evecsvp)
use modmain
implicit none
! arguments
integer, intent(in) :: ikp
real(8), intent(in) :: vmt(lmmaxvr,nrcmtmax,natmtot),vir(ngtot)
real(8), intent(in) :: bmt(lmmaxvr,nrcmtmax,natmtot,ndmag),bir(ngtot,ndmag)
complex(8), intent(inout) :: evecsvp(nstsv,nstsv)
! local variables
integer ik,jk,ist1,ist2,ist3,ispn
integer iv(3),igk,ifg,ig,iq,igq0
integer lwork,info
real(8) cfq,v(3),t1
complex(8) zrho01,zrho02,z1,z2
! automatic arrays
integer idx(nstsv)
complex(8) sfacgq0(natmtot)
! allocatable arrays
real(8), allocatable :: vgqc(:,:),gqc(:),tpgqc(:,:)
real(8), allocatable :: jlgqr(:,:,:),jlgq0r(:,:,:)
real(8), allocatable :: rwork(:)
complex(8), allocatable :: apwalm(:,:,:,:),evecfv(:,:),evecsv(:,:)
complex(8), allocatable :: h(:,:),c(:,:),kmat(:,:),z(:),work(:)
complex(8), allocatable :: ylmgq(:,:),sfacgq(:,:)
complex(8), allocatable :: wfmt1(:,:,:,:,:),wfmt2(:,:,:,:,:)
complex(8), allocatable :: wfir1(:,:,:),wfir2(:,:,:)
complex(8), allocatable :: zrhomt(:,:,:,:),zrhoir(:,:)
complex(8), allocatable :: zvclmt(:,:,:),zvclir(:)
! external functions
complex(8) zfinp
external zfinp
!$OMP CRITICAL
write(*,'("Info(eveqnhf): ",I6," of ",I6," k-points")') ikp,nkpt
!$OMP END CRITICAL
! allocate local arrays
allocate(vgqc(3,ngvec),gqc(ngvec),tpgqc(2,ngvec))
allocate(jlgqr(0:lnpsd,ngvec,nspecies),jlgq0r(0:lmaxvr,nrcmtmax,nspecies))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv(nmatmax,nstfv),evecsv(nstsv,nstsv))
allocate(h(nstsv,nstsv),c(nstsv,nstsv))
allocate(ylmgq(lmmaxvr,ngvec),sfacgq(ngvec,natmtot))
allocate(wfmt1(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfmt2(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfir1(ngtot,nspinor,nstsv),wfir2(ngtot,nspinor,nstsv))
allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot,nstsv),zrhoir(ngtot,nstsv))
! coefficient of long-range term
cfq=0.5d0*(omega/pi)**2
! get the eigenvectors from file for input k-point
call getevecfv(filext,vkl(:,ikp),vgkl(:,:,:,ikp),evecfv)
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
allocate(z(ngkmax))
t1=1.d0/sqrt(omega)
do ist1=1,nstsv
  do ispn=1,nspinor
    call zcopy(ngk(1,ikp),wfir1(:,ispn,ist1),1,z,1)
    wfir1(:,ispn,ist1)=0.d0
    do igk=1,ngk(1,ikp)
      ifg=igfft(igkig(igk,1,ikp))
      wfir1(ifg,ispn,ist1)=t1*z(igk)
    end do
    call zfftifc(3,ngridg,1,wfir1(:,ispn,ist1))
  end do
end do
deallocate(z)
!---------------------------------!
!     kinetic matrix elements     !
!---------------------------------!
allocate(kmat(nstsv,nstsv))
call getkmat(ikp,kmat)
call zgemm('N','N',nstsv,nstsv,nstsv,zone,kmat,nstsv,evecsvp,nstsv,zzero,c, &
 nstsv)
call zgemm('C','N',nstsv,nstsv,nstsv,zone,evecsvp,nstsv,c,nstsv,zone,h,nstsv)
deallocate(kmat)
!---------------------------------------------------------!
!     Coulomb valence-valence-valence matrix elements     !
!---------------------------------------------------------!
c(:,:)=0.d0
! start loop over non-reduced k-point set
do ik=1,nkptnr
! find the equivalent reduced k-point
  jk=ivkik(ivk(1,ik),ivk(2,ik),ivk(3,ik))
! determine q-vector
  iv(:)=ivk(:,ikp)-ivk(:,ik)
  iv(:)=modulo(iv(:),ngridk(:))
  iq=iqmap(iv(1),iv(2),iv(3))
  v(:)=vkc(:,ikp)-vkc(:,ik)
  do ig=1,ngvec
! determine G+q-vectors
    vgqc(:,ig)=vgc(:,ig)+v(:)
! G+q-vector length and (theta, phi) coordinates
    call sphcrd(vgqc(:,ig),gqc(ig),tpgqc(:,ig))
! spherical harmonics for G+q-vector
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
! find the matching coefficients
  call match(ngk(1,ik),gkc(:,1,ik),tpgkc(:,:,1,ik),sfacgk(:,:,1,ik),apwalm)
! get the eigenvectors from file for non-reduced k-point
  call getevecfv(filext,vkl(:,ik),vgkl(:,:,1,ik),evecfv)
  call getevecsv(filext,vkl(:,ik),evecsv)
! calculate the wavefunctions for all states
  call genwfsv(.false.,.false.,nstsv,idx,ngk(1,ik),igkig(:,1,ik),apwalm, &
   evecfv,evecsv,wfmt2,ngtot,wfir2)
  do ist3=1,nstsv
    if (abs(occsv(ist3,jk)).lt.epsocc) cycle
! calculate the complex overlap densities for all states (T. McQueen)
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
    do ist1=1,nstsv
      call genzrho(.true.,.true.,wfmt2(:,:,:,:,ist3),wfir2(:,:,ist3), &
       wfmt1(:,:,:,:,ist1),wfir1(:,:,ist1),zrhomt(:,:,:,ist1),zrhoir(:,ist1))
    end do
!$OMP END DO
!$OMP END PARALLEL
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zvclmt,zvclir,zrho01,zrho02) &
!$OMP PRIVATE(ist1,z1,z2,t1)
!$OMP DO
    do ist2=1,nstsv
      allocate(zvclmt(lmmaxvr,nrcmtmax,natmtot),zvclir(ngtot))
! calculate the Coulomb potential
      call genzvclmt(nrcmt,nrcmtinr,nrcmtmax,rcmt,nrcmtmax,zrhomt(:,:,:,ist2), &
       zvclmt)
      call zpotcoul(nrcmt,nrcmtinr,nrcmtmax,rcmt,igq0,gqc,jlgqr,ylmgq,sfacgq, &
       zrhoir(:,ist2),nrcmtmax,zvclmt,zvclir,zrho02)
      do ist1=1,ist2
        z1=zfinp(zrhomt(:,:,:,ist1),zrhoir(:,ist1),zvclmt,zvclir)
! compute the density coefficient of the smallest G+q-vector
        call zrhogp(jlgq0r,ylmgq(:,igq0),sfacgq0,zrhomt(:,:,:,ist1), &
         zrhoir(:,ist1),zrho01)
        z2=cfq*wiq2(iq)*conjg(zrho01)*zrho02
        t1=occsv(ist3,jk)/occmax
        c(ist1,ist2)=c(ist1,ist2)-t1*(wkptnr*z1+z2)
      end do
      deallocate(zvclmt,zvclir)
    end do
!$OMP END DO
!$OMP END PARALLEL
  end do
! end loop over non-reduced k-point set
end do
deallocate(vgqc,tpgqc,gqc,jlgqr,jlgq0r)
deallocate(ylmgq,sfacgq,apwalm,evecfv)
deallocate(wfmt1,wfmt2,wfir1,wfir2)
deallocate(zrhomt,zrhoir)
! scale the Coulomb matrix elements in the case of a hybrid functional
if (hybrid) c(:,:)=hybridc*c(:,:)
! add the Coulomb matrix elements to Hamiltonian
h(:,:)=h(:,:)+c(:,:)
!----------------------------------------------!
!     diagonalise Hartree-Fock Hamiltonian     !
!----------------------------------------------!
allocate(rwork(3*nstsv))
lwork=2*nstsv
allocate(work(lwork))
call zheev('V','U',nstsv,h,nstsv,evalsv(:,ikp),work,lwork,rwork,info)
deallocate(rwork,work)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(eveqnhf): diagonalisation of the Hartree-Fock Hamiltonian &
   &failed")')
  write(*,'(" for k-point ",I8)') ikp
  write(*,'(" ZHEEV returned INFO = ",I8)') info
  write(*,*)
  stop
end if
! apply unitary transformation to second-variational states
evecsv(:,:)=evecsvp(:,:)
call zgemm('N','N',nstsv,nstsv,nstsv,zone,evecsv,nstsv,h,nstsv,zzero,evecsvp, &
 nstsv)
deallocate(evecsv,h,c)
return
end subroutine

