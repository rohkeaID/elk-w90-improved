
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine exxengyk(ikp)
use modmain
implicit none
! arguments
integer, intent(in) :: ikp
! local variables
integer jkp,ik,jk
integer nst1,nst2,ist,jst
integer is,ia,ias
integer nrc,nrci,m,npc
integer iv(3),ig,iq,igq0
real(8) ex,cfq,v(3),tp(2),t1
complex(8) zgq0,z1
! automatic arrays
integer idx(nstsv)
! allocatable arrays
real(8), allocatable :: vgqc(:,:),gqc(:),jlgqrmt(:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: evecfv(:,:),evecsv(:,:)
complex(8), allocatable :: ylmgq(:,:),sfacgq(:,:)
complex(8), allocatable :: wfmt1(:,:,:,:),wfir1(:,:,:)
complex(8), allocatable :: wfmt2(:,:,:,:),wfir2(:,:,:)
complex(8), allocatable :: wfcr(:,:),zfmt(:)
complex(8), allocatable :: zrhomt(:,:),zrhoir(:)
complex(8), allocatable :: zvclmt(:,:),zvclir(:)
! external functions
complex(8) zfinp,zfmtinp
external zfinp,zfmtinp
! allocate local arrays
allocate(vgqc(3,ng2gk),gqc(ng2gk))
allocate(jlgqrmt(0:lnpsd,ng2gk,nspecies))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv(nmatmax,nstfv),evecsv(nstsv,nstsv))
allocate(ylmgq(lmmaxo,ng2gk),sfacgq(ng2gk,natmtot))
allocate(wfmt2(npcmtmax,natmtot,nspinor,nstsv))
allocate(wfir2(ngtot,nspinor,nstsv))
allocate(zrhomt(npcmtmax,natmtot),zrhoir(ngtot))
allocate(zvclmt(npcmtmax,natmtot),zvclir(ngtot))
! coefficient for long-range term
cfq=0.5d0*(omega/pi)**2
! get the eigenvectors from file for input reduced k-point
call getevecfv(filext,ikp,vkl(:,ikp),vgkl(:,:,:,ikp),evecfv)
call getevecsv(filext,ikp,vkl(:,ikp),evecsv)
! find the matching coefficients
call match(ngk(1,ikp),gkc(:,1,ikp),tpgkc(:,:,1,ikp),sfacgk(:,:,1,ikp),apwalm)
! equivalent reduced input k-point
jkp=ivkik(ivk(1,ikp),ivk(2,ikp),ivk(3,ikp))
! count and index the occupied states
nst1=0
do ist=1,nstsv
  if (evalsv(ist,jkp).lt.efermi) then
    nst1=nst1+1
    idx(nst1)=ist
  end if
end do
! calculate the wavefunctions for occupied states of the input k-point
allocate(wfmt1(npcmtmax,natmtot,nspinor,nst1))
allocate(wfir1(ngtot,nspinor,nst1))
call genwfsv(.false.,.false.,nst1,idx,ngk(1,ikp),igkig(:,1,ikp),apwalm,evecfv, &
 evecsv,wfmt1,ngtot,wfir1)
! zero the local exchange energy variable
ex=0.d0
! start loop over non-reduced k-point set
do ik=1,nkptnr
! equivalent reduced k-point
  jk=ivkik(ivk(1,ik),ivk(2,ik),ivk(3,ik))
! determine q-vector
  iv(:)=ivk(:,ikp)-ivk(:,ik)
  iv(:)=modulo(iv(:),ngridk(:))
  iq=iqmap(iv(1),iv(2),iv(3))
  v(:)=vkc(:,ikp)-vkc(:,ik)
  do ig=1,ng2gk
! determine the G+q-vectors
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
! find the matching coefficients
  call match(ngk(1,ik),gkc(:,1,ik),tpgkc(:,:,1,ik),sfacgk(:,:,1,ik),apwalm)
! get the eigenvectors from file for non-reduced k-points
  call getevecfv(filext,0,vkl(:,ik),vgkl(:,:,1,ik),evecfv)
  call getevecsv(filext,0,vkl(:,ik),evecsv)
! count and index the occupied states
  nst2=0
  do jst=1,nstsv
    if (evalsv(jst,jk).lt.efermi) then
      nst2=nst2+1
      idx(nst2)=jst
    end if
  end do
! calculate the wavefunctions for occupied states
  call genwfsv(.false.,.false.,nst2,idx,ngk(1,ik),igkig(:,1,ik),apwalm,evecfv, &
   evecsv,wfmt2,ngtot,wfir2)
!--------------------------------------------!
!    valence-valence-valence contribution    !
!--------------------------------------------!
  do jst=1,nst2
    do ist=1,nst1
! calculate the complex overlap density
      call genzrho(.true.,.true.,wfmt2(:,:,:,jst),wfir2(:,:,jst), &
       wfmt1(:,:,:,ist),wfir1(:,:,ist),zrhomt,zrhoir)
! calculate the Coulomb potential
      call genzvclmt(nrcmt,nrcmti,nrcmtmax,rcmt,npcmtmax,zrhomt,zvclmt)
      call zpotcoul(nrcmt,nrcmti,npcmt,npcmti,nrcmtmax,rcmt,ng2gk,igq0,gqc, &
       ng2gk,jlgqrmt,ylmgq,sfacgq,zrhoir,npcmtmax,zvclmt,zvclir,zgq0)
      z1=zfinp(zrhomt,zrhoir,zvclmt,zvclir)
      t1=cfq*wiq2(iq)*(dble(zgq0)**2+aimag(zgq0)**2)
      ex=ex-0.5d0*occmax*wkpt(ikp)*(wkptnr*dble(z1)+t1)
    end do
  end do
! end loop over non-reduced k-point set
end do
deallocate(vgqc,gqc,jlgqrmt)
deallocate(evecfv,evecsv)
deallocate(apwalm,ylmgq,sfacgq)
deallocate(wfmt2,wfir2)
!-----------------------------------------!
!    valence-core-valence contribution    !
!-----------------------------------------!
allocate(wfcr(npcmtmax,2),zfmt(npcmtmax))
! begin loops over atoms and species
do is=1,nspecies
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  npc=npcmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do jst=1,nstsp(is)
      if (spcore(jst,is)) then
        do m=-ksp(jst,is),ksp(jst,is)-1
! generate the core wavefunction in spherical coordinates (pass in m-1/2)
          call wavefcr(.false.,lradstp,is,ia,jst,m,npcmtmax,wfcr)
          do ist=1,nst1
! calculate the complex overlap density in spherical harmonics
            if (spinpol) then
              call zrho2(npc,wfcr(:,1),wfcr(:,2),wfmt1(:,ias,1,ist), &
               wfmt1(:,ias,2,ist),zfmt)
            else
              call zrho1(npc,wfcr(:,1),wfmt1(:,ias,1,ist),zfmt)
            end if
            call zfsht(nrc,nrci,zfmt,zrhomt(:,ias))
! calculate the Coulomb potential
            call zpotclmt(nrc,nrci,rcmt(:,is),zrhomt(:,ias),zvclmt(:,ias))
            z1=zfmtinp(nrc,nrci,rcmt(:,is),r2cmt(:,is),zrhomt(:,ias), &
             zvclmt(:,ias))
            ex=ex-occmax*wkpt(ikp)*dble(z1)
          end do
! end loop over m
        end do
! end loop over jst
      end if
    end do
! end loops over atoms and species
  end do
end do
! add to global exchange energy
!$OMP CRITICAL
engyx=engyx+ex
!$OMP END CRITICAL
deallocate(wfmt1,wfir1,wfcr,zfmt)
deallocate(zrhomt,zrhoir,zvclmt,zvclir)
return

contains

subroutine zrho1(n,x,y,z)
implicit none
integer, intent(in) :: n
complex(8), intent(in) :: x(n),y(n)
complex(8), intent(out) :: z(n)
z(:)=conjg(x(:))*y(:)
return
end subroutine

subroutine zrho2(n,x1,x2,y1,y2,z)
implicit none
integer, intent(in) :: n
complex(8), intent(in) :: x1(n),x2(n),y1(n),y2(n)
complex(8), intent(out) :: z(n)
z(:)=conjg(x1(:))*y1(:)+conjg(x2(:))*y2(:)
return
end subroutine

end subroutine


