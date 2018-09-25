
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine exxengyk(ikp)
use modmain
implicit none
! arguments
integer, intent(in) :: ikp
! local variables
integer iq,ik,jk,i,m
integer nst1,nst2,ist,jst
integer is,ia,ias
integer nrc,nrci,npc
integer iv(3),ig
real(8) ex,vc(3),tp(2)
complex(8) z1
! automatic arrays
integer idx(nstsv)
! allocatable arrays
real(8), allocatable :: vgqc(:,:),gqc(:),gclgq(:),jlgqrmt(:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:),evecfv(:,:),evecsv(:,:)
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
allocate(vgqc(3,ngvc),gqc(ngvc),gclgq(ngvc))
allocate(jlgqrmt(0:lnpsd,ngvc,nspecies))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv(nmatmax,nstfv),evecsv(nstsv,nstsv))
allocate(ylmgq(lmmaxo,ngvc),sfacgq(ngvc,natmtot))
allocate(wfmt2(npcmtmax,natmtot,nspinor,nstsv))
allocate(wfir2(ngtc,nspinor,nstsv))
allocate(zrhomt(npcmtmax,natmtot),zrhoir(ngtc))
allocate(zvclmt(npcmtmax,natmtot),zvclir(ngtc))
! get the eigenvectors from file for input reduced k-point
call getevecfv(filext,ikp,vkl(:,ikp),vgkl(:,:,:,ikp),evecfv)
call getevecsv(filext,ikp,vkl(:,ikp),evecsv)
! find the matching coefficients
call match(ngk(1,ikp),gkc(:,1,ikp),tpgkc(:,:,1,ikp),sfacgk(:,:,1,ikp),apwalm)
! count and index the occupied states
nst1=0
do ist=1,nstsv
  if (evalsv(ist,ikp).lt.efermi) then
    nst1=nst1+1
    idx(nst1)=ist
  end if
end do
! calculate the wavefunctions for occupied states of the input k-point
allocate(wfmt1(npcmtmax,natmtot,nspinor,nst1),wfir1(ngtc,nspinor,nst1))
call genwfsv(.false.,.false.,nst1,idx,ngdc,igfc,ngk(1,ikp),igkig(:,1,ikp), &
 apwalm,evecfv,evecsv,wfmt1,ngtc,wfir1)
! zero the local exchange energy variable
ex=0.d0
! start loop over non-reduced k-point set
do ik=1,nkptnr
! equivalent reduced k-point
  jk=ivkik(ivk(1,ik),ivk(2,ik),ivk(3,ik))
! determine the q-vector
  iv(:)=ivk(:,ikp)-ivk(:,ik)
  iv(:)=modulo(iv(:),ngridk(:))
! check if the q-point is in user-defined set
  iv(:)=iv(:)*ngridq(:)
  do i=1,3
    if (modulo(iv(i),ngridk(i)).ne.0) goto 10
  end do
  iv(:)=iv(:)/ngridk(:)
  iq=iqmap(iv(1),iv(2),iv(3))
  vc(:)=vkc(:,ikp)-vkc(:,ik)
  do ig=1,ngvc
! determine the G+q-vectors
    vgqc(:,ig)=vgc(:,ig)+vc(:)
! G+q-vector length and (theta, phi) coordinates
    call sphcrd(vgqc(:,ig),gqc(ig),tp)
! spherical harmonics for G+q-vectors
    call genylm(lmaxo,tp,ylmgq(:,ig))
  end do
! structure factors for G+q
  call gensfacgp(ngvc,vgqc,ngvc,sfacgq)
! generate the regularised Coulomb Green's function in G+q-space
  call gengclgq(.true.,iq,ngvc,gqc,gclgq)
! compute the required spherical Bessel functions
  call genjlgprmt(lnpsd,ngvc,gqc,ngvc,jlgqrmt)
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
  call genwfsv(.false.,.false.,nst2,idx,ngdc,igfc,ngk(1,ik),igkig(:,1,ik), &
   apwalm,evecfv,evecsv,wfmt2,ngtc,wfir2)
!--------------------------------------------!
!    valence-valence-valence contribution    !
!--------------------------------------------!
  do jst=1,nst2
    do ist=1,nst1
! calculate the complex overlap density
      call genzrho(.true.,.true.,ngtc,wfmt2(:,:,:,jst),wfir2(:,:,jst), &
       wfmt1(:,:,:,ist),wfir1(:,:,ist),zrhomt,zrhoir)
! calculate the Coulomb potential
      call genzvclmt(nrcmt,nrcmti,nrcmtmax,rcmt,npcmtmax,zrhomt,zvclmt)
      call zpotcoul(nrcmt,nrcmti,npcmt,npcmti,nrcmtmax,rcmt,ngdc,igfc,ngvc, &
       gqc,gclgq,ngvc,jlgqrmt,ylmgq,sfacgq,zrhoir,npcmtmax,zvclmt,zvclir)
      z1=zfinp(zrhomt,zrhoir,zvclmt,zvclir)
      ex=ex-0.5d0*occmax*wkpt(ikp)*wqptnr*dble(z1)
    end do
  end do
10 continue
! end loop over non-reduced k-point set
end do
deallocate(vgqc,gqc,gclgq,jlgqrmt)
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
!$OMP CRITICAL(exxengyk_)
engyx=engyx+ex
!$OMP END CRITICAL(exxengyk_)
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

