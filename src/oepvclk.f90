
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine oepvclk(ikp,vclcv,vclvv)
use modmain
implicit none
! arguments
integer, intent(in) :: ikp
complex(8), intent(out) :: vclcv(ncrmax,natmtot,nstsv)
complex(8), intent(out) :: vclvv(nstsv,nstsv)
! local variables
integer ik,jk,nst,ist1,ist2,ist3
integer is,ia,ias,nrc,nrci,npc
integer iv(3),ig,iq,igq0
integer ic,jc,m1,m2
real(8) v(3),tp(2),cfq
complex(8) zgq01,zgq02,z1,z2
! automatic arrays
integer idx(nstsv)
complex(8) sfacgq0(natmtot)
! allocatable arrays
real(8), allocatable :: vgqc(:,:),gqc(:)
real(8), allocatable :: jlgqrmt(:,:,:),jlgq0r(:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: evecfv(:,:),evecsv(:,:)
complex(8), allocatable :: ylmgq(:,:),sfacgq(:,:)
complex(8), allocatable :: wfmt1(:,:,:,:),wfmt2(:,:,:,:)
complex(8), allocatable :: wfir1(:,:,:),wfir2(:,:,:)
complex(8), allocatable :: wfcr1(:,:),wfcr2(:,:)
complex(8), allocatable :: zrhomt1(:,:,:),zrhomt2(:,:),zrhoir1(:,:)
complex(8), allocatable :: zvclmt(:,:),zvclir(:),zfmt(:)
! external functions
complex(8) zfinp,zfmtinp
external zfinp,zfmtinp
! allocate local arrays
allocate(vgqc(3,ng2gk),gqc(ng2gk))
allocate(jlgqrmt(0:lnpsd,ng2gk,nspecies),jlgq0r(0:lmaxo,nrcmtmax,nspecies))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv(nmatmax,nstfv),evecsv(nstsv,nstsv))
allocate(ylmgq(lmmaxo,ng2gk),sfacgq(ng2gk,natmtot))
allocate(wfmt1(npcmtmax,natmtot,nspinor,nstsv))
allocate(wfmt2(npcmtmax,natmtot,nspinor,nstsv))
allocate(wfir1(ngtot,nspinor,nstsv),wfir2(ngtot,nspinor,nstsv))
allocate(wfcr1(npcmtmax,2),wfcr2(npcmtmax,2))
allocate(zrhomt1(npcmtmax,natmtot,nstsv),zrhomt2(npcmtmax,nstcr))
allocate(zrhoir1(ngtot,nstsv))
allocate(zvclmt(npcmtmax,natmtot),zvclir(ngtot))
allocate(zfmt(npcmtmax))
! factor for long-range term
cfq=0.5d0*(omega/pi)**2
! zero the Coulomb matrix elements
vclcv(:,:,:)=0.d0
vclvv(:,:)=0.d0
! get the eigenvectors from file for input reduced k-point
call getevecfv(filext,ikp,vkl(:,ikp),vgkl(:,:,:,ikp),evecfv)
call getevecsv(filext,ikp,vkl(:,ikp),evecsv)
! find the matching coefficients
call match(ngk(1,ikp),gkc(:,1,ikp),tpgkc(:,:,1,ikp),sfacgk(:,:,1,ikp),apwalm)
! index to all states
do ist1=1,nstsv
  idx(ist1)=ist1
end do
! calculate the wavefunctions for all states of the input k-point
call genwfsv(.false.,.false.,nstsv,idx,ngk(1,ikp),igkig(:,1,ikp),apwalm, &
 evecfv,evecsv,wfmt1,ngtot,wfir1)
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
! determine G+q-vectors
    vgqc(:,ig)=vgc(:,ig)+v(:)
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
! get the eigenvectors from file for non-reduced k-points
  call getevecfv(filext,0,vkl(:,ik),vgkl(:,:,1,ik),evecfv)
  call getevecsv(filext,0,vkl(:,ik),evecsv)
! count and index occupied states
  nst=0
  do ist3=1,nstsv
    if (evalsv(ist3,jk).lt.efermi) then
      nst=nst+1
      idx(nst)=ist3
    end if
  end do
! calculate the wavefunctions for occupied states
  call genwfsv(.false.,.false.,nst,idx,ngk(1,ik),igkig(:,1,ik),apwalm,evecfv, &
   evecsv,wfmt2,ngtot,wfir2)
  do ist3=1,nst
! compute the complex overlap densities for all valence-valence states
    do ist1=1,nstsv
      call genzrho(.true.,.true.,wfmt2(:,:,:,ist3),wfir2(:,:,ist3), &
       wfmt1(:,:,:,ist1),wfir1(:,:,ist1),zrhomt1(:,:,ist1),zrhoir1(:,ist1))
    end do
! compute the complex overlap densities for all valence-core states
    jc=0
    do is=1,nspecies
      nrc=nrcmt(is)
      nrci=nrcmti(is)
      npc=npcmt(is)
      do ia=1,natoms(is)
        ias=idxas(ia,is)
        do ist1=1,nstsp(is)
          if (spcore(ist1,is)) then
            do m1=-ksp(ist1,is),ksp(ist1,is)-1
              jc=jc+1
! generate the core wavefunction in spherical coordinates (pass in m-1/2)
              call wavefcr(.false.,lradstp,is,ia,ist1,m1,npcmtmax,wfcr1)
              if (spinpol) then
                call zrho2(npc,wfmt2(:,ias,1,ist3),wfmt2(:,ias,2,ist3), &
                 wfcr1(:,1),wfcr1(:,2),zfmt)
              else
                call zrho1(npc,wfmt2(:,ias,1,ist3),wfcr1(:,1),zfmt)
              end if
! convert to spherical harmonics
              call zfsht(nrc,nrci,zfmt,zrhomt2(:,jc))
            end do
          end if
        end do
      end do
    end do
    do ist2=1,nstsv
      if (evalsv(ist2,ikp).gt.efermi) then
! calculate the Coulomb potential
        call genzvclmt(nrcmt,nrcmti,nrcmtmax,rcmt,npcmtmax,zrhomt1(:,:,ist2), &
         zvclmt)
        call zpotcoul(nrcmt,nrcmti,npcmt,npcmti,nrcmtmax,rcmt,ng2gk,igq0,gqc, &
         ng2gk,jlgqrmt,ylmgq,sfacgq,zrhoir1(:,ist2),npcmtmax,zvclmt,zvclir, &
         zgq02)
!----------------------------------------------!
!     valence-valence-valence contribution     !
!----------------------------------------------!
        do ist1=1,nstsv
          if (evalsv(ist1,ikp).lt.efermi) then
            z1=zfinp(zrhomt1(:,:,ist1),zrhoir1(:,ist1),zvclmt,zvclir)
! compute the density coefficient of the smallest G+q-vector
            call zftgq0(jlgq0r,ylmgq(:,igq0),sfacgq0,zrhomt1(:,:,ist1), &
             zrhoir1(:,ist1),zgq01)
            z2=cfq*wiq2(iq)*(conjg(zgq01)*zgq02)
            vclvv(ist1,ist2)=vclvv(ist1,ist2)-(wkptnr*z1+z2)
          end if
        end do
!-------------------------------------------!
!     core-valence-valence contribution     !
!-------------------------------------------!
        jc=0
        do is=1,nspecies
          nrc=nrcmt(is)
          nrci=nrcmti(is)
          do ia=1,natoms(is)
            ias=idxas(ia,is)
            ic=0
            do ist1=1,nstsp(is)
              if (spcore(ist1,is)) then
                do m1=-ksp(ist1,is),ksp(ist1,is)-1
                  ic=ic+1
                  jc=jc+1
                  z1=zfmtinp(nrc,nrci,rcmt(:,is),r2cmt(:,is),zrhomt2(:,jc), &
                   zvclmt(:,ias))
                  vclcv(ic,ias,ist2)=vclcv(ic,ias,ist2)-wkptnr*z1
                end do
! end loop over ist1
              end if
            end do
! end loops over atoms and species
          end do
        end do
! end loop over ist2
      end if
    end do
! end loop over ist3
  end do
! end loop over non-reduced k-point set
end do
! begin loops over atoms and species
do is=1,nspecies
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  npc=npcmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ist3=1,nstsp(is)
      if (spcore(ist3,is)) then
        do m1=-ksp(ist3,is),ksp(ist3,is)-1
! generate the core wavefunction in spherical coordinates (pass in m-1/2)
          call wavefcr(.false.,lradstp,is,ia,ist3,m1,npcmtmax,wfcr1)
! compute the complex overlap densities for the core-valence states
          do ist1=1,nstsv
            if (spinpol) then
              call zrho2(npc,wfcr1(:,1),wfcr1(:,2),wfmt1(:,ias,1,ist1), &
               wfmt1(:,ias,2,ist1),zfmt)
            else
              call zrho1(npc,wfcr1(:,1),wfmt1(:,ias,1,ist1),zfmt)
            end if
            call zfsht(nrc,nrci,zfmt,zrhomt1(:,ias,ist1))
          end do
! compute the complex overlap densities for the core-core states
          ic=0
          do ist1=1,nstsp(is)
            if (spcore(ist1,is)) then
              do m2=-ksp(ist1,is),ksp(ist1,is)-1
                ic=ic+1
                call wavefcr(.false.,lradstp,is,ia,ist1,m2,npcmtmax,wfcr2)
                call zrho2(npc,wfcr1(:,1),wfcr1(:,2),wfcr2(:,1),wfcr2(:,2),zfmt)
                call zfsht(nrc,nrci,zfmt,zrhomt2(:,ic))
              end do
            end if
          end do
          do ist2=1,nstsv
            if (evalsv(ist2,ikp).gt.efermi) then
! calculate the Coulomb potential
              call zpotclmt(nrc,nrci,rcmt(:,is),zrhomt1(:,ias,ist2),zvclmt)
!-------------------------------------------!
!     valence-core-valence contribution     !
!-------------------------------------------!
              do ist1=1,nstsv
                if (evalsv(ist1,ikp).lt.efermi) then
                  z1=zfmtinp(nrc,nrci,rcmt(:,is),r2cmt(:,is), &
                   zrhomt1(:,ias,ist1),zvclmt)
                  vclvv(ist1,ist2)=vclvv(ist1,ist2)-z1
                end if
              end do
!----------------------------------------!
!     core-core-valence contribution     !
!----------------------------------------!
              ic=0
              do ist1=1,nstsp(is)
                if (spcore(ist1,is)) then
                  do m2=-ksp(ist1,is),ksp(ist1,is)-1
                    ic=ic+1
                    z1=zfmtinp(nrc,nrci,rcmt(:,is),r2cmt(:,is),zrhomt2(:,ic), &
                     zvclmt)
                    vclcv(ic,ias,ist2)=vclcv(ic,ias,ist2)-z1
                  end do
! end loop over ist1
                end if
              end do
! end loop over ist2
            end if
          end do
! end loops over ist3 and m1
        end do
      end if
    end do
! end loops over atoms and species
  end do
end do
deallocate(vgqc,gqc,jlgqrmt,jlgq0r)
deallocate(apwalm,evecfv,evecsv,ylmgq,sfacgq)
deallocate(wfmt1,wfmt2,wfir1,wfir2,wfcr1,wfcr2)
deallocate(zrhomt1,zrhomt2,zrhoir1)
deallocate(zvclmt,zvclir,zfmt)
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
!EOC

