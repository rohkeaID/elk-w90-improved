
! Copyright (C) 2017-18 Arsenii Gerasimov, Yaroslav Kvashnin and Lars Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modw90overlap

contains

!BOP
! !ROUTINE: genw90overlap
! !INTERFACE:
subroutine genw90overlap(wfmt,wfir,nproj,wfmtq,wfirq,omn,phase)
! !USES:
use modmain
use modrandom
use modw90
! !INPUT/OUTPUT PARAMETERS:
!   wfmt  : muffin-tin wavefunction for all states at k-points in
!           spherical coordinates (in,complex(npcmtmax,natmtot,nspinor,wann_nband))
!   wfir  : interstitial wavefunction for all states at k-points in
!           real-space (in,complex(ngtot,nspinor,wann_nband))
!   nproj : number of bands for Wannierization or number of projections per spin
!           (in,integer)
!   wfmtq : muffin-tin wavefunction for all states at (k+b)-points in
!           spherical coordinates (in,complex(npcmtmax,natmtot,nspinor,nproj))
!   wfirq : interstitial wavefunction for all states at (k+b)-points in
!           real-space (in,complex(ngtot,nspinor,nproj))
!   omn   : general overlap matrices
!           (inout,complex(wann_nband,nspinor,nproj,nspinor))
!   phase : phase factor exp(ikr) in the muffin-tins
!           (in,complex(npcmtmax,natmtot),optional)
! !DESCRIPTION:
!   Calculates both $M_{mn}$ and $A_{mn}$ overlap integrals required for
!   the Wannier90.
!
!   The overlaps between cell periodic part of the Bloch states are given by
!   $$M_{mn}^{\bf k,\bf b}=\langle u_{m{\bf k}}|u_{n{\bf k}+{\bf b}}\rangle,$$
!   where
!   $|u_{n{\bf k}}\rangle$ denotes the Bloch states of $n$-th band
!   at each $\bf k$-point.
!
!   A starting guess projection of the Bloch states onto trial localised
!   orbitals are given by
!   $$A_{mn}^{(\bf k )}=\langle \psi_{m{\bf k}}|g_{n}\rangle,$$
!   where
!   $|\psi_{n\bf{k}}\rangle$ denotes the Bloch states and $|g_{n}\rangle$ -
!   trial localised orbitals.
!
! !REVISION HISTORY:
!   Created October 2017 (Arsenii Gerasimov)
!EOP
!BOC
implicit none
! input/output
complex(8),           intent(in   ) :: wfmt(npcmtmax,natmtot,nspinor,wann_nband)
complex(8),           intent(in   ) :: wfir(ngtot,nspinor,wann_nband)
integer,              intent(in   ) :: nproj
complex(8),           intent(in   ) :: wfmtq(npcmtmax,natmtot,nspinor,nproj)
complex(8),           intent(in   ) :: wfirq(ngtot,nspinor,nproj)
complex(8),           intent(inout) :: omn(wann_nband,nspinor,nproj,nspinor)
complex(8), optional, intent(in   ) :: phase(npcmtmax,natmtot)
! local variables
integer    :: ist,jst,ispin,jspin,i,irc,is,ia,ias,nrc,nrci,npc
real(8)    :: t1real,t1img
complex(8) :: t1
! automatic arrays
complex(8) :: fr(nrcmtmax)
! external functions
real(8)       fintgt
complex(8)    zdotu
external      fintgt, zdotu
! allocatable arrays
complex(8), allocatable :: wfmtq1(:,:,:,:)
complex(8), allocatable :: zrhomt(:,:),zrhoir(:)
complex(8), allocatable :: omnmt(:,:,:,:),omnir(:,:,:,:)
complex(8), allocatable :: cfunir_complex(:)
!-------------------------------------------------------------------------------

! Generate the wavefunctions at muffin-tins for all states at k+q
allocate(wfmtq1(npcmtmax,natmtot,nspinor,nproj))
if(present(phase))then ! Mmn case
  do ist = 1, nproj
    do ispin = 1,nspinor
      do is = 1, nspecies
        npc=npcmt(is)
        do ia = 1,natoms(is)
          ias = idxas(ia,is)
          wfmtq1(1:npc,ias,ispin,ist) = conjg(phase(1:npc,ias)) &
                                              *wfmtq(1:npc,ias,ispin,ist)
        end do
      end do
    end do
  end do
else ! Amn case
  wfmtq1 = wfmtq
endif

! Calculate matrix elements
allocate(zrhomt(npcmtmax,natmtot),zrhoir(ngtot))
allocate(omnmt(wann_nband,nspinor,nproj,nspinor))
allocate(omnir(wann_nband,nspinor,nproj,nspinor))
allocate(cfunir_complex(ngtot))

omnmt = zzero; omnir = zzero
cfunir_complex = zone*cfunir

do jst = 1,nproj
  do ist = 1,wann_nband
    do jspin = 1,nspinor
      do ispin = 1,nspinor

        call genzrho(.true.,.false.,ngtot,wfmt(:,:,ispin,ist)  , &
                                          wfir(:,ispin,ist)    , &
                                          wfmtq1(:,:,jspin,jst), &
                                          wfirq(:,jspin,jst)   , &
                                          zrhomt,zrhoir)

        do ias=1,natmtot
          is=idxis(ias)
          nrc=nrcmt(is)
          nrci=nrcmti(is)
          i=1
          do irc=1,nrci
            fr(irc)=zrhomt(i,ias)*r2cmt(irc,is)
            i=i+lmmaxi
          enddo
          do irc=nrci+1,nrc
            fr(irc)=zrhomt(i,ias)*r2cmt(irc,is)
            i=i+lmmaxo
          enddo
          t1real = fintgt(-1,nrc,rcmt(:,is),dble(fr))
          t1img = fintgt(-1,nrc,rcmt(:,is),aimag(fr))
          t1 = cmplx(t1real,t1img,kind=8)
          omnmt(ist,ispin,jst,jspin) = omnmt(ist,ispin,jst,jspin) + &
                                                          fourpi*y00*t1
        enddo

        ! Find the interstitial Omn
        t1 = zdotu(ngtot,zrhoir,1,cfunir_complex,1)
        omnir(ist,ispin,jst,jspin) = t1 * omega / dble(ngtot)

      enddo
    enddo
  enddo
  ! Total calculated Omn
  omn(:,:,jst,:) = omnmt(:,:,jst,:) + omnir(:,:,jst,:)
enddo

deallocate(cfunir_complex,omnmt,omnir,zrhomt,zrhoir,wfmtq1)

end subroutine genw90overlap
!EOC

end module modw90overlap
