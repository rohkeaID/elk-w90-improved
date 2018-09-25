! Copyright (C) 2017 Arsenii Gerasimov, Yaroslav Kvashnin and Lars Nordstrom
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !MODULE: w90overlap - had to make a module, as we have an optional var. in a routine
! !ROUTINE: genw90overlap
! !INTERFACE:
module modw90overlap
  contains
    subroutine genw90overlap(wfmt,wfir,nproj,ns,wfmtq,wfirq,omn,phase)
      ! !USES:
      use modmain
      use modrandom
      use modw90
      ! !DESCRIPTION:
      !   Calculates the general Omn (Mmn and Amn) overlap matrices required for Wannier90.
      !
      ! !REVISION HISTORY:
      !   Created June 2017 (Arsenii Gerasimov)
      !EOP
      !BOC
      implicit none
      ! input/output
      complex(8), dimension(npcmtmax,natmtot,nspinor,wann_nband), intent(in   ) :: wfmt
      complex(8), dimension(ngtot,nspinor,wann_nband),            intent(in   ) :: wfir
      integer,                                                    intent(in   ) :: nproj
      integer,                                                    intent(in   ) :: ns
      complex(8), dimension(npcmtmax,natmtot,nspinor,nproj),      intent(in   ) :: wfmtq
      complex(8), dimension(ngtot,nspinor,nproj),                 intent(in   ) :: wfirq
      complex(8), dimension(wann_nband,nspinor,nproj,ns),         intent(inout) :: omn
      complex(8), dimension(npcmtmax,natmtot),          optional, intent(in   ) :: phase

      ! local variables
      integer    :: ist,jst,ispin,jspin,i,irc
      complex(8) :: matel
      real(8)    :: norm
      integer    :: is,ia,ias
      integer    :: nrc,nrci,npc
      real(8)    :: t1real, t1img
      complex(8) :: t1
      ! automatic arrays
      complex(8) :: fr(nrcmtmax)
      ! external functions
      real(8)       fintgt
      complex(8)    zdotu
      external      fintgt, zdotu
      ! allocatable arrays
      complex(8), dimension(:,:,:,:), allocatable :: wfmtq1
      complex(8), dimension(:,:),     allocatable :: zrhomt
      complex(8), dimension(:),       allocatable :: zrhoir
      complex(8), dimension(:,:,:,:), allocatable :: omnmt
      complex(8), dimension(:,:,:,:), allocatable :: omnir
      complex(8), dimension(:),       allocatable :: cfunir_complex
!-------------------------------------------------------------------------------

      ! generate the wavefunctions for all states at k+q
      allocate(wfmtq1(npcmtmax,natmtot,nspinor,nproj))

      if(present(phase))then ! Mmn case
        do ist = 1, nproj
          do ispin = 1,nspinor
            do is = 1, nspecies
              npc=npcmt(is)
              do ia = 1,natoms(is)
                ias = idxas(ia,is)
                wfmtq1(1:npc,ias,ispin,ist) = conjg(phase(1:npc,ias))*wfmtq(1:npc,ias,ispin,ist)
              end do
            end do
          end do
        end do
      else ! Amn case
        wfmtq1 = wfmtq
      endif

      !Calculate matrix elements
      allocate(zrhomt(npcmtmax,natmtot),zrhoir(ngtot))
      allocate(omnmt(wann_nband,nspinor,nproj,ns))
      allocate(omnir(wann_nband,nspinor,nproj,ns))
      allocate(cfunir_complex(ngtot))

      omnmt = cmplx(0.0d0,0.0d0,kind=8)
      omnir = cmplx(0.0d0,0.0d0,kind=8)
      cfunir_complex = cmplx(cfunir,0.0d0,kind=8)

      do jst = 1,nproj
        if((.not. present(phase)) .and. wann_proj_isrand(jst)) then ! Amn case
          ! if this projection should be random (at present, any non-atom centred projection!)
          norm = 0.0d0
          do ispin = 1,nspinor
            do ist = 1,wann_nband
              matel = cmplx(randomu(),randomu(),kind=8)
              norm = norm + abs(matel)
              omn(ist,ispin,jst,1) = matel
            end do
          end do
          omn(:,:,jst,1) = omn(:,:,jst,1) / norm
        else ! Amn and Mmn case
          do ist = 1,wann_nband
            do jspin = 1,ns
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
                  omnmt(ist,ispin,jst,jspin) = omnmt(ist,ispin,jst,jspin) + fourpi*y00*t1
                enddo

                ! find the interstitial omn
                t1 = zdotu(ngtot,zrhoir,1,cfunir_complex,1)
                omnir(ist,ispin,jst,jspin) = t1 * omega / dble(ngtot)

              enddo
            enddo
          enddo
        endif
        ! total calculated charge
        omn(:,:,jst,:) = omnmt(:,:,jst,:) + omnir(:,:,jst,:)
      enddo


      deallocate(cfunir_complex)
      deallocate(omnmt,omnir)
      deallocate(zrhomt,zrhoir)
      deallocate(wfmtq1)


    end subroutine genw90overlap
    !EOC
end module modw90overlap
