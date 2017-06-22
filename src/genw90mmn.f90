! Copyright (C) 2015 Jon Lafuente and Manh Duc Le
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genw90mmn
! !INTERFACE:
subroutine genw90mmn(ikp,vqpl,phase,gqc,ylmgq,sfacgq,matel_mmn,wfmt,wfir)
! !USES:
use modmain
use modw90
! !DESCRIPTION:
!   Calculates the Mmn overlap matrices required for Wannier90.
!
! !REVISION HISTORY:
!   Created March 2015 (Jon Lafuente and Manh Duc Le)
!EOP
!BOC
implicit none
! input/output
integer, intent(in) :: ikp
real(8), intent(in) :: vqpl(3)
complex(8), intent(in) :: phase(lmmaxvr,nrcmtmax,natmtot)
real(8), intent(in) :: gqc(ngrf)
complex(8), intent(in) :: ylmgq(lmmaxvr,ngrf)
complex(8), intent(in) :: sfacgq(ngrf,natmtot)
complex(8), intent(inout) :: matel_mmn(wann_nband,4,wann_nband,4)
complex(8), intent(in) :: wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,wann_nband)
complex(8), intent(in) :: wfir(ngtot,nspinor,wann_nband)
! local variables
integer ist,jst,ispin,jspin,i  
integer is,ia,ias,nrc,nrci
! automatic arrays
real(8) vkql(3)
! allocatable arrays
complex(8), allocatable :: wfmtq(:,:,:,:,:),wfirq(:,:,:), wfmtq1(:,:,:,:,:)
complex(8), allocatable :: zrhomt(:,:,:),zrhoir(:)
complex(8), allocatable :: zrhoig(:,:)

! k+b-vector in lattice coordinates
vkql(:)=vkl(:,ikp)+vqpl(:)

! equivalent reduced k-points for k and k+q (This gives problems now, GS must be runned with reducek=0. Should be fixed.)
!call findkpt(vkl(:,ikp),isym,ekp)
!call findkpt(vkql,isym,ekpq)

! generate the wavefunctions for all states at k+q
allocate(wfmtq1(lmmaxvr,nrcmtmax,natmtot,nspinor,wann_nband))
allocate(wfirq(ngtot,nspinor,wann_nband))
call genwfsvp(.false.,.false.,wann_nband,wann_bands,vkql,wfmtq1,ngtot,wfirq)

allocate(wfmtq(lmmaxvr,nrcmtmax,natmtot,nspinor,wann_nband))
do ist=1,wann_nband
  do is=1,nspecies
    nrc=nrcmt(is)
    nrci=nrcmtinr(is)
    do ia=1,natoms(is)
      ias=idxas(ia,is)
        do ispin=1,nspinor
          ! multiply by exp(-ib.r) 
          call genzrmt1(nrc,nrci,phase(:,:,ias),wfmtq1(:,:,ias,ispin,ist),wfmtq(:,:,ias,ispin,ist))
        end do
    end do
  end do
end do

!Calculate matrix elements
do ist=1,wann_nband
  allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot),zrhoir(ngtot))
  allocate(zrhoig(ngrf,4))
    do jst=1,wann_nband
    i=0
      do jspin=1,nspinor
        do ispin=1,nspinor
          i=i+1

          call genzrho(.true.,.false.,wfmt(1:lmmaxvr,1:nrcmtmax,1:natmtot,ispin,ist), &
                                      wfir(1:ngtot,ispin,ist)                       , &
                                      wfmtq(1:lmmaxvr,1:nrcmtmax,1:natmtot,jspin,jst), &
                                      wfirq(1:ngtot,jspin,jst)                       , &
                                      zrhomt,zrhoir)

          !zrhoig for all G rf (this should be optimized to calculate zrhoig only for G=0)
          call zftzf(1,gqc,ylmgq,ngrf,sfacgq,zrhomt,zrhoir,zrhoig(:,i))

          !Choose G=0, actual G vector already included in bqvec (therefore in wavefunction)
          matel_mmn(ist,ispin,jst,jspin) = zrhoig(1,i) * omega

        enddo
      enddo
    enddo
  deallocate(zrhomt,zrhoir)
  deallocate(zrhoig)

enddo

deallocate(wfmtq,wfirq)


end subroutine genw90mmn
!EOC

