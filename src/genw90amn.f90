! Copyright (C) 2015 Jon Lafuente and Manh Duc Le
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genw90amn
! !INTERFACE:
subroutine genw90amn(wfmt,wfir,wfmtq,wfirq,gqc,ylmgq,sfacgq,matel_amn)
! !USES:
use modmain
use modrandom
use modw90
! !DESCRIPTION:
!   Calculates the Amn overlap matrices between Bloch states and 
!   trial orbitals required for Wannier90.
!EOP
!BOC
implicit none
! local variables
complex(8), intent(in) :: wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,wann_nband)
complex(8), intent(in) :: wfir(ngtot,nspinor,wann_nband)
complex(8), intent(in) :: wfmtq(lmmaxvr,nrcmtmax,natmtot,nspinor,wann_nproj)
complex(8), intent(in) :: wfirq(ngtot,nspinor,wann_nproj)
real(8), intent(in) :: gqc(ngrf)
complex(8), intent(in) :: ylmgq(lmmaxvr,ngrf)
complex(8), intent(in) :: sfacgq(ngrf,natmtot)
complex(8), intent(inout) :: matel_amn(wann_nband,4,wann_nproj,4)
! local variables
integer :: ist,jst,ispin,jspin
complex(8) :: matel
real(8) :: norm
! allocatable arrays
complex(8), allocatable :: zrhomt(:,:,:),zrhoir(:)
complex(8), allocatable :: zrhoig(:)

allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot),zrhoir(ngtot))
allocate(zrhoig(ngrf))
do jst=1,wann_nproj
  if(wann_proj_isrand(jst)) then
! if this projection should be random (at present, any non-atom centred projection!)
    norm=0.
    do jspin=1,nspinor
      do ispin=1,nspinor
        do ist=1,wann_nband
          matel=cmplx(randomu(),randomu(),kind=8)
          norm=norm+abs(matel)
          matel_amn(ist,ispin,jst,jspin)=matel
        end do
      end do
    end do
    do jspin=1,nspinor
      do ispin=1,nspinor
        do ist=1,wann_nband
          matel_amn(ist,ispin,jst,jspin)=matel_amn(ist,ispin,jst,jspin)/norm
        end do
      end do
    end do
  else
! otherwise calculate the overlap matrix element (integral) directly.
    do ist=1,wann_nband
      do jspin=1,nspinor
        do ispin=1,nspinor
          call genzrho(.true.,.false.,wfmt(1:lmmaxvr,1:nrcmtmax,1:natmtot,ispin,ist), &
                                      wfir(1:ngtot,ispin,ist)                       , &
                                      wfmtq(1:lmmaxvr,1:nrcmtmax,1:natmtot,jspin,jst), &
                                      wfirq(1:ngtot,jspin,jst)                       , &
                                      zrhomt,zrhoir)
          !zrhoig for all G rf (this should be optimized to calculate zrhoig only for G=0)
          call zftzf(1,gqc,ylmgq,ngrf,sfacgq,zrhomt,zrhoir,zrhoig)
          !Choose G=0, actual G vector already included in bqvec (therefore in wavefunction)
          matel_amn(ist,ispin,jst,jspin) = zrhoig(1) * omega
        enddo
      enddo
    enddo
  end if
enddo
deallocate(zrhomt,zrhoir)
deallocate(zrhoig)

end subroutine
!EOC
