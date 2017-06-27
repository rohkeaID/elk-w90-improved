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
integer is,ia,ias,nrc,nrci
integer nr
real(8) t1real, t1img
complex(8) t1
! automatic arrays
complex(8) fr(nrcmtmax)
! external functions
real(8) fintgt
complex(8) zdotu
external fintgt,zdotu
! allocatable arrays
complex(8), allocatable :: zrhomt(:,:,:),zrhoir(:)
complex(8), allocatable :: zrhoig(:)
complex(8), allocatable :: amnmttot(:,:,:,:)
complex(8), allocatable :: amnir(:,:,:,:)
complex(8), allocatable :: cfunir_complex(:)

allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot),zrhoir(ngtot))
allocate(zrhoig(ngrf))
allocate(amnmttot(wann_nband,4,wann_nproj,4))
allocate(amnir(wann_nband,4,wann_nproj,4))
allocate(cfunir_complex(ngtot))

amnmttot = cmplx(0.0d0,0.0d0,kind=8)
cfunir_complex = cmplx(cfunir,0.0d0,kind=8)

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
    matel_amn(:,:,jst,:)=matel_amn(:,:,jst,:)/norm
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
          ! !zrhoig for all G rf (this should be optimized to calculate zrhoig only for G=0)
          ! call zftzf(1,gqc,ylmgq,ngrf,sfacgq,zrhomt,zrhoir,zrhoig)
          ! !Choose G=0, actual G vector already included in bqvec (therefore in wavefunction)
          ! matel_amn(ist,ispin,jst,jspin) = zrhoig(1) * omega

          ! find the muffin-tin amn
          do ias=1,natmtot
            is=idxis(ias)
            nr=nrcmt(is)
            fr(1:nr)=zrhomt(1,1:nr,ias)*r2cmt(1:nr,is)
            t1real=fintgt(-1,nrcmt(is),rcmt(:,is),dble(fr))
            t1img=fintgt(-1,nrcmt(is),rcmt(:,is),aimag(fr))
            t1=cmplx(t1real,t1img,kind=8)
            amnmttot(ist,ispin,jst,jspin)=amnmttot(ist,ispin,jst,jspin)+fourpi*y00*t1
          end do

          ! find the interstitial amn
          t1=zdotu(ngtot,zrhoir,1,cfunir_complex,1)
          amnir(ist,ispin,jst,jspin)=t1*omega/dble(ngtot)

        enddo
      enddo
    enddo
    ! total calculated amn
    matel_amn(:,:,jst,:)=amnmttot(:,:,jst,:)+amnir(:,:,jst,:)
  end if
enddo

deallocate(cfunir_complex)
deallocate(amnmttot,amnir)
deallocate(zrhomt,zrhoir)
deallocate(zrhoig)

end subroutine
!EOC
