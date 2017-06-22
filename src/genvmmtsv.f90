
! Copyright (C) 2015 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genvmmtsv(wfmt,vmat)
use modmain
use moddftu
implicit none
! arguments
complex(8), intent(in) :: wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv)
complex(8), intent(inout) :: vmat(nstsv,nstsv)
! local variables
integer ist,jst,ispn,jspn
integer is,ias,nrc,nrci
integer ld,nm,lm,l
! allocatable arrays
complex(8), allocatable :: wfmt1(:,:,:,:),wfmt2(:,:)
! external functions
complex(8) zfmtinp
external zfmtinp
if (.not.tvmatmt) return
ld=lmmaxdm*nspinor
allocate(wfmt1(lmmaxvr,nrcmtmax,nspinor,nstsv),wfmt2(lmmaxvr,nrcmtmax))
! loop over atoms
do ias=1,natmtot
  is=idxis(ias)
  nrc=nrcmt(is)
  nrci=nrcmtinr(is)
! convert wavefunctions to spherical harmonics
  do ist=1,nstsv
    do ispn=1,nspinor
      call zfsht(nrc,nrci,wfmt(:,:,ias,ispn,ist),wfmt1(:,:,ispn,ist))
      wfmt1(lmmaxinr+1:,1:nrci,ispn,ist)=0.d0
    end do
  end do
! loop over second-variational states
  do jst=1,nstsv
! loop over spins
    do ispn=1,nspinor
      call zfmtzero(nrc,nrci,wfmt2)
      do jspn=1,nspinor
        do l=0,lmaxdm
          if (tvmmt(l,ias)) then
            nm=2*l+1
            lm=idxlm(l,-l)
            call zgemm('N','N',nm,nrc,nm,zone,vmatmt(lm,ispn,lm,jspn,ias),ld, &
             wfmt1(lm,1,jspn,jst),lmmaxvr,zone,wfmt2(lm,1),lmmaxvr)
          end if
        end do
      end do
! compute the matrix elements
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
      do ist=1,nstsv
        vmat(ist,jst)=vmat(ist,jst)+zfmtinp(nrc,nrci,rcmt(:,is),r2cmt(:,is), &
         wfmt1(:,:,ispn,ist),wfmt2)
      end do
!$OMP END DO
!$OMP END PARALLEL
    end do
  end do
end do
deallocate(wfmt1,wfmt2)
return
end subroutine

