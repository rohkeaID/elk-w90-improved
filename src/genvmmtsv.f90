
! Copyright (C) 2015 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genvmmtsv(wfmt,vmat)
use modmain
use moddftu
implicit none
! arguments
complex(8), intent(in) :: wfmt(npcmtmax,natmtot,nspinor,nstsv)
complex(8), intent(inout) :: vmat(nstsv,nstsv)
! local variables
integer ist,jst,ispn,jspn
integer is,ias,ld
integer nrc,nrci,nrco
integer nm,l,lm
integer npc,npci,i
! allocatable arrays
complex(8), allocatable :: wfmt1(:,:,:),wfmt2(:)
! external functions
complex(8) zfmtinp
external zfmtinp
if (.not.tvmatmt) return
ld=lmmaxdm*nspinor
allocate(wfmt1(npcmtmax,nspinor,nstsv),wfmt2(npcmtmax))
! loop over atoms
do ias=1,natmtot
  is=idxis(ias)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  nrco=nrc-nrci
  npc=npcmt(is)
  npci=npcmti(is)
! convert wavefunctions to spherical harmonics
  do ist=1,nstsv
    do ispn=1,nspinor
      call zfsht(nrc,nrci,wfmt(:,ias,ispn,ist),wfmt1(:,ispn,ist))
    end do
  end do
! loop over second-variational states
  do jst=1,nstsv
! loop over spins
    do ispn=1,nspinor
      wfmt2(1:npc)=0.d0
      do jspn=1,nspinor
        do l=0,lmaxdm
          if (tvmmt(l,ias)) then
            nm=2*l+1
            lm=idxlm(l,-l)
            if (l.le.lmaxi) then
              call zgemm('N','N',nm,nrci,nm,zone,vmatmt(lm,ispn,lm,jspn,ias), &
               ld,wfmt1(lm,jspn,jst),lmmaxi,zone,wfmt2(lm),lmmaxi)
            end if
            i=npci+lm
            call zgemm('N','N',nm,nrco,nm,zone,vmatmt(lm,ispn,lm,jspn,ias),ld, &
             wfmt1(i,jspn,jst),lmmaxo,zone,wfmt2(i),lmmaxo)
          end if
        end do
      end do
! compute the matrix elements
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
      do ist=1,nstsv
        vmat(ist,jst)=vmat(ist,jst)+zfmtinp(nrc,nrci,rcmt(:,is),r2cmt(:,is), &
         wfmt1(:,ispn,ist),wfmt2)
      end do
!$OMP END DO
!$OMP END PARALLEL
    end do
  end do
end do
deallocate(wfmt1,wfmt2)
return
end subroutine

