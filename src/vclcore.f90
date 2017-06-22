
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine vclcore(wfmt,vmat)
use modmain
implicit none
! arguments
complex(8), intent(in) :: wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv)
complex(8), intent(inout) :: vmat(nstsv,nstsv)
! local variables
integer ist1,ist2,ist3,m
integer is,ia,ias,nrc,nrci
complex(8) z1
! allocatable arrays
complex(8), allocatable :: zrhomt(:,:,:)
complex(8), allocatable :: wfcr(:,:,:),zfmt(:,:)
! external functions
complex(8) zfmtinp
external zfmtinp
allocate(zrhomt(lmmaxvr,nrcmtmax,nstsv))
allocate(wfcr(lmmaxvr,nrcmtmax,2))
vmat(:,:)=0.d0
do is=1,nspecies
  nrc=nrcmt(is)
  nrci=nrcmtinr(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ist3=1,nstsp(is)
      if (spcore(ist3,is)) then
        do m=-ksp(ist3,is),ksp(ist3,is)-1
! generate the core wavefunction in spherical coordinates (pass in m-1/2)
          call wavefcr(.false.,lradstp,is,ia,ist3,m,nrcmtmax,wfcr)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(zfmt)
!$OMP DO
          do ist1=1,nstsv
            allocate(zfmt(lmmaxvr,nrcmtmax))
! calculate the complex overlap density in spherical harmonics
            if (spinpol) then
              call genzrmt2(nrc,nrci,wfcr(:,:,1),wfcr(:,:,2), &
               wfmt(:,:,ias,1,ist1),wfmt(:,:,ias,2,ist1),zfmt)
            else
              call genzrmt1(nrc,nrci,wfcr(:,:,1),wfmt(:,:,ias,1,ist1),zfmt)
            end if
            call zfsht(nrc,nrci,zfmt,zrhomt(:,:,ist1))
            deallocate(zfmt)
          end do
!$OMP END DO
!$OMP END PARALLEL
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(zfmt,ist1,z1)
!$OMP DO
          do ist2=1,nstsv
            allocate(zfmt(lmmaxvr,nrcmtmax))
            call zpotclmt(nrc,nrci,rcmt(:,is),zrhomt(:,:,ist2),zfmt)
            do ist1=1,ist2
              z1=zfmtinp(nrc,nrci,rcmt(:,is),r2cmt(:,is),zrhomt(:,:,ist1),zfmt)
              vmat(ist1,ist2)=vmat(ist1,ist2)-z1
            end do
            deallocate(zfmt)
          end do
!$OMP END DO
!$OMP END PARALLEL
        end do
      end if
    end do
  end do
end do
! set the lower triangular part of the matrix
do ist1=1,nstsv
  do ist2=1,ist1-1
    vmat(ist1,ist2)=conjg(vmat(ist2,ist1))
  end do
end do
! scale the Coulomb matrix elements in the case of a hybrid functional
if (hybrid) vmat(:,:)=hybridc*vmat(:,:)
deallocate(zrhomt,wfcr)
return
end subroutine

