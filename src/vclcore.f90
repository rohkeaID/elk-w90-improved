
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine vclcore(wfmt,vmat)
use modmain
use modomp
implicit none
! arguments
complex(8), intent(in) :: wfmt(npcmtmax,natmtot,nspinor,nstsv)
complex(8), intent(inout) :: vmat(nstsv,nstsv)
! local variables
integer ist1,ist2,ist3
integer is,ia,ias,m,nthd
integer nrc,nrci,npc
complex(8) z1
! allocatable arrays
complex(8), allocatable :: zrhomt(:,:),wfcr(:,:),zfmt(:)
! external functions
complex(8) zfmtinp
external zfmtinp
allocate(zrhomt(npcmtmax,nstsv),wfcr(npcmtmax,2))
do is=1,nspecies
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  npc=npcmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ist3=1,nstsp(is)
      if (spcore(ist3,is)) then
        do m=-ksp(ist3,is),ksp(ist3,is)-1
! generate the core wavefunction in spherical coordinates (pass in m-1/2)
          call wavefcr(.false.,lradstp,is,ia,ist3,m,npcmtmax,wfcr)
          call omp_hold(nstsv,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zfmt) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
          do ist1=1,nstsv
            allocate(zfmt(npcmtmax))
! calculate the complex overlap density in spherical harmonics
            if (spinpol) then
              call zrho2(npc,wfcr(:,1),wfcr(:,2),wfmt(:,ias,1,ist1), &
               wfmt(:,ias,2,ist1),zfmt)
            else
              call zrho1(npc,wfcr(:,1),wfmt(:,ias,1,ist1),zfmt)
            end if
            call zfsht(nrc,nrci,zfmt,zrhomt(:,ist1))
            deallocate(zfmt)
          end do
!$OMP END DO
!$OMP END PARALLEL
          call omp_free(nthd)
          call omp_hold(nstsv,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zfmt,ist1,z1) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
          do ist2=1,nstsv
            allocate(zfmt(npcmtmax))
            call zpotclmt(nrc,nrci,rcmt(:,is),zrhomt(:,ist2),zfmt)
            do ist1=1,ist2
              z1=zfmtinp(nrc,nrci,rcmt(:,is),r2cmt(:,is),zrhomt(:,ist1),zfmt)
              vmat(ist1,ist2)=vmat(ist1,ist2)-z1
            end do
            deallocate(zfmt)
          end do
!$OMP END DO
!$OMP END PARALLEL
          call omp_free(nthd)
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

