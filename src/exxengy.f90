
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine exxengy
use modmain
use modmpi
use modomp
implicit none
! local variables
integer ik,ist,jst,is,ia
integer nrc,nrci,npc
integer m1,m2,nthd
complex(8) z1
! allocatable arrays
complex(8), allocatable :: wfcr1(:,:),wfcr2(:,:)
complex(8), allocatable :: zrhomt(:),zvclmt(:),zfmt(:)
! external functions
complex(8) zfmtinp
external zfmtinp
allocate(wfcr1(npcmtmax,2),wfcr2(npcmtmax,2))
allocate(zrhomt(npcmtmax),zvclmt(npcmtmax),zfmt(npcmtmax))
! zero the exchange energy
engyx=0.d0
!--------------------------------------------------!
!     val-val-val and val-cr-val contributions     !
!--------------------------------------------------!
call omp_hold(nkpt/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
!$OMP CRITICAL(exxengy_)
  write(*,'("Info(exxengy): ",I6," of ",I6," k-points")') ik,nkpt
!$OMP END CRITICAL(exxengy_)
  call exxengyk(ik)
end do
!$OMP END DO
!$OMP END PARALLEL
call omp_free(nthd)
! add energies from each process and redistribute
call mpi_allreduce(mpi_in_place,engyx,1,mpi_double_precision,mpi_sum,mpicom, &
 ierror)
!-----------------------------------!
!    core-core-core contribution    !
!-----------------------------------!
! begin loops over atoms and species
do is=1,nspecies
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  npc=npcmt(is)
  do ia=1,natoms(is)
    do jst=1,nstsp(is)
      if (spcore(jst,is)) then
        do m2=-ksp(jst,is),ksp(jst,is)-1
! generate the core wavefunction in spherical coordinates (pass in m-1/2)
          call wavefcr(.false.,lradstp,is,ia,jst,m2,npcmtmax,wfcr2)
          do ist=1,nstsp(is)
            if (spcore(ist,is)) then
              do m1=-ksp(ist,is),ksp(ist,is)-1
                call wavefcr(.false.,lradstp,is,ia,ist,m1,npcmtmax,wfcr1)
! calculate the complex overlap density
                call zrho2(npc,wfcr1(:,1),wfcr1(:,2),wfcr2(:,1),wfcr2(:,2),zfmt)
                call zfsht(nrc,nrci,zfmt,zrhomt)
! calculate the Coulomb potential
                call zpotclmt(nrc,nrci,rcmt(:,is),zrhomt,zvclmt)
                z1=zfmtinp(nrc,nrci,rcmt(:,is),r2cmt(:,is),zrhomt,zvclmt)
                engyx=engyx-0.5d0*dble(z1)
              end do
! end loop over ist
            end if
          end do
        end do
! end loop over jst
      end if
    end do
! end loops over atoms and species
  end do
end do
deallocate(wfcr1,wfcr2,zrhomt,zvclmt,zfmt)
return

contains

subroutine zrho2(n,x1,x2,y1,y2,z)
implicit none
integer, intent(in) :: n
complex(8), intent(in) :: x1(n),x2(n),y1(n),y2(n)
complex(8), intent(out) :: z(n)
z(:)=conjg(x1(:))*y1(:)+conjg(x2(:))*y2(:)
return
end subroutine

end subroutine

