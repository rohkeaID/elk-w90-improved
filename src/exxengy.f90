
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine exxengy
use modmain
use modmpi
implicit none
! local variables
integer ik,ist,jst,is,ia
integer nrc,nrci,m1,m2
complex(8) z1
! allocatable arrays
complex(8), allocatable :: wfcr1(:,:,:),wfcr2(:,:,:)
complex(8), allocatable :: zrhomt(:,:),zvclmt(:,:),zfmt(:,:)
! external functions
complex(8) zfmtinp
external zfmtinp
allocate(wfcr1(lmmaxvr,nrcmtmax,2),wfcr2(lmmaxvr,nrcmtmax,2))
allocate(zrhomt(lmmaxvr,nrcmtmax),zvclmt(lmmaxvr,nrcmtmax))
allocate(zfmt(lmmaxvr,nrcmtmax))
! zero the exchange energy
engyx=0.d0
!--------------------------------------------------!
!     val-val-val and val-cr-val contributions     !
!--------------------------------------------------!
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
!$OMP CRITICAL
  write(*,'("Info(exxengy): ",I6," of ",I6," k-points")') ik,nkpt
!$OMP END CRITICAL
  call exxengyk(ik)
end do
!$OMP END DO
!$OMP END PARALLEL
! add energies from each process and redistribute
call mpi_allreduce(mpi_in_place,engyx,1,mpi_double_precision,mpi_sum, &
 mpi_comm_kpt,ierror)
!-----------------------------------!
!    core-core-core contribution    !
!-----------------------------------!
! begin loops over atoms and species
do is=1,nspecies
  nrc=nrcmt(is)
  nrci=nrcmtinr(is)
  do ia=1,natoms(is)
    do jst=1,nstsp(is)
      if (spcore(jst,is)) then
        do m2=-ksp(jst,is),ksp(jst,is)-1
! generate the core wavefunction in spherical coordinates (pass in m-1/2)
          call wavefcr(.false.,lradstp,is,ia,jst,m2,nrcmtmax,wfcr2)
          do ist=1,nstsp(is)
            if (spcore(ist,is)) then
              do m1=-ksp(ist,is),ksp(ist,is)-1
                call wavefcr(.false.,lradstp,is,ia,ist,m1,nrcmtmax,wfcr1)
! calculate the complex overlap density
                call genzrmt2(nrc,nrci,wfcr1(:,:,1),wfcr1(:,:,2),wfcr2(:,:,1), &
                 wfcr2(:,:,2),zfmt)
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
end subroutine

