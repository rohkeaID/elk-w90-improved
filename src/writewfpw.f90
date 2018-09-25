
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writewfpw
use modmain
use modpw
use modmpi
use modomp
implicit none
! local variables
integer ik,recl,nthd
! allocatable arrays
complex(8), allocatable :: wfpw(:,:,:)
! initialise global variables
call init0
call init1
call init4
! read density and potentials from file
call readstate
! read Fermi energy from file
call readfermi
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! delete existing WFPW.OUT
if (mp_mpi) then
  open(170,file='WFPW.OUT')
  close(170,status='DELETE')
end if
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
! determine the record length and open WFPW.OUT
allocate(wfpw(nhkmax,nspinor,nstsv))
inquire(iolength=recl) vkl(:,1),nhkmax,nspinor,nstsv,wfpw
deallocate(wfpw)
open(170,file='WFPW.OUT',form='UNFORMATTED',access='DIRECT',recl=recl)
! begin parallel loop over k-points
call omp_hold(nkpt/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(wfpw) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  allocate(wfpw(nhkmax,nspinor,nstsv))
!$OMP CRITICAL(writewfpw_)
  write(*,'("Info(writewfpw): ",I6," of ",I6," k-points")') ik,nkpt
!$OMP END CRITICAL(writewfpw_)
! generate the plane wave wavefunctions
  call genwfpw(vkl(:,ik),ngk(:,ik),igkig(:,:,ik),vgkl(:,:,:,ik), &
   vgkc(:,:,:,ik),gkc(:,:,ik),tpgkc(:,:,:,ik),sfacgk(:,:,:,ik),nhk(:,ik), &
   vhkc(:,:,:,ik),hkc(:,:,ik),tphkc(:,:,:,ik),sfachk(:,:,:,ik),wfpw)
!$OMP CRITICAL(u170)
  write(170,rec=ik) vkl(:,ik),nhkmax,nspinor,nstsv,wfpw
!$OMP END CRITICAL(u170)
  deallocate(wfpw)
end do
!$OMP END DO
!$OMP END PARALLEL
call omp_free(nthd)
close(170)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
if (mp_mpi) then
  write(*,*)
  write(*,'("Info(writewfpw): plane wave wavefunctions written to WFPW.OUT")')
  write(*,'(" for all H+k-vectors up to |H+k| < hkmax")')
end if
return
end subroutine

