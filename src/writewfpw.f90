
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writewfpw
use modmain
use modpw
use modmpi
implicit none
! local variables
integer ik,recl
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
  open(50,file='WFPW.OUT')
  close(50,status='DELETE')
end if
! synchronise MPI processes
call mpi_barrier(mpi_comm_kpt,ierror)
! determine the record length and open WFPW.OUT
allocate(wfpw(nhkmax,nspinor,nstsv))
inquire(iolength=recl) vkl(:,1),nhkmax,nspinor,nstsv,wfpw
deallocate(wfpw)
open(50,file='WFPW.OUT',action='WRITE',form='UNFORMATTED',access='DIRECT', &
 recl=recl)
! begin parallel loop over k-points
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(wfpw)
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  allocate(wfpw(nhkmax,nspinor,nstsv))
!$OMP CRITICAL
  write(*,'("Info(writewfpw): ",I6," of ",I6," k-points")') ik,nkpt
!$OMP END CRITICAL
! generate the plane wave wavefunctions
  call genwfpw(vkl(:,ik),ngk(:,ik),igkig(:,:,ik),vgkl(:,:,:,ik), &
   vgkc(:,:,:,ik),gkc(:,:,ik),tpgkc(:,:,:,ik),sfacgk(:,:,:,ik),nhk(:,ik), &
   vhkc(:,:,:,ik),hkc(:,:,ik),tphkc(:,:,:,ik),sfachk(:,:,:,ik),wfpw)
!$OMP CRITICAL
  write(50,rec=ik) vkl(:,ik),nhkmax,nspinor,nstsv,wfpw
!$OMP END CRITICAL
  deallocate(wfpw)
end do
!$OMP END DO
!$OMP END PARALLEL
close(50)
! synchronise MPI processes
call mpi_barrier(mpi_comm_kpt,ierror)
if (mp_mpi) then
  write(*,*)
  write(*,'("Info(writewfpw): plane wave wavefunctions written to WFPW.OUT")')
  write(*,'(" for all H+k-vectors up to |H+k| < hkmax")')
end if
return
end subroutine

