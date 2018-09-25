
! Copyright (C) 2016 A. Davydov, A. Sanna, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine gwsefm
use modmain
use modgw
use modmpi
use modomp
implicit none
! local variables
integer ik,nthd
! allocatable arrays
real(8), allocatable :: vmt(:,:),vir(:)
real(8), allocatable :: bmt(:,:,:),bir(:,:)
complex(8), allocatable :: se(:,:,:)
! initialise universal variables
call init0
call init1
call init2
call init3
! read density and potentials from file
call readstate
! generate the core wavefunctions and densities
call gencore
! read Fermi energy from file
call readfermi
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! get the eigenvalues and occupancies from file
do ik=1,nkpt
  call getevalsv(filext,ik,vkl(:,ik),evalsv(:,ik))
  call getoccsv(filext,ik,vkl(:,ik),occsv(:,ik))
end do
! write the momentum matrix elements in the second-variational basis to file
call genpmat(.false.,.true.)
! generate the inverse dielectric function and write to file
call epsinv
! compute the matrix elements of -V_xc and -B_xc
allocate(vmt(npcmtmax,natmtot),vir(ngtot))
if (spinpol) then
  allocate(bmt(npcmtmax,natmtot,ndmag),bir(ngtot,ndmag))
end if
call gwlocal(vmt,vir,bmt,bir)
if (mp_mpi) write(*,*)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
! loop over reduced k-point set
call omp_hold(nkpt/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(se) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  allocate(se(nstsv,nstsv,0:nwfm))
!$OMP CRITICAL(gwsefm_)
  write(*,'("Info(gwsefm): ",I6," of ",I6," k-points")') ik,nkpt
!$OMP END CRITICAL(gwsefm_)
! determine the self-energy at the fermionic frequencies for current k-point
  call gwsefmk(ik,vmt,vir,bmt,bir,se)
! write the self-energy to file
  call putgwsefm(ik,se)
  deallocate(se)
end do
!$OMP END DO
!$OMP END PARALLEL
call omp_free(nthd)
if (mp_mpi) then
  write(*,*)
  write(*,'("Info(gwsefm): GW self-energy at the fermionic frequencies &
   &written to GWSEFM.OUT")')
end if
deallocate(vmt,vir)
if (spinpol) deallocate(bmt,bir)
return
end subroutine

