
! Copyright (C) 2014 K. Krieger, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine tddft
use modmain
use modtddft
use moddftu
use modmpi
implicit none
! local variables
integer ik,itimes0
real(8) t1
! allocatable arrays
complex(8), allocatable :: evecfv(:,:,:),evecsv(:,:)
if (symtype.ne.0) then
  write(*,*)
  write(*,'("Error(tddft): use nosym = .true. for the ground-state run")')
  write(*,*)
  stop
end if
t1=sum(abs(vkloff(:)))
if (t1.gt.epslat) then
  write(*,*)
  write(*,'("Warning(tddft): non-zero vkloff may cause inaccuracies")')
end if
! initialise global variables
call init0
call init1
! read the charge density and potentials from file
call readstate
! generate the core wavefunctions and densities
call gencore
! read Fermi energy from file
call readfermi
! find the linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! get the eigenvalues and occupation numbers from file
do ik=1,nkpt
  call getevalsv(filext,vkl(:,ik),evalsv(:,ik))
  call getoccsv(filext,vkl(:,ik),occsv(:,ik))
end do
! DFT+U
if (dftu.ne.0) then
  call gendmatmt
  call genvmatmt
end if
! generate the kinetic matrix elements in the second-variational basis
call genkmat(.false.,.false.)
! write the momentum matrix elements in first- and second-variational basis
call genpmat(.true.,.true.)
! read time-dependent A-field from file
call readafieldt
! write the power density to file
if (mp_mpi) call writeafpdt
! copy OCCSV.OUT, EVECFV.OUT and EVECSV.OUT to _TD.OUT extension
if (mp_mpi.and.(task.eq.460)) then
  do ik=1,nkpt
    call putoccsv('_TD.OUT',ik,occsv(:,ik))
    allocate(evecfv(nmatmax,nstfv,nspnfv))
    call getevecfv('.OUT',vkl(:,ik),vgkl(:,:,:,ik),evecfv)
    call putevecfv('_TD.OUT',ik,evecfv)
    deallocate(evecfv)
    allocate(evecsv(nstsv,nstsv))
    call getevecsv('.OUT',vkl(:,ik),evecsv)
    call putevecsv('_TD.OUT',ik,evecsv)
    deallocate(evecsv)
  end do
end if
! set global file extension
filext='_TD.OUT'
! synchronise MPI processes
call mpi_barrier(mpi_comm_kpt,ierror)
itimes0=0
! restart if required
if (task.eq.461) call readtimes(itimes0)
!---------------------------------!
!    main loop over time steps    !
!---------------------------------!
if (mp_mpi) write(*,*)
do itimes=itimes0+1,ntimes-1
  if (mp_mpi) then
    write(*,'("Info(tddft): time step ",I8," of ",I8,",   t = ",G18.10)') &
     itimes,ntimes,times(itimes)
  end if
! generate the density and magnetisation at current time step
  call rhomag
! compute the total current
  call current
! DFT+U
  if (dftu.ne.0) then
    call gendmatmt
    call genvmatmt
  end if
! compute the time-dependent Kohn-Sham potentials and magnetic fields
  call potkst
! evolve the wavefunctions across a single time step
  call timestep
  if (mp_mpi) then
! write TDDFT output
    call writetddft
! write the time step to file
    call writetimes
  end if
! synchronise MPI processes
  call mpi_barrier(mpi_comm_kpt,ierror)
end do
filext='.OUT'
return
end subroutine

