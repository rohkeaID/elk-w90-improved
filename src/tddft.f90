
! Copyright (C) 2014 K. Krieger, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine tddft
use modmain
use modtddft
use moddftu
use modmpi
use modomp
use modstore
use modtest
implicit none
! local variables
logical exist
integer ik,itimes0
! allocatable arrays
real(8), allocatable :: evalfv(:,:)
complex(8), allocatable :: evecfv(:,:,:),evecsv(:,:)
if (tshift) then
  write(*,*)
  write(*,'("Error(tddft): use tshift = .false. for the ground-state run")')
  write(*,*)
  stop
end if
! average force can be non-zero (allow for translation of atomic basis)
tfav0_=tfav0
tfav0=.false.
! initialise global variables
call init0
call init1
! read the charge density and potentials from file
call readstate
! generate the first- and second-variational eigenvectors and eigenvalues for
! the k-point set reduced with the symmetries which leave A(t) invariant for all
! time steps
call genvsig
call gencore
call readfermi
call linengy
call genapwfr
call genlofr
call olprad
call hmlrad
call gensocfr
call genevfsv
call occupy
! DFT+U
if (dftu.ne.0) then
  call gendmatmt
  call genvmatmt
end if
! generate the kinetic matrix elements in the second-variational basis
call genkmat(.false.,.false.)
! write the momentum matrix elements in first- and second-variational basis
call genpmat(.true.,.true.)
! write the power density to file
if (mp_mpi) call writeafpdt
! copy EVALFV.OUT, EVECFV.OUT, OCCSV.OUT and EVECSV.OUT to _TD.OUT extension
if (mp_mpi.and.(task.eq.460)) then
  allocate(evalfv(nstfv,nspnfv),evecfv(nmatmax,nstfv,nspnfv))
  allocate(evecsv(nstsv,nstsv))
  do ik=1,nkpt
    call getevalfv('.OUT',ik,vkl(:,ik),evalfv)
    call putevalfv('_TD.OUT',ik,evalfv)
    call getevecfv('.OUT',ik,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
    call putevecfv('_TD.OUT',ik,evecfv)
    call putoccsv('_TD.OUT',ik,occsv(:,ik))
    call getevecsv('.OUT',ik,vkl(:,ik),evecsv)
! randomise eigenvectors at t=0 if required
    call rndevsv(rndevt0,evecsv)
    call putevecsv('_TD.OUT',ik,evecsv)
  end do
  deallocate(evalfv,evecfv,evecsv)
end if
! set global file extension
filext='_TD.OUT'
! output the new k-point set to file
if (mp_mpi) call writekpts
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
itimes0=0
! restart if required
if (task.eq.461) call readtimes(itimes0)
! set the stop signal to .false.
tstop=.false.
!---------------------------------!
!    main loop over time steps    !
!---------------------------------!
if (mp_mpi) write(*,*)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
do itimes=itimes0+1,ntimes-1
  if (mp_mpi) then
    write(*,'("Info(tddft): time step ",I8," of ",I8,",   t = ",G18.10)') &
     itimes,ntimes,times(itimes)
  end if
! reset the OpenMP thread variables
  call omp_reset
! evolve the wavefunctions across a single time step
  call timestep
! generate the density and magnetisation at current time step
  call rhomag
! compute the total current
  call curden(afieldt(:,itimes))
! compute the time-dependent Kohn-Sham potentials and magnetic fields
  call potkst
! DFT+U
  if (dftu.ne.0) then
    call gendmatmt
    call genvmatmt
  end if
! compute the atomic forces if required
  if (tforce) call force
  if (mp_mpi) then
! write TDDFT output
    call writetddft
! write the time step to file
    call writetimes
! check for STOP file
    inquire(file='STOP',exist=exist)
    if (exist) then
      open(50,file='STOP')
      close(50,status='DELETE')
      tstop=.true.
    end if
  end if
! broadcast tstop from master process to all other processes
  call mpi_bcast(tstop,1,mpi_logical,0,mpicom,ierror)
  if (tstop) exit
end do
filext='.OUT'
tfav0=tfav0_
! write the total current of the last step to test file
call writetest(460,'total current of last time step',nv=3,tol=1.d-4,rva=curtot)
return
end subroutine

