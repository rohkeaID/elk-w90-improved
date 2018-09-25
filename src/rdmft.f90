
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rdmft
! !INTERFACE:
subroutine rdmft
! !USES:
use modmain
use modrdm
use modmpi
! !DESCRIPTION:
!  Main routine for one-body reduced density matrix functional theory (RDMFT).
!
! !REVISION HISTORY:
!   Created 2008 (Sharma)
!EOP
!BOC
implicit none
! initialise global variables
call init0
call init1
! generate q-point set and gclq array
call init2
! read density and potentials from file
call readstate
! Fourier transform Kohn-Sham potential to G-space
call genvsig
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
! compute the overlap radial integrals
call olprad
! compute the Hamiltonian radial integrals
call hmlrad
! generate the spin-orbit coupling radial functions
call gensocfr
! compute the kinetic energy of the core
call energykncr
! generate the first- and second-variational eigenvectors and eigenvalues
call genevfsv
! find the occupation numbers
call occupy
! generate the kinetic matrix elements in the first-variational basis
call genkmat(.true.,.false.)
! open information files (MPI master process only)
if (mp_mpi) then
  open(60,file='RDM_INFO.OUT',form='FORMATTED')
  open(61,file='RDMN_ENERGY.OUT',form='FORMATTED')
  open(62,file='RDMC_ENERGY.OUT',form='FORMATTED')
  if (spinpol) then
    open(63,file='RDMN_MOMENT.OUT',form='FORMATTED')
    open(64,file='RDMC_MOMENT.OUT',form='FORMATTED')
  end if
  open(65,file='RDM_ENERGY.OUT',form='FORMATTED')
! write out general information to RDM_INFO.OUT
  call writeinfo(60)
  write(60,*)
  write(60,'("+------------------------------+")')
  write(60,'("| Self-consistent loop started |")')
  write(60,'("+------------------------------+")')
end if
! begin main RDMFT self-consistent loop
do iscl=1,rdmmaxscl
  if (mp_mpi) then
    write(60,*)
    write(60,'("+--------------------+")')
    write(60,'("| Loop number : ",I4," |")') iscl
    write(60,'("+--------------------+")')
    flush(60)
    write(*,*)
    write(*,'("Info(rdmft): self-consistent loop number : ",I4)') iscl
  end if
! synchronise MPI processes
  call mpi_barrier(mpicom,ierror)
! minimisation over natural orbitals
  call rdmminc
! minimisation over occupation number
  call rdmminn
! compute the RDMFT 'eigenvalues'
  call rdmeval
  if (mp_mpi) then
    call rdmwriteengy(60)
    call writechg(60)
    call writeeval
! write out the total energy
    write(65,'(G18.10)') engytot
    flush(65)
  end if
end do
if (mp_mpi) then
  write(60,*)
  write(60,'("+------------------------------+")')
  write(60,'("| Self-consistent loop stopped |")')
  write(60,'("+------------------------------+")')
! write density to STATE.OUT
  call writestate
  write(60,*)
  write(60,'("Wrote STATE.OUT")')
  write(60,*)
  write(60,'("+----------------------------+")')
  write(60,'("| Elk version ",I1.1,".",I1.1,".",I2.2," stopped |")') version
  write(60,'("+----------------------------+")')
! close information files
  close(60)
  close(61)
  close(62)
  if (spinpol) then
    close(63)
    close(64)
  end if
  close(65)
end if
return
end subroutine
!EOC

