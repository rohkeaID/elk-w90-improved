
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rdmminn
! !INTERFACE:
subroutine rdmminn
! !USES:
use modmain
use modrdm
use modmpi
! !DESCRIPTION:
!   Minimizes the total energy w.r.t. occupation numbers. The steepest-descent
!   algorithm is used.
!
! !REVISION HISTORY:
!   Created 2008 (Sharma)
!EOP
!BOC
implicit none
! local variables
integer it,n
if (maxitn.lt.1) return
! write the Coulomb matrix elements to file
call writevclijji
! synchronise MPI processes
call mpi_barrier(mpi_comm_kpt,ierror)
! calculate derivative of kinetic energy w.r.t. evecsv
call rdmdkdc
! begin iteration loop
do it=1,maxitn
  if (mp_mpi) then
    if (mod(it,10).eq.0) then
      write(*,'("Info(rdmminn): iteration ",I4," of ",I4)') it,maxitn
    end if
  end if
! generate the density and magnetisation
  call rhomag
! calculate the Coulomb potential
  call potcoul
! calculate Coulomb potential matrix elements
  call genvmat(vclmt,vclir,vclmat)
! update occupation numbers and write to file (MPI master process only)
  if (mp_mpi) call rdmvaryn
! broadcast occupation numbers to all other processes
  n=nstsv*nkpt
  call mpi_bcast(occsv,n,mpi_double_precision,0,mpi_comm_kpt,ierror)
! calculate the energy
  call rdmenergy
! write energy to file
  if (mp_mpi) then
    write(61,'(I6,G18.10)') it,engytot
    call flushifc(61)
  end if
! end iteration loop
end do
if (mp_mpi) then
  write(60,*)
  write(60,'("Occupation number minimisation done")')
  write(61,*)
  if (spinpol) write(63,*)
end if
return
end subroutine
!EOC
