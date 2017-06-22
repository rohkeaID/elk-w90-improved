
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rdmminc
! !INTERFACE:
subroutine rdmminc
! !USES:
use modmain
use modrdm
use modmpi
! !DESCRIPTION:
!   Minimizes the total energy with respect to the second-variational
!   coefficients {\tt evecsv}. The steepest-descent algorithm is used.
!
! !REVISION HISTORY:
!   Created 2008 (Sharma)
!EOP
!BOC
implicit none
! local variables
integer it
if (maxitc.lt.1) return
! begin iteration loop
do it=1,maxitc
  if (mp_mpi) then
    write(*,'("Info(rdmminc): iteration ",I4," of ",I4)') it,maxitc
  end if
! generate the density and magnetisation
  call rhomag
! calculate the Coulomb potential
  call potcoul
! calculate Coulomb potential matrix elements
  call genvmat(vclmt,vclir,vclmat)
! calculate derivative of kinetic energy w.r.t. evecsv
  call rdmdkdc
! write the Coulomb matrix elements to file
  call writevclijjk
! synchronise MPI processes
  call mpi_barrier(mpi_comm_kpt,ierror)
! update evecsv, orthogonalise and write to file (MPI master process only)
  if (mp_mpi) call rdmvaryc
! synchronise MPI processes
  call mpi_barrier(mpi_comm_kpt,ierror)
! calculate the energy
  call rdmenergy
! write energy to file
  if (mp_mpi) then
    write(62,'(I6,G18.10)') it,engytot
    call flushifc(62)
  end if
! end iteration loop
end do
if (mp_mpi) then
  write(60,*)
  write(60,'("Natural orbital minimisation done")')
  write(62,*)
  if (spinpol) write(64,*)
end if
return
end subroutine
!EOC
