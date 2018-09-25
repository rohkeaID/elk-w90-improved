
! Copyright (C) 2017 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gndstulr
use modmain
use modulr
implicit none
! local variables
integer ik0
! allocatable arrays
complex(8), allocatable :: evecu(:,:)
! initialise global variables
call init0
call init1
call initulr
! read the charge density from file
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
! generate the first- and second-variational eigenvectors and eigenvalues
call genevfsv
! generate the index to the ultra long-range states
call genidxulr
! set the stop signal to .false.
tstop=.false.
! begin the self-consistent loop
do iscl=1,maxscl
! zero the density
  rhormt(:,:,:)=0.d0
  rhorir(:,:)=0.d0
! zero the magnetisation
  if (spinpol) then
    magrmt(:,:,:,:)=0.d0
    magrir(:,:,:)=0.d0
  end if
! loop over original k-points
  do ik0=1,nkpt0
    allocate(evecu(nstulr,nstulr))
! solve the ultra long-range eigenvalue equation
    call eveqnulr(ik0)

!*********
occulr(:,:)=1.d0
evecu(:,:)=1.d0
!*********


! add to the density, magnetisation and current
    call rhomaguk(ik0,evecu)
    deallocate(evecu)
  end do
end do

return
end subroutine

