
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modscdft

! number of normal Kohn-Sham states to use in the BdG equations
integer nbdg
! size of the BdG matrix (2*nbdg)
integer nmbdg
! energy window around the Fermi energy containing the BdG states
real(8) ewbdg
! index from the BdG states to the normal Kohn-Sham states
integer, allocatable :: idxbdg(:,:)
! eigenvalues of the BdG Hamiltonian
real(8), allocatable :: evalbdg(:)
! BdG inversion algorithm mixing parameter
real(8) taubdg
! magnitude of random numbers used to initialise the anomalous density
real(8) rndachi

end module

