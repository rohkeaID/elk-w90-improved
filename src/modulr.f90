
! Copyright (C) 2017 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modulr

!-----------------------------!
!     ultracell variables     !
!-----------------------------!
! ultracell lattice vectors stored column-wise
real(8) avecu(3,3)
! ultracell reciprocal lattice vectors
real(8) bvecu(3,3)
! ultracell volume and Brillouin zone volume
real(8) omegau,omegabzu
! original number of k-points
integer nkpt0
! kappa-point cut-off
real(8) kpamax
! kappa-point grid sizes
integer ngridkpa(3)
! number of kappa-points
integer nkpa

!------------------------------!
!     G+Q-vector variables     !
!------------------------------!
! |G+Q| for all G+Q-vectors
real(8), allocatable :: gqc(:,:)

!---------------------------------------------------!
!     ultra long-range densities and potentials     !
!---------------------------------------------------!
! R-dependent densities and potentials
! charge density
real(8), allocatable :: rhormt(:,:,:),rhorir(:,:)
! magnetisation vector field
real(8), allocatable :: magrmt(:,:,:,:),magrir(:,:,:)
! Kohn-Sham potential
real(8), allocatable :: vsrmt(:,:,:),vsrir(:,:)
! phase factor functions exp(iQ.r) in each muffin-tin
complex(8), allocatable :: expmtq(:,:,:)
! Q-dependent densities and potentials
complex(8), allocatable :: rhoqmt(:,:,:),rhoqir(:,:)
complex(8), allocatable :: magqmt(:,:,:,:),magqir(:,:,:)

!----------------------------------------------!
!     eigenvalue and eigenvector variables     !
!----------------------------------------------!
! energy cut-off around Fermi energy of the second-variational states used for
! the ultra long-range basis
real(8) emaxulr
! maximum number of second-variational states over all k-points
integer nsvukmax
! number of long-range second-variational states for each k-point
integer, allocatable :: nsvuk(:)
! index from long-range basis states to second-variational states
integer, allocatable :: istuk(:,:)
! number of ultra long-range states
integer nstulr
! long-range occupation numbers
real(8), allocatable :: occulr(:,:)

end module

