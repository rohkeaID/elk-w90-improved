
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modpw

!----------------------------!
!     H-vector variables     !
!----------------------------!
! reduceh is .true. if the H-vectors are reduced with the crystal symmetries
logical reduceh
! H-vector cut-off for interstitial potential and density
real(8) hmaxvr
! H-vector grid sizes
integer ngridh(3)
! total number of H-vectors
integer nhtot
! integer grid intervals for each direction
integer inthv(2,3)
! number of H-vectors with |H| < hmaxvr
integer nhvec
! H-vector integer coordinates (i1,i2,i3)
integer, allocatable :: ivh(:,:)
! H-vector multiplicity after symmetry reduction
integer, allocatable :: mulh(:)
! H-vectors in Cartesian coordinates
real(8), allocatable :: vhc(:,:)
! length of H-vectors
real(8), allocatable :: hc(:)
! H-vector transformation matrix
real(8) vhmat(3,3)

!------------------------------!
!     H+k-vector variables     !
!------------------------------!
! maximum |H+k| cut-off for plane wave
real(8) hkmax
! number of H+k-vectors for plane waves
integer, allocatable :: nhk(:,:)
! maximum number of H+k-vectors over all k-points
integer nhkmax
! index from H+k-vectors to G-vectors
integer, allocatable :: ihkig(:,:,:)
! H+k-vectors in lattice coordinates
real(8), allocatable :: vhkl(:,:,:,:)
! H+k-vectors in Cartesian coordinates
real(8), allocatable :: vhkc(:,:,:,:)
! length of H+k-vectors
real(8), allocatable :: hkc(:,:,:)
! (theta, phi) coordinates of H+k-vectors
real(8), allocatable :: tphkc(:,:,:,:)
! structure factors for the H+k-vectors
complex(8), allocatable :: sfachk(:,:,:,:)

end module

