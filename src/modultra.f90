
! Copyright (C) 2016 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modultra

!-----------------------------!
!     ultracell variables     !
!-----------------------------!
! ultracell is .true. if an ultracell is to be used
logical ultracell
! kappa-point grid sizes
integer ngridkpa(3)
! ultracell lattice vectors stored column-wise
real(8) avecu(3,3)
! ultracell reciprocal lattice vectors
real(8) bvecu(3,3)
! ultracell volume and Brillouin zone volume
real(8) omegau,omegabzu
! original number of k-points
integer nkpt0
! number of kappa-points
integer nkpa
! integer grid intervals for kappa-points
integer intkpa(2,3)
! locations of kappa-points on integer grid
integer, allocatable :: ivkpa(:,:)
! index of kappa = 0 point
integer ikpa0
! kappa-points in ultracell lattice coordinates
real(8), allocatable :: vkpalu(:,:)
! kappa-points in unit cell lattice coordinates
real(8), allocatable :: vkpal(:,:)
! kappa-points in Cartesian coordinates
real(8), allocatable :: vkpac(:,:)
! integer grid intervals for the Q-points
integer intq(2,3)
! map from (i1,i2,i3) to Q-vector index
integer, allocatable :: ivqiq(:,:,:)
! scalar potential
real(8), allocatable :: vsu(:),vsuq(:)
! vector potential
real(8), allocatable :: asu(:,:),asuq(:,:)
! magnetic field
real(8), allocatable :: bsu(:,:),bsuq(:,:)

end module

