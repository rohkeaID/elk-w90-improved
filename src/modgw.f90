
! Copyright (C) 2016 A. Davydov, A. Sanna, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

module modgw

! maximum Matsubara frequency for the GW calculation
real(8) wmaxgw
! maximum number of Matsubara frequencies
integer nwgw
! integer grid intervals for Matsubara frequencies
integer intwgw(2)
! map from frequency index to FFT array
integer, allocatable :: iwfft(:)
! maximum fermionic Matsubara frequency index to be used for the GW calculation
integer nwfm
! maximum bosonic frequency index
integer nwbs
! imaginary frequencies used for the GW calculation
real(8), allocatable :: wgw(:)
! diagonal approximations for the GW self-energy
!  0 : Sigma and W_c are full matrices
!  1 : Sigma is taken to be diagonal
!  2 : Sigma and W_c are taken to be diagonal
integer gwdiag

end module

