
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modtddft

!-----------------------------------------!
!     TDDFT linear response variables     !
!-----------------------------------------!
! exchange-correlation kernel type
integer fxctype(3)
! parameters for long-range correction (LRC) kernel
real(8) fxclrc(2)
! number of independent spin components of the f_xc spin tensor
integer nscfxc

!---------------------------------------------!
!     TDDFT real-time evolution variables     !
!---------------------------------------------!
! number of laser pulses defining the time-dependent A-field
integer npulse
! laser pulse parameters: vector amplitude, peak time, full-width at
! half-maximum, frequency and phase
real(8), allocatable :: pulse(:,:)
! number of A-field ramps
integer nramp
! ramp parameters: vector amplitude, ramp start time, linear and quadratic
! coefficients
real(8), allocatable :: ramp(:,:)
! total simulation time
real(8) tstime
! time step length
real(8) dtimes
! number of time-steps
integer ntimes
! current time-step
integer itimes
! time steps
real(8), allocatable :: times(:)
! time-dependent A-field
real(8), allocatable :: afieldt(:,:)
! number of time steps after which observables are written to file
integer ntswrite
! the following variables are .true. if the corresponding quantities are to be
! written every ntswrite time steps
logical tdrho1d,tdrho2d,tdrho3d
logical tdmag2d,tdmag3d

end module

