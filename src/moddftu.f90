
! Copyright (C) 2008 F. Bultmark, F. Cricchio and L. Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module moddftu
use modmain

!-----------------------------------------------------------!
!     muffin-tin density and potential matrix variables     !
!-----------------------------------------------------------!
! maximum angular momentum for muffin-tin density matrix
integer, parameter :: lmaxdm=3
integer, parameter :: lmmaxdm=(lmaxdm+1)**2
! density matrix in each muffin-tin
complex(8), allocatable :: dmatmt(:,:,:,:,:)
! potential matrix in each muffin-tin
complex(8), allocatable :: vmatmt(:,:,:,:,:)
! tvmatmt is .true. if the potential matrices are calculated
logical tvmatmt
! tvmmt is .true. if the potential matrix for that l and atom is non-zero
logical, allocatable :: tvmmt(:,:)

!-------------------------!
!     DFT+U variables     !
!-------------------------!
! type of DFT+U to use (0 = none)
integer dftu
! input type for DFT+U calculation (1:5)
integer inpdftu
! maximum number of DFT+U entries
integer, parameter :: maxdftu=40
! number of DFT+U entries
integer ndftu
! species and angular momentum for each entry
integer idftu(2,maxdftu)
! U and J values for each entry
real(8) ujdu(2,maxdftu)
! interpolation constant alpha for each atom and entry [PRB 67, 153106 (2003)]
real(8), allocatable :: alphadu(:,:)
! readadu is .true. if alphadu is to be read from file
logical readadu
! DFT+U energy for each atom and entry
real(8), allocatable :: engyadu(:,:)
! energy from the DFT+U correction
real(8) engydu
! Slater parameters
real(8) fdu(0:2*lmaxdm,maxdftu)
! Racah parameters
real(8) edu(0:lmaxdm,maxdftu)
! screening length of Yukawa potential to calculate Slater integrals
real(8) lambdadu(maxdftu)
! energies to calculate radial functions for Slater integrals
real(8), allocatable :: fdue(:,:)
! radial functions to calculate Slater integrals
real(8), allocatable :: fdufr(:,:,:)
! fixed value of U for which screening length has to be determined
real(8) udufix(maxdftu)
! initial values of screening length if U is fixed
real(8) lambdadu0(maxdftu)

!---------------------------------!
!     tensor moment variables     !
!---------------------------------!
! tmwrite is .true. if tensor moments are written out at every s.c. loop
logical tmwrite
! fixed tensor moment type
!  0      : none
!  2 (-2) : fixed 2-index tensor moment (or just lowering the symmetry)
!  3 (-3) : fixed 3-index tensor moment (or just lowering the symmetry)
integer ftmtype
! number of fixed tensor moment entries
integer ntmfix
! tensor moment indices for each entry: is, ia, l, n and k, p, x, y for the
! 2-index tensor or k, p, r, t for the 3-index tensor
integer, allocatable :: itmfix(:,:)
! tensor component
complex(8), allocatable :: tmfix(:)
! spatial and spin rotation matrices of tensor
real(8), allocatable :: rtmfix(:,:,:,:)
! density matrices corresponding to the fixed tensor moments
complex(8), allocatable :: dmftm(:,:,:,:,:)
! fixed tensor moment potential matrix
complex(8), allocatable :: vmftm(:,:,:,:,:)
! fixed tensor moment step size
real(8) tauftm
! number of self-consistent loops after which FTM field is updated
integer ftmstep

end module

