
! Copyright (C) 2015 Jon Lafuente and Manh Duc Le
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modw90

!-------------------------------------------------!
!     Wannier90 interface variables               !
!-------------------------------------------------!
! conversion factor from a.u. (Bohr) to Angstrom (needed 4: libwannier.a)
real(8),      parameter   :: bohr2angstrom = 0.52917721067
real(8),      parameter   :: hartree2ev    = 27.2113850560
! total number of atoms to pass to wannier90
integer                      wann_natoms
! number of bands to pass to wannier90
integer                      wann_nband
! atom symbols to pass to wannier90
character(4), allocatable :: wann_atomsymb(:)
! atom Cartesian coordinates in Angstrom
real(8),      allocatable :: wann_atompos(:,:)
! index of bands to pass to wannier90
integer,      allocatable :: wann_bands(:)
! list of neighbouring k-points
integer,      allocatable :: wann_nnkp(:,:)
! number of iterations for the minimisation in Wannier90 iterations
integer                      wann_numiter

! wannier projections definitions. Uses same syntax as wannier90 library
! number of lines in projection block
integer                      wann_projlines
! number of lines in wann_input block
integer                      wann_inputlines
! the projection block to write to .win
character(256)               wann_projstr(256)
! arbitrary input for wannier90 .win
character(256)               wann_input(256)
! name of all files, produced by interface; default: ELK
character(256)               wann_seedname
! k-point grid sizes for Wannier calculations
integer                      wann_ngridk(3)
! number of projections per spin
integer                      wann_nproj
! Ang/Bohr units (cart coords)
character(4)                 wann_proj_units
! if proj has weight in this MT
logical,      allocatable :: wann_proj_haswt(:,:)
! use a random s-type gaussian
logical,      allocatable :: wann_proj_isrand(:)
! the radial part on coarse mesh
real(8),      allocatable :: wann_projulr(:,:,:,:)
! the coefficients in Ylm basis
complex(8),   allocatable :: wann_projclm(:,:,:)
! maximum number of nearest neighbours per k-point
integer                   :: num_nnmax = 12

! total number of nearest neighbours for each k-point
integer                      wann_nntot
! list of nearest neighbours for each k-point
integer,      allocatable :: nnlist(:,:)
! (see the notation of Wannier90)
integer,      allocatable :: nncell(:,:,:)
! number of bands in first-principles calc. (excluding eg. semi-core states)
integer                   :: wann_nband_total
! number of wannier functions to calculate
integer                      wann_nwf
! coordinates of sites
real(8),      allocatable :: wann_proj_site(:,:)
! l angular momentum number
integer,      allocatable :: wann_proj_l(:)
! m azimuth number
integer,      allocatable :: wann_proj_m(:)
! type of radial function
integer,      allocatable :: wann_proj_radial(:)
! local z-axis
real(8),      allocatable :: wann_proj_zaxis(:,:)
! local x-axis
real(8),      allocatable :: wann_proj_xaxis(:,:)
! Z/a of radial function
real(8),      allocatable :: wann_proj_zona(:)
! k-points independant list of bands to exclude from the calculation
integer,      allocatable :: wann_proj_exclude_bands_lib(:)
! u(+1) or d(-1) in quantdir
integer,      allocatable :: wann_proj_spin(:)
! local spin quantisation dir
real(8),      allocatable :: wann_proj_quantdir(:,:)

end module
