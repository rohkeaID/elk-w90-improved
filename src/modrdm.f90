
! Copyright (C) 2008 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modrdm
use modmain

!--------------------------------------------------------------------!
!     reduced density matrix functional theory (RDMFT) variables     !
!--------------------------------------------------------------------!
! Coulomb potential matrix elements
complex(8), allocatable :: vclmat(:,:,:)
! derivative of kinetic energy w.r.t. natural orbital coefficients
complex(8), allocatable :: dkdc(:,:,:)
! step size for occupation numbers
real(8) taurdmn
! step size for natural orbital coefficients
real(8) taurdmc
! xc functional
integer rdmxctype
! maximum number of self-consistent loops
integer rdmmaxscl
! maximum number of iterations for occupation number optimisation
integer maxitn
! maximum number of iteration for natural orbital optimisation
integer maxitc
! exponent for the power and hybrid functional
real(8) rdmalpha
! mixing for the hybrid functional
real(8) rdmbeta
! temperature
real(8) rdmtemp
! entropy
real(8) rdmentrpy

end module
