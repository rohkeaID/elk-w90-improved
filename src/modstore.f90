
! Copyright (C) 2012 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!---------------------------------------------------------!
!     variables for storing original input parameters     !
!---------------------------------------------------------!

module modstore
use modmain

real(8) avec0(3,3)
real(8) bvec0(3,3),binv0(3,3)
real(8) omega0
logical tshift0
logical primcell0
integer natoms0(maxspecies)
integer natmtot0
real(8) atposl0(3,maxatoms,maxspecies)
real(8) atposc0(3,maxatoms,maxspecies)
real(8) rmtdelta0
integer ngridg0(3),ngtot0
integer, allocatable :: ivg0(:,:),igfft0(:)
logical spinpol0,spinorb0,cmagz0
real(8) bfieldc00(3)
real(8) bfcmt00(3,maxatoms,maxspecies)
real(8) reducebf0
integer fsmtype0
real(8) momfix0(3)
real(8) mommtfix0(3,maxatoms,maxspecies)
logical tforce0
logical autokpt0
integer ngridk0(3)
real(8) vkloff0(3)
integer lmaxinr0
logical trimvg0

end module

