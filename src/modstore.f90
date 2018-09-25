
! Copyright (C) 2012 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!---------------------------------------------------------!
!     variables for storing original input parameters     !
!---------------------------------------------------------!

module modstore
use modmain

real(8) avec_(3,3)
real(8) bvec_(3,3),binv_(3,3)
real(8) omega_
logical tshift_
logical primcell_
integer natoms_(maxspecies)
integer natmtot_
integer idxis_(maxatoms*maxspecies)
real(8) atposl_(3,maxatoms,maxspecies)
real(8) atposc_(3,maxatoms,maxspecies)
integer ngridg_(3),ngtot_
integer, allocatable :: ivg_(:,:),igfft_(:)
logical spinpol_,spinorb_,cmagz_,spinsprl_
real(8) bfieldc0_(3)
real(8) bfcmt0_(3,maxatoms,maxspecies)
real(8) reducebf_
integer fsmtype_
real(8) momfix_(3)
real(8) mommtfix_(3,maxatoms,maxspecies)
logical tforce_
logical autokpt_
integer ngridk_(3)
real(8) vkloff_(3)
integer lmaxi_
logical tfav0_
real(8) vqlss_(3)
integer msmooth_

end module

