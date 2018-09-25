
! Copyright (C) 2017 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gengclg
use modmain
implicit none
! local variables
if (allocated(gclg)) deallocate(gclg)
allocate(gclg(ngvec))
gclg(1)=0.d0
gclg(2:ngvec)=fourpi/gc(2:ngvec)**2
return
end subroutine

