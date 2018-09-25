
! Copyright (C) 2017 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

real(8) function polefit(zfm,c)
use modgw
implicit none
! arguments
complex(8), intent(in) :: zfm(0:nwfm)
complex(8), intent(in) :: c(2*npole+1)
! local variables
integer iw
real(8) sum
complex(8) z1
! external functions
complex(8) zfpole
external zfpole
sum=0.d0
do iw=0,nwfm
  z1=zfm(iw)-zfpole(c,wfm(iw))
  sum=sum+dble(z1)**2+aimag(z1)**2
end do
polefit=sum
return
end function

