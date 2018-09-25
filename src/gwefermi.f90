
! Copyright (C) 2018 P. Elliott, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gwefermi
use modmain
use modgw
implicit none

! local variables
integer, parameter :: maxit=1000
integer ik,it
integer ist,iw
real(8) e
!****** clean up

do it=1,maxit
  do ik=1,nkpt
! loop over the fermionic Matsubara frequencies
    do iw=-nwfm,nwfm,2
! compute the diagonal matrix G_s
      do ist=1,nstsv
        e=evalsv(ist,ik)-efermi

      end do
    end do
  end do
end do



return
end subroutine

