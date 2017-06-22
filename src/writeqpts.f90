
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writeqpts
use modmain
implicit none
! local variables
integer iq
open(50,file='QPOINTS.OUT',action='WRITE',form='FORMATTED')
write(50,'(I6," : nqpt; q-point, vql, wqpt below")') nqpt
do iq=1,nqpt
  write(50,'(I6,4G18.10)') iq,vql(:,iq),wqpt(iq)
end do
close(50)
return
end subroutine

