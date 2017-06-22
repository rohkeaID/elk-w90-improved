
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine putevalfv(fext,ik,evalfv)
use modmain
implicit none
! arguments
character(*), intent(in) :: fext
integer, intent(in) :: ik
real(8), intent(in) :: evalfv(nstfv,nspnfv)
! local variables
integer recl
! find the record length
inquire(iolength=recl) vkl(:,ik),nstfv,nspnfv,evalfv
!$OMP CRITICAL
open(70,file='EVALFV'//trim(fext),action='WRITE',form='UNFORMATTED', &
 access='DIRECT',recl=recl)
write(70,rec=ik) vkl(:,ik),nstfv,nspnfv,evalfv
close(70)
!$OMP END CRITICAL
return
end subroutine

