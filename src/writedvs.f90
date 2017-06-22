
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writedvs(fext)
use modmain
use modphonon
implicit none
! arguments
character(*), intent(in) :: fext
! local variables
integer is
open(50,file='DVS'//trim(fext),action='WRITE',form='UNFORMATTED')
write(50) version
write(50) nspecies
write(50) lmmaxvr
do is=1,nspecies
  write(50) natoms(is)
  write(50) nrmt(is)
end do
write(50) ngridg
write(50) dvsmt,dvsir
close(50)
return
end subroutine

