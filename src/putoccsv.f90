
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine putoccsv(fext,ik,occsvp)
use modmain
implicit none
! arguments
character(*), intent(in) :: fext
integer, intent(in) :: ik
real(8), intent(in) :: occsvp(nstsv)
! local variables
integer recl
! find the record length
inquire(iolength=recl) vkl(:,ik),nstsv,occsvp
!$OMP CRITICAL(u130)
open(130,file='OCCSV'//trim(fext),form='UNFORMATTED',access='DIRECT',recl=recl)
write(130,rec=ik) vkl(:,ik),nstsv,occsvp
close(130)
!$OMP END CRITICAL(u130)
return
end subroutine

