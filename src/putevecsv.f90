
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine putevecsv(fext,ik,evecsv)
use modmain
implicit none
! arguments
character(*), intent(in) :: fext
integer, intent(in) :: ik
complex(8), intent(in) :: evecsv(nstsv,nstsv)
! local variables
integer recl
! find the record length
inquire(iolength=recl) vkl(:,ik),nstsv,evecsv
!$OMP CRITICAL
open(70,file=trim(scrpath)//'EVECSV'//trim(fext),action='WRITE', &
 form='UNFORMATTED',access='DIRECT',recl=recl)
write(70,rec=ik) vkl(:,ik),nstsv,evecsv
close(70)
!$OMP END CRITICAL
return
end subroutine


