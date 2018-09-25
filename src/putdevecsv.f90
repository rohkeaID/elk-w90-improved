
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine putdevecsv(ik,devecsv)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: ik
complex(8), intent(in) :: devecsv(nstsv,nstsv)
! local variables
integer recl
character(256) fext
! find the record length
inquire(iolength=recl) vkl(:,ik),nstsv,devecsv
! construct the phonon file extension
call phfext(iqph,isph,iaph,ipph,fext)
!$OMP CRITICAL(u226)
open(226,file=trim(scrpath)//'DEVECSV'//trim(fext),form='UNFORMATTED', &
 access='DIRECT',recl=recl)
write(226,rec=ik) vkl(:,ik),nstsv,devecsv
close(226)
!$OMP END CRITICAL(u226)
return
end subroutine

