
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
integer recl,i
character(256) fname
! find the record length
inquire(iolength=recl) vkl(:,ik),nstsv,evecsv
fname=trim(scrpath)//'EVECSV'//trim(fext)
!$OMP CRITICAL(u126)
do i=1,2
  open(126,file=trim(fname),form='UNFORMATTED',access='DIRECT',recl=recl,err=10)
  write(126,rec=ik,err=10) vkl(:,ik),nstsv,evecsv
  close(126)
  exit
10 continue
  if (i.eq.2) then
    write(*,*)
    write(*,'("Error(putevecsv): unable to write to ",A)') trim(fname)
    write(*,*)
    stop
  end if
  close(126)
end do
!$OMP END CRITICAL(u126)
return
end subroutine

