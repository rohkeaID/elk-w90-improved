
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine sstask(fnum,fext)
use modmain
implicit none
! arguments
integer, intent(in) :: fnum
character(*), intent(out) :: fext
! local variables
logical exist
do iqss=1,nqpt
! construct the spin-spiral file extension
  call ssfext(iqss,fext)
! determine if the SS file exists
  inquire(file='SS'//trim(fext),exist=exist)
  if (.not.exist) then
    open(fnum,file='SS'//trim(fext),action='WRITE',form='FORMATTED')
    return
  end if
end do
fext='.OUT'
iqss=0
write(*,*)
write(*,'("Info(sstask): nothing more to do")')
return
end subroutine

