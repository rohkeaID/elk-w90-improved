
! Copyright (C) 2011 A. Sanna and E. K. U. Gross
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine readalpha2f(w,a2f)
use modmain
implicit none
! arguments
real(8), intent(out) :: w(nwplot)
real(8), intent(out) :: a2f(nwplot)
! local variables
integer iw,iostat
open(50,file='ALPHA2F.OUT',action='READ',form='FORMATTED',status='OLD', &
 iostat=iostat)
if (iostat.ne.0) then
  write(*,*)
  write(*,'("Error(readalpha2f): error opening ALPHA2F.OUT")')
  write(*,*)
  stop
end if
do iw=1,nwplot
  read(50,*,iostat=iostat) w(iw),a2f(iw)
  if (iostat.ne.0) then
    write(*,*)
    write(*,'("Error(readalpha2f): error reading from ALPHA2F.OUT")')
    write(*,'(" for frequency ",I6)') iw
    write(*,*)
    stop
  end if
end do
close(50)
return
end subroutine

