
! Copyright (C) 2017 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine readcurrent(jt)
use modmain
use modtddft
implicit none
! arguments
real(8), intent(out) :: jt(3,ntimes)
! local variables
integer its,iostat
real(8) times_,t1
! initialise global variables
call init0
call init1
open(50,file='CURRENT_TD.OUT',form='FORMATTED',status='OLD',iostat=iostat)
if (iostat.ne.0) then
  write(*,*)
  write(*,'("Error(readcurrent): error opening CURRENT_TD.OUT")')
  write(*,*)
  stop
end if
do its=1,ntimes-1
  read(50,*) times_,jt(:,its)
  t1=abs(times(its)-times_)
  if (t1.gt.1.d-10) then
    write(*,*)
    write(*,'("Error(readcurrent): time step mismatch : ",G18.10)')
    write(*,'(" internal       : ",G18.10)') times(its)
    write(*,'(" CURRENT_TD.OUT : ",G18.10)') times_
    write(*,*)
    stop
  end if
end do
close(50)
! the final current is not computed so set it equal to the penultimate
jt(:,ntimes)=jt(:,ntimes-1)
return
end subroutine

