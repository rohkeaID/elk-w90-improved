
! Copyright (C) 2015 K. Krieger, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine readtimes(itimes0)
use modmain
use modtddft
implicit none
! arguments
integer, intent(out) :: itimes0
! local variables
integer iostat
real(8) times_,t1
open(50,file='TIMESTEP.OUT',action='READ',form='FORMATTED',status='OLD', &
 iostat=iostat)
if (iostat.ne.0) then
  write(*,*)
  write(*,'("Error(readtimes): error opening TIMESTEP.OUT")')
  write(*,*)
  stop
end if
read(50,*,iostat=iostat) itimes0,times_
if (iostat.ne.0) then
  write(*,*)
  write(*,'("Error(readtimes): error reading time step from TIMESTEP.OUT")')
  write(*,*)
  stop
end if
if ((itimes0.lt.1).or.(itimes0.gt.ntimes)) then
  write(*,*)
  write(*,'("Error(readtimes): invalid itimes : ",I8)') itimes0
  write(*,*)
  stop
end if
t1=abs(times(itimes0)-times_)
if (t1.gt.1.d-8) then
  write(*,*)
  write(*,'("Error(readtimes): differing time step")')
  write(*,'(" current      : ",G18.10)') times(itimes0)
  write(*,'(" TIMESTEP.OUT : ",G18.10)') times_
  write(*,*)
  stop
end if
return
end subroutine

