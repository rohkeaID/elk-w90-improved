
! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

real(8) function fmin_nm(id,rd,x)
implicit none
! arguments
integer, intent(in) :: id(*)
real(8), intent(in) :: rd(*),x(*)
! external functions
real(8) polefit
external polefit
select case(id(1))
case(1)
  fmin_nm=polefit(rd,x)
case default
  write(*,*)
  write(*,'("Error(fmin_nm): function type not defined : ",I8)') id(1)
  write(*,*)
  stop
end select
return
end function

