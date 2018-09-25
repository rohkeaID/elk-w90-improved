
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine wxcplot
use modmain
implicit none
! initialise universal variables
call init0
if (xcgrad.ne.4) then
  write(*,*)
  write(*,'("Error(wxcplot): tau-DFT not in use")')
  write(*,*)
  stop
end if
! read the density and potentials from file
call readstate
! write the potential plots to file
select case(task)
case(341)
  open(50,file='WXC1D.OUT',form='FORMATTED')
  open(51,file='WLINES.OUT',form='FORMATTED')
  call plot1d(50,51,1,wxcmt,wxcir)
  close(50)
  close(51)
  write(*,*)
  write(*,'("Info(wxcplot):")')
  write(*,'(" 1D tau-DFT exchange-correlation potential written to WXC1D.OUT")')
  write(*,'(" vertex location lines written to WLINES.OUT")')
case(342)
  open(50,file='WXC2D.OUT',form='FORMATTED')
  call plot2d(50,1,wxcmt,wxcir)
  close(50)
  write(*,*)
  write(*,'("Info(wxcplot):")')
  write(*,'(" 2D tau-DFT exchange-correlation potential written to WXC2D.OUT")')
case(343)
  open(50,file='WXC3D.OUT',form='FORMATTED')
  call plot3d(50,1,wxcmt,wxcir)
  close(50)
  write(*,*)
  write(*,'("Info(wxcplot):")')
  write(*,'(" 3D tau-DFTexchange-correlation potential written to WXC3D.OUT")')
end select
return
end subroutine

