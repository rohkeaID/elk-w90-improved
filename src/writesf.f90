
! Copyright (C) 2009 J. K. Dewhurst, S. Sharma and E. K. U. Gross
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine writesf
use modmain
implicit none
! local variables
integer iw
! width of plotting interval in units of swidth
real(8), parameter :: swf=10.d0
real(8) dw,w,x
! external functions
real(8) sdelta,stheta
external sdelta,stheta
open(50,file='SDELTA.OUT',action='WRITE',form='FORMATTED')
open(51,file='STHETA.OUT',action='WRITE',form='FORMATTED')
dw=(2.d0*swf*swidth)/dble(nwplot-1)
do iw=1,nwplot
  w=-swf*swidth+dw*dble(iw-1)
  x=w/swidth
  write(50,'(2G18.10)') w,sdelta(stype,x)/swidth
  write(51,'(2G18.10)') w,stheta(stype,x)
end do
close(50)
close(51)
write(*,*)
write(*,'("Info(writesf): smooth Dirac delta and Heaviside functions written")')
write(*,'(" SDELTA.OUT and STHETA.OUT, respectively")')
return
end subroutine

