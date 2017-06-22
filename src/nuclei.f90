
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine nuclei
use modmain
implicit none
! local variables
integer is,ir,irc
! external functions
real(8) radnucl
external radnucl
do is=1,nspecies
! approximate nuclear radius
  rnucl(is)=radnucl(spzn(is))
! nuclear volume
  vnucl(is)=(4.d0/3.d0)*pi*rnucl(is)**3
! number of radial mesh points to nuclear radius
  nrnucl(is)=1
  do ir=1,nrmt(is)
    if (rsp(ir,is).gt.rnucl(is)) then
      nrnucl(is)=ir
      exit
    end if
  end do
! number of coarse radial mesh points to nuclear radius
  nrcnucl(is)=1
  do irc=1,nrcmt(is)
    if (rcmt(irc,is).gt.rnucl(is)) then
      nrcnucl(is)=irc
      exit
    end if
  end do
end do
return
end subroutine

