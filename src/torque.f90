
! Copyright (C) 2017 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine torque
use modmain
implicit none
! local variables
integer idm
real(8) torq(3)
! allocatable arrays
real(8), allocatable :: rvfmt(:,:,:),rvfir(:,:)
! external functions
real(8) rfint
external rfint
! initialise universal variables
call init0
if (.not.ncmag) then
  torq(:)=0.d0
  goto 10
end if
! read magnetisation and exchange-correlation magnetic field from file
call readstate
! compute m(r) x B_xc(r)
allocate(rvfmt(npmtmax,natmtot,3),rvfir(ngtot,3))
call rvfcross(magmt,magir,bxcmt,bxcir,rvfmt,rvfir)
! integrate to find the total torque
do idm=1,ndmag
  torq(idm)=rfint(rvfmt(:,:,idm),rvfir(:,idm))
end do
10 continue
write(*,*)
write(*,'("Info(torque):")')
write(*,'(" Total torque exerted by B_xc on the magnetisation :")')
write(*,'(3G18.10)') torq
return
end subroutine

