
! Copyright (C) 2017 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine acgwse(ist,jst,swfm,wr,swr)
use modmain
use modgw
implicit none
! arguments
integer ist,jst
complex(8), intent(in) :: swfm(nstsv,nstsv,0:nwfm)
real(8), intent(in) :: wr(nwplot)
complex(8), intent(out) :: swr(nstsv,nstsv,nwplot)
! allocatable arrays
complex(8), allocatable :: zfm(:),zwr(:),zr(:)
allocate(zfm(0:nwfm),zwr(nwplot),zr(nwplot))
zfm(:)=swfm(ist,jst,:)
zwr(:)=wr(:)
select case(actype)
case(1)
! fit a multipole model
  call acpole(zfm,zwr,zr)
case(10)
! stabilised Pade approximant
  call pades(nspade,swidth,nwfm+1,wfm,zfm,nwplot,zwr,zr)
case default
  write(*,*)
  write(*,'("Error(acgwse): actype not defined : ",I8)') actype
  write(*,*)
  stop
end select
swr(ist,jst,:)=zr(:)
deallocate(zfm,zwr,zr)
return
end subroutine

