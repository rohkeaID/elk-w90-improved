
! Copyright (C) 2017 A. Davydov, A. Sanna, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine getgwsefm(ik,se)
use modmain
use modgw
implicit none
! arguments
integer, intent(in) :: ik
complex(8), intent(out) :: se(nstsv,nstsv,0:nwfm)
! local variables
integer recl,nstsv_,nwfm_
real(8) vkl_(3),t1
! find the record length
inquire(iolength=recl) vkl_,nstsv_,nwfm_,se
!$OMP CRITICAL(u280)
open(280,file='GWSEFM.OUT',form='UNFORMATTED',access='DIRECT',recl=recl)
read(280,rec=ik) vkl_,nstsv_,nwfm_,se
close(280)
!$OMP END CRITICAL(u280)
t1=abs(vkl(1,ik)-vkl_(1))+abs(vkl(2,ik)-vkl_(2))+abs(vkl(3,ik)-vkl_(3))
if (t1.gt.epslat) then
  write(*,*)
  write(*,'("Error(getgwsefm): differing vectors for k-point ",I8)') ik
  write(*,'(" current    : ",3G18.10)') vkl(:,ik)
  write(*,'(" GWSEFM.OUT : ",3G18.10)') vkl_
  write(*,*)
  stop
end if
if (nstsv.ne.nstsv_) then
  write(*,*)
  write(*,'("Error(getgwsefm): differing nstsv for k-point ",I8)') ik
  write(*,'(" current    : ",I8)') nstsv
  write(*,'(" GWSEFM.OUT : ",I8)') nstsv_
  write(*,*)
  stop
end if
if (nwfm.ne.nwfm_) then
  write(*,*)
  write(*,'("Error(getgwsefm): differing nwfm for k-point ",I8)') ik
  write(*,'(" current    : ",I8)') nwfm
  write(*,'(" GWSEFM.OUT : ",I8)') nwfm_
  write(*,*)
  stop
end if
return
end subroutine

