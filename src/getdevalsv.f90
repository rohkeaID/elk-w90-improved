
! Copyright (C) 2015 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine getdevalsv(ik,iq,is,ia,ip,devalsv)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: ik,iq,is,ia,ip
real(8), intent(out) :: devalsv(nstsv)
! local variables
integer recl,nstsv_
real(8) vkl_(3),t1
character(256) fext,fname
! construct the phonon file extension
call phfext(iq,is,ia,ip,fext)
! construct filename
fname='DEVALSV'//trim(fext)
! find the record length
inquire(iolength=recl) vkl_,nstsv_,devalsv
!$OMP CRITICAL
open(70,file=trim(fname),action='READ',form='UNFORMATTED',access='DIRECT', &
 recl=recl)
read(70,rec=ik) vkl_,nstsv_,devalsv
close(70)
!$OMP END CRITICAL
t1=abs(vkl(1,ik)-vkl_(1))+abs(vkl(2,ik)-vkl_(2))+abs(vkl(3,ik)-vkl_(3))
if (t1.gt.epslat) then
  write(*,*)
  write(*,'("Error(getdevalsv): differing vectors for k-point ",I8)') ik
  write(*,'(" current : ",3G18.10)') vkl(:,ik)
  write(*,'(" ",A," : ",3G18.10)') trim(fname),vkl_
  write(*,*)
  stop
end if
if (nstsv.ne.nstsv_) then
  write(*,*)
  write(*,'("Error(getdevalsv): differing nstsv for k-point ",I8)') ik
  write(*,'(" current : ",I8)') nstsv
  write(*,'(" ",A," : ",I8)') trim(fname),nstsv_
  write(*,*)
  stop
end if
return
end subroutine

