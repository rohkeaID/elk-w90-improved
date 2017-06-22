
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine getdevalfv(ik,iq,is,ia,ip,devalfv)
use modmain
implicit none
! arguments
integer, intent(in) :: ik,iq,is,ia,ip
complex(8), intent(out) :: devalfv(nstfv,nspnfv)
! local variables
integer recl,nstfv_,nspnfv_
real(8) vkl_(3),t1
character(256) fext,fname
! construct the phonon file extension
call phfext(iq,is,ia,ip,fext)
! construct filename
fname='DEVALFV'//trim(fext)
! find the record length
inquire(iolength=recl) vkl_,nstfv_,nspnfv_,devalfv
!$OMP CRITICAL
open(70,file=trim(fname),action='READ',form='UNFORMATTED',access='DIRECT', &
 recl=recl)
read(70,rec=ik) vkl_,nstfv_,nspnfv_,devalfv
close(70)
!$OMP END CRITICAL
t1=abs(vkl(1,ik)-vkl_(1))+abs(vkl(2,ik)-vkl_(2))+abs(vkl(3,ik)-vkl_(3))
if (t1.gt.epslat) then
  write(*,*)
  write(*,'("Error(getdevalfv): differing vectors for k-point ",I8)') ik
  write(*,'(" current : ",3G18.10)') vkl(:,ik)
  write(*,'(" ",A," : ",3G18.10)') trim(fname),vkl_
  write(*,*)
  stop
end if
if (nstfv.ne.nstfv_) then
  write(*,*)
  write(*,'("Error(getdevalfv): differing nstfv for k-point ",I8)') ik
  write(*,'(" current : ",I8)') nstfv
  write(*,'(" ",A," : ",I8)') trim(fname),nstfv_
  write(*,*)
  stop
end if
if (nspnfv.ne.nspnfv_) then
  write(*,*)
  write(*,'("Error(getdevalfv): differing nspnfv for k-point ",I8)') ik
  write(*,'(" current : ",I8)') nspnfv
  write(*,'(" ",A," : ",I8)') trim(fname),nspnfv_
  write(*,*)
  stop
end if
return
end subroutine

