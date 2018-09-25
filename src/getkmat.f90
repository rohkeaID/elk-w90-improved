
! Copyright (C) 2007-2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine getkmat(ik,kmat)
use modmain
implicit none
! arguments
integer, intent(in) :: ik
complex(8), intent(out) :: kmat(nstsv,nstsv)
! local variables
integer recl,nstsv_,i
real(8) vkl_(3),t1
! find the record length
inquire(iolength=recl) vkl_,nstsv_,kmat
!$OMP CRITICAL(u140)
do i=1,2
  open(140,file='KMAT.OUT',form='UNFORMATTED',access='DIRECT',recl=recl,err=10)
  read(140,rec=ik,err=10) vkl_,nstsv_,kmat
  exit
10 continue
  if (i.eq.2) then
    write(*,*)
    write(*,'("Error(getkmat): unable to read from KMAT.OUT")')
    write(*,*)
    stop
  end if
  close(140)
end do
!$OMP END CRITICAL(u140)
t1=abs(vkl(1,ik)-vkl_(1))+abs(vkl(2,ik)-vkl_(2))+abs(vkl(3,ik)-vkl_(3))
if (t1.gt.epslat) then
  write(*,*)
  write(*,'("Error(getkmat): differing vectors for k-point ",I8)') ik
  write(*,'(" current  : ",3G18.10)') vkl(:,ik)
  write(*,'(" KMAT.OUT : ",3G18.10)') vkl_
  write(*,*)
  stop
end if
if (nstsv.ne.nstsv_) then
  write(*,*)
  write(*,'("Error(getkmat): differing nstsv for k-point ",I8)') ik
  write(*,'(" current  : ",I8)') nstsv
  write(*,'(" KMAT.OUT : ",I8)') nstsv_
  write(*,*)
  stop
end if
return
end subroutine

