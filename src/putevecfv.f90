
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine putevecfv(fext,ik,evecfv)
use modmain
implicit none
! arguments
character(*), intent(in) :: fext
integer, intent(in) :: ik
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv)
! local variables
integer recl,i
character(256) fname
! find the record length
inquire(iolength=recl) vkl(:,ik),nmatmax,nstfv,nspnfv,evecfv
fname=trim(scrpath)//'EVECFV'//trim(fext)
!$OMP CRITICAL(u122)
do i=1,2
  open(122,file=trim(fname),form='UNFORMATTED',access='DIRECT',recl=recl,err=10)
  write(122,rec=ik,err=10) vkl(:,ik),nmatmax,nstfv,nspnfv,evecfv
  close(122)
  exit
10 continue
  if (i.eq.2) then
    write(*,*)
    write(*,'("Error(putevecfv): unable to write to ",A)') trim(fname)
    write(*,*)
    stop
  end if
  close(122)
end do
!$OMP END CRITICAL(u122)
return
end subroutine

