
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine putdevecfv(ik,devecfv)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: ik
complex(8), intent(in) :: devecfv(nmatmax,nstfv,nspnfv)
! local variables
integer recl
character(256) fext
! find the record length
inquire(iolength=recl) vkl(:,ik),nmatmax,nstfv,nspnfv,devecfv
! construct the phonon file extension
call phfext(iqph,isph,iaph,ipph,fext)
!$OMP CRITICAL
open(70,file=trim(scrpath)//'DEVECFV'//trim(fext),action='WRITE', &
 form='UNFORMATTED',access='DIRECT',recl=recl)
write(70,rec=ik) vkl(:,ik),nmatmax,nstfv,nspnfv,devecfv
close(70)
!$OMP END CRITICAL
return
end subroutine

