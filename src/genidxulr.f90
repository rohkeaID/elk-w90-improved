
! Copyright (C) 2017 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genidxulr

!***** check

use modmain
use modulr
implicit none
! local variables
integer ik0,ik,ist,i
! find maximum number of second-variational states for ultra long-range
nsvukmax=0
do ik0=1,nkpt0
! central k-point
  ik=(ik0-1)*nkpa+1
  i=0
  do ist=1,nstsv
    if (abs(evalsv(ist,ik)-efermi).lt.emaxulr) i=i+1
  end do
  if (i.gt.nsvukmax) nsvukmax=i
end do
if (nsvukmax.eq.0) then
  write(*,*)
  write(*,'("Error(genidxulr): no second-variational states available for &
   &ultra long-range calculation")')
  write(*,'("Increase emaxulr")')
  write(*,*)
  stop
end if
if (allocated(nsvuk)) deallocate(nsvuk)
allocate(nsvuk(nkpt0))
if (allocated(istuk)) deallocate(istuk)
allocate(istuk(nsvukmax,nkpt0))
! find number of and index to second-variational states for ultra long-range
do ik0=1,nkpt0
  ik=(ik0-1)*nkpa+1
  i=0
  do ist=1,nstsv
    if (abs(evalsv(ist,ik)-efermi).lt.emaxulr) then
      i=i+1
      istuk(i,ik0)=ist
    end if
  end do
  nsvuk(ik0)=i
end do
return
end subroutine

