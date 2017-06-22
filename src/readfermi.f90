
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: readfermi
! !INTERFACE:
subroutine readfermi
! !USES:
use modmain
! !DESCRIPTION:
!   Reads the Fermi energy from the file {\tt EFERMI.OUT}.
!
! !REVISION HISTORY:
!   Created March 2005 (JKD)
!EOP
!BOC
implicit none
! local variables
integer iostat
open(50,file='EFERMI'//trim(filext),action='READ',form='FORMATTED', &
 status='OLD',iostat=iostat)
if (iostat.ne.0) then
  write(*,*)
  write(*,'("Error(readfermi): error opening ",A)') 'EFERMI'//trim(filext)
  write(*,*)
  stop
end if
read(50,*,iostat=iostat) efermi
if (iostat.ne.0) then
  write(*,*)
  write(*,'("Error(readfermi): error reading Fermi energy from EFERMI.OUT")')
  write(*,*)
  stop
end if
close(50)
return
end subroutine
!EOC

