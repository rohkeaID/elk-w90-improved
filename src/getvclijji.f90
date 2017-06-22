
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: getvclijji
! !INTERFACE:
subroutine getvclijji(ikp,vclijji)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ikp     : k-point from non-reduced set (in,integer)
!   vclijji : Coulomb matrix elements (out,real(nstsv,nstsv,nkpt))
! !DESCRIPTION:
!   Retrieves Coulomb matrix elements of the type $(i-jj-i)$ from the file
!   {\tt VCLIJJI.OUT}.
!
! !REVISION HISTORY:
!   Created 2009 (Sharma)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ikp
real(8), intent(out) :: vclijji(nstsv,nstsv,nkpt)
! local variables
integer recl,iostat
! determine record length
inquire(iolength=recl) vclijji
!$OMP CRITICAL
open(95,file='VCLIJJI.OUT',action='READ',form='UNFORMATTED',access='DIRECT', &
 recl=recl,iostat=iostat)
if (iostat.ne.0) then
  write(*,*)
  write(*,'("Error(getvclijji): error opening file VCLIJJI.OUT")')
  write(*,*)
  stop
end if
read(95,rec=ikp) vclijji
close(95)
!$OMP END CRITICAL
return
end subroutine
!EOC

