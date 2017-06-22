
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rdmwriteengy
! !INTERFACE:
subroutine rdmwriteengy(fnum)
! !USES:
use modmain
use modrdm
! !INPUT/OUTPUT PARAMETERS:
!   fnum : file number for writing output (in,integer)
! !DESCRIPTION:
!   Writes all contributions to the total energy to file.
!
! !REVISION HISTORY:
!   Created 2008 (Sharma)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: fnum
write(fnum,*)
write(fnum,'("Energies :")')
write(fnum,'(" electronic kinetic",T30,": ",G18.10)') engykn
write(fnum,'(" core electron kinetic",T30,": ",G18.10)') engykncr
write(fnum,'(" Coulomb",T30,": ",G18.10)') engyvcl
write(fnum,'(" Madelung",T30,": ",G18.10)') engymad
write(fnum,'(" exchange-correlation",T30,": ",G18.10)') engyx
if (rdmtemp.gt.0.d0) then
  write(fnum,'(" entropy",T30,": ",G18.10)') rdmentrpy
end if
write(fnum,'(" total energy",T30,": ",G18.10)') engytot
call flushifc(fnum)
return
end subroutine
!EOC
