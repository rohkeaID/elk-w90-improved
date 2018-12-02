! Copyright (C) 2015 Jon Lafuente and Manh Duc Le
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: runw90
! !INTERFACE:
subroutine runw90
! !USES:
! !DESCRIPTION:
!
!
! !REVISION HISTORY:
!   Created December 2017 (Arsenii Gerasimov)
!EOP
!BOC
implicit none
!-------------------------------------------------------------------------------
write(*,*)
write(*,'("Error(Wannier90 interface): 603 task is currently unavailable")')
write(*,*)

stop

end subroutine runw90
!EOC
