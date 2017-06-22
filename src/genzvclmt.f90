
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genzvclmt(nr,nri,ld1,r,ld2,zrhomt,zvclmt)
use modmain
implicit none
! arguments
integer, intent(in) :: nr(nspecies),nri(nspecies)
integer, intent(in) :: ld1
real(8), intent(in) :: r(ld1,nspecies)
integer, intent(in) :: ld2
complex(8), intent(in) :: zrhomt(lmmaxvr,ld2,natmtot)
complex(8), intent(out) :: zvclmt(lmmaxvr,ld2,natmtot)
! local variables
integer is,ias
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(is)
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  call zpotclmt(nr(is),nri(is),r(:,is),zrhomt(:,:,ias),zvclmt(:,:,ias))
end do
!$OMP END DO
!$OMP END PARALLEL
return
end subroutine

