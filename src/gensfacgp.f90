
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: gensfacgp
! !INTERFACE:
subroutine gensfacgp(ngp,vgpc,ld,sfacgp)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ngp    : number of G+p-vectors (in,integer)
!   vgpc   : G+p-vectors in Cartesian coordinates (in,real(3,*))
!   ld     : leading dimension (in,integer)
!   sfacgp : structure factors of G+p-vectors (out,complex(ld,natmtot))
! !DESCRIPTION:
!   Generates the atomic structure factors for a set of ${\bf G+p}$-vectors:
!   $$ S_{\alpha}({\bf G+p})=\exp(i({\bf G+p})\cdot{\bf r}_{\alpha}), $$
!   where ${\bf r}_{\alpha}$ is the position of atom $\alpha$.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ngp
real(8), intent(in) :: vgpc(3,ngp)
integer, intent(in) :: ld
complex(8), intent(out) :: sfacgp(ld,natmtot)
! local variables
integer is,ia,ias,igp
real(8) v(3),t1
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(is,ia,v,igp,t1)
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  ia=idxia(ias)
  v(:)=atposc(:,ia,is)
  do igp=1,ngp
    t1=vgpc(1,igp)*v(1)+vgpc(2,igp)*v(2)+vgpc(3,igp)*v(3)
    sfacgp(igp,ias)=cmplx(cos(t1),sin(t1),8)
  end do
end do
!$OMP END DO
!$OMP END PARALLEL
return
end subroutine
!EOC

