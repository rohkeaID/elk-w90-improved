
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: symveca
! !INTERFACE:
subroutine symveca(vca)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   vca   : vectors in Cartesian coordinates for all atoms (in,real(3,natmtot))
! !DESCRIPTION:
!   Symmetrises a 3-vector at each atomic site by rotating and averaging over
!   equivalent atoms. Depending on {\tt tspin}, either the spatial or spin part
!   of the crystal symmetries are used.
!
! !REVISION HISTORY:
!   Created June 2004 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(inout) :: vca(3,natmtot)
! local variables
integer is,ia,ja,ias,jas
integer isym,ls
real(8) v(3),t1
! automatic arrays
real(8) vs(3,natmtot)
! make symmetric average
vs(:,:)=0.d0
do isym=1,nsymcrys
  ls=lsplsymc(isym)
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      ja=ieqatom(ia,is,isym)
      jas=idxas(ja,is)
      call r3mv(symlatc(:,:,ls),vca(:,jas),v)
      vs(:,ias)=vs(:,ias)+v(:)
    end do
  end do
end do
! normalise
t1=1.d0/dble(nsymcrys)
vca(:,:)=t1*vs(:,:)
return
end subroutine
!EOC

