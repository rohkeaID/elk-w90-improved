
! Copyright (C) 2002-2012 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: findngkmax
! !INTERFACE:
subroutine findngkmax(nkpt,vkc,nspnfv,vqcss,ngvec,vgc,gkmax,ngkmax)
! !INPUT/OUTPUT PARAMETERS:
!   nkpt   : number of k-points (in,integer)
!   vkc    : k-point vectors in Cartesian coordinates (in,real(3,nkpt))
!   nspnfv : number of first-variational spin components: 1 normal case, 2 for
!            spin-spiral case (in,integer)
!   vqcss  : spin-spiral q-vector, not referenced if nspnfv=1 (in,integer)
!   ngvec  : number of G-vectors (in,integer)
!   vgc    : G-vectors in Cartesian coordinates (in,real(3,ngvec))
!   gkmax  : maximum allowed |G+k| (in,real)
!   ngkmax : maximum number of G+k-vectors over all k-points (out,integer)
! !DESCRIPTION:
!   Determines the largest number of ${\bf G+k}$-vectors with length less than
!   {\tt gkmax} over all the $k$-points. This variable is used for allocating
!   arrays.
!
! !REVISION HISTORY:
!   Created October 2004 (JKD)
!   Modified, August 2012 (JKD)
!   Removed modmain and added arguments, September 2012 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: nkpt
real(8), intent(in) :: vkc(3,nkpt)
integer, intent(in) :: nspnfv
real(8), intent(in) :: vqcss(3)
integer, intent(in) :: ngvec
real(8), intent(in) :: vgc(3,ngvec)
real(8), intent(in) :: gkmax
integer, intent(out) :: ngkmax
! local variables
integer ispn,ik,n,ig
real(8) v1(3),v2(3),t0,t1
t0=gkmax**2
ngkmax=0
do ispn=1,nspnfv
  do ik=1,nkpt
    if (nspnfv.eq.2) then
! spin-spiral case
      if (ispn.eq.1) then
        v1(:)=vkc(:,ik)+0.5d0*vqcss(:)
      else
        v1(:)=vkc(:,ik)-0.5d0*vqcss(:)
      end if
    else
      v1(:)=vkc(:,ik)
    end if
    n=0
    do ig=1,ngvec
      v2(:)=vgc(:,ig)+v1(:)
      t1=v2(1)**2+v2(2)**2+v2(3)**2
      if (t1.lt.t0) n=n+1
    end do
    ngkmax=max(ngkmax,n)
  end do
end do
return
end subroutine
!EOC
