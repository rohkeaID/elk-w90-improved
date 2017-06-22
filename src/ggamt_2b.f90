
! Copyright (C) 2009 T. McQueen and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: ggamt_2b
! !INTERFACE:
subroutine ggamt_2b(is,g2rho,gvrho,vx,vc,dxdg2,dcdg2)
! !USES:
use modmain
! !DESCRIPTION:
!   Spin-unpolarised version of {\tt ggamt\_sp\_2b}.
!
! !REVISION HISTORY:
!   Created November 2009 (JKD and TMcQ)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: is
real(8), intent(in) :: g2rho(lmmaxvr,nrmtmax)
real(8), intent(in) :: gvrho(lmmaxvr,nrmtmax,3)
real(8), intent(inout) :: vx(lmmaxvr,nrmtmax)
real(8), intent(inout) :: vc(lmmaxvr,nrmtmax)
real(8), intent(in) :: dxdg2(lmmaxvr,nrmtmax)
real(8), intent(in) :: dcdg2(lmmaxvr,nrmtmax)
! local variables
integer nr,nri,i
! allocatable arrays
real(8), allocatable :: rfmt1(:,:),rfmt2(:,:),grfmt(:,:,:)
allocate(rfmt1(lmmaxvr,nrmtmax),rfmt2(lmmaxvr,nrmtmax))
allocate(grfmt(lmmaxvr,nrmtmax,3))
nr=nrmt(is)
nri=nrmtinr(is)
!------------------!
!     exchange     !
!------------------!
! convert dxdg2 to spherical harmonics
call rfsht(nr,nri,1,dxdg2,1,rfmt1)
! compute grad dxdg2
call gradrfmt(nr,nri,rsp(:,is),rfmt1,nrmtmax,grfmt)
! (grad dxdg2).(grad rho) in spherical coordinates
rfmt1(:,1:nr)=0.d0
do i=1,3
  call rbsht(nr,nri,1,grfmt(:,:,i),1,rfmt2)
  rfmt1(:,1:nr)=rfmt1(:,1:nr)+rfmt2(:,1:nr)*gvrho(:,1:nr,i)
end do
vx(:,1:nr)=vx(:,1:nr)-2.d0*(rfmt1(:,1:nr)+dxdg2(:,1:nr)*g2rho(:,1:nr))
!---------------------!
!     correlation     !
!---------------------!
! convert dcdg2 to spherical harmonics
call rfsht(nr,nri,1,dcdg2,1,rfmt1)
! compute grad dcdg2
call gradrfmt(nr,nri,rsp(:,is),rfmt1,nrmtmax,grfmt)
! (grad dcdg2).(grad rho) in spherical coordinates
rfmt1(:,1:nr)=0.d0
do i=1,3
  call rbsht(nr,nri,1,grfmt(:,:,i),1,rfmt2)
  rfmt1(:,1:nr)=rfmt1(:,1:nr)+rfmt2(:,1:nr)*gvrho(:,1:nr,i)
end do
vc(:,1:nr)=vc(:,1:nr)-2.d0*(rfmt1(:,1:nr)+dcdg2(:,1:nr)*g2rho(:,1:nr))
deallocate(rfmt1,rfmt2,grfmt)
return
end subroutine
!EOC

