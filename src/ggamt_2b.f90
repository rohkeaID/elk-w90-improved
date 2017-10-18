
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
real(8), intent(in) :: g2rho(npmtmax),gvrho(npmtmax,3)
real(8), intent(inout) :: vx(npmtmax),vc(npmtmax)
real(8), intent(in) :: dxdg2(npmtmax),dcdg2(npmtmax)
! local variables
integer nr,nri,np,i
! allocatable arrays
real(8), allocatable :: rfmt1(:),rfmt2(:),grfmt(:,:)
allocate(rfmt1(npmtmax),rfmt2(npmtmax),grfmt(npmtmax,3))
nr=nrmt(is)
nri=nrmti(is)
np=npmt(is)
!------------------!
!     exchange     !
!------------------!
! convert dxdg2 to spherical harmonics
call rfsht(nr,nri,dxdg2,rfmt1)
! compute grad dxdg2
call gradrfmt(nr,nri,rsp(:,is),rfmt1,npmtmax,grfmt)
! (grad dxdg2).(grad rho) in spherical coordinates
rfmt1(1:np)=0.d0
do i=1,3
  call rbsht(nr,nri,grfmt(:,i),rfmt2)
  rfmt1(1:np)=rfmt1(1:np)+rfmt2(1:np)*gvrho(1:np,i)
end do
vx(1:np)=vx(1:np)-2.d0*(rfmt1(1:np)+dxdg2(1:np)*g2rho(1:np))
!---------------------!
!     correlation     !
!---------------------!
! convert dcdg2 to spherical harmonics
call rfsht(nr,nri,dcdg2,rfmt1)
! compute grad dcdg2
call gradrfmt(nr,nri,rsp(:,is),rfmt1,npmtmax,grfmt)
! (grad dcdg2).(grad rho) in spherical coordinates
rfmt1(1:np)=0.d0
do i=1,3
  call rbsht(nr,nri,grfmt(:,i),rfmt2)
  rfmt1(1:np)=rfmt1(1:np)+rfmt2(1:np)*gvrho(1:np,i)
end do
vc(1:np)=vc(1:np)-2.d0*(rfmt1(1:np)+dcdg2(1:np)*g2rho(1:np))
deallocate(rfmt1,rfmt2,grfmt)
return
end subroutine
!EOC

