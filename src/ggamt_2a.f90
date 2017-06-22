
! Copyright (C) 2009 T. McQueen and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: ggamt_2a
! !INTERFACE:
subroutine ggamt_2a(ias,g2rho,gvrho,grho2)
! !USES:
use modmain
! !DESCRIPTION:
!   Spin-unpolarised version of {\tt ggamt\_sp\_2a}.
!
! !REVISION HISTORY:
!   Created November 2009 (JKD and TMcQ)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ias
real(8), intent(out) :: g2rho(lmmaxvr,nrmtmax)
real(8), intent(out) :: gvrho(lmmaxvr,nrmtmax,3)
real(8), intent(out) :: grho2(lmmaxvr,nrmtmax)
! local variables
integer is,nr,nri,i
! allocatable arrays
real(8), allocatable :: rfmt(:,:),grfmt(:,:,:)
allocate(rfmt(lmmaxvr,nrmtmax),grfmt(lmmaxvr,nrmtmax,3))
is=idxis(ias)
nr=nrmt(is)
nri=nrmtinr(is)
! compute grad^2 rho in spherical coordinates
call grad2rfmt(nr,nri,rsp(:,is),rhomt(:,:,ias),rfmt)
call rbsht(nr,nri,1,rfmt,1,g2rho)
! compute grad rho in spherical coordinates
call gradrfmt(nr,nri,rsp(:,is),rhomt(:,:,ias),nrmtmax,grfmt)
do i=1,3
  call rbsht(nr,nri,1,grfmt(:,:,i),1,gvrho(:,:,i))
end do
! (grad rho)^2
grho2(:,1:nr)=gvrho(:,1:nr,1)**2+gvrho(:,1:nr,2)**2+gvrho(:,1:nr,3)**2
deallocate(rfmt,grfmt)
return
end subroutine
!EOC

