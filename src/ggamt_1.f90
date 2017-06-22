
! Copyright (C) 2009 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: ggamt_1
! !INTERFACE:
subroutine ggamt_1(ias,grho,g2rho,g3rho)
! !USES:
use modmain
! !DESCRIPTION:
!   Spin-unpolarised version of {\tt ggamt\_sp\_1}.
!
! !REVISION HISTORY:
!   Created November 2009 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ias
real(8), intent(out) :: grho(lmmaxvr,nrmtmax)
real(8), intent(out) :: g2rho(lmmaxvr,nrmtmax)
real(8), intent(out) :: g3rho(lmmaxvr,nrmtmax)
! local variables
integer is,nr,nri,i
! allocatable arrays
real(8), allocatable :: grfmt(:,:,:),gvrho(:,:,:),rfmt(:,:)
allocate(grfmt(lmmaxvr,nrmtmax,3),gvrho(lmmaxvr,nrmtmax,3))
allocate(rfmt(lmmaxvr,nrmtmax))
is=idxis(ias)
nr=nrmt(is)
nri=nrmtinr(is)
! |grad rho|
call gradrfmt(nr,nri,rsp(:,is),rhomt(:,:,ias),nrmtmax,grfmt)
do i=1,3
  call rbsht(nr,nri,1,grfmt(:,:,i),1,gvrho(:,:,i))
end do
grho(:,1:nr)=sqrt(gvrho(:,1:nr,1)**2+gvrho(:,1:nr,2)**2+gvrho(:,1:nr,3)**2)
! grad^2 rho in spherical coordinates
call grad2rfmt(nr,nri,rsp(:,is),rhomt(:,:,ias),rfmt)
call rbsht(nr,nri,1,rfmt,1,g2rho)
! (grad rho).(grad |grad rho|)
call rfsht(nr,nri,1,grho,1,rfmt)
call gradrfmt(nr,nri,rsp(:,is),rfmt,nrmtmax,grfmt)
g3rho(:,1:nr)=0.d0
do i=1,3
  call rbsht(nr,nri,1,grfmt(:,:,i),1,rfmt)
  g3rho(:,1:nr)=g3rho(:,1:nr)+gvrho(:,1:nr,i)*rfmt(:,1:nr)
end do
deallocate(grfmt,gvrho,rfmt)
return
end subroutine
!EOC

