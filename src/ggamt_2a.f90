
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
real(8), intent(out) :: g2rho(npmtmax),gvrho(npmtmax,3),grho2(npmtmax)
! local variables
integer is,nr,nri,np,i
! allocatable arrays
real(8), allocatable :: rfmt(:),grfmt(:,:)
allocate(rfmt(npmtmax),grfmt(npmtmax,3))
is=idxis(ias)
nr=nrmt(is)
nri=nrmti(is)
np=npmt(is)
! compute grad^2 rho in spherical coordinates
call grad2rfmt(nr,nri,rsp(:,is),rhomt(:,ias),rfmt)
call rbsht(nr,nri,rfmt,g2rho)
! compute grad rho in spherical coordinates
call gradrfmt(nr,nri,rsp(:,is),rhomt(:,ias),npmtmax,grfmt)
do i=1,3
  call rbsht(nr,nri,grfmt(:,i),gvrho(:,i))
end do
! (grad rho)^2
grho2(1:np)=gvrho(1:np,1)**2+gvrho(1:np,2)**2+gvrho(1:np,3)**2
deallocate(rfmt,grfmt)
return
end subroutine
!EOC

