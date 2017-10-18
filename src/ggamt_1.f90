
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
real(8), intent(out) :: grho(npmtmax),g2rho(npmtmax),g3rho(npmtmax)
! local variables
integer is,nr,nri,np,i
! allocatable arrays
real(8), allocatable :: grfmt(:,:),gvrho(:,:),rfmt(:)
allocate(grfmt(npmtmax,3),gvrho(npmtmax,3),rfmt(npmtmax))
is=idxis(ias)
nr=nrmt(is)
nri=nrmti(is)
np=npmt(is)
! |grad rho|
call gradrfmt(nr,nri,rsp(:,is),rhomt(:,ias),npmtmax,grfmt)
do i=1,3
  call rbsht(nr,nri,grfmt(:,i),gvrho(:,i))
end do
grho(1:np)=sqrt(gvrho(1:np,1)**2+gvrho(1:np,2)**2+gvrho(1:np,3)**2)
! grad^2 rho in spherical coordinates
call grad2rfmt(nr,nri,rsp(:,is),rhomt(:,ias),rfmt)
call rbsht(nr,nri,rfmt,g2rho)
! (grad rho).(grad |grad rho|)
call rfsht(nr,nri,grho,rfmt)
call gradrfmt(nr,nri,rsp(:,is),rfmt,npmtmax,grfmt)
g3rho(1:np)=0.d0
do i=1,3
  call rbsht(nr,nri,grfmt(:,i),rfmt)
  g3rho(1:np)=g3rho(1:np)+gvrho(1:np,i)*rfmt(1:np)
end do
deallocate(grfmt,gvrho,rfmt)
return
end subroutine
!EOC

