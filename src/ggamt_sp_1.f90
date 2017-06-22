
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: ggamt_sp_1
! !INTERFACE:
subroutine ggamt_sp_1(is,rhoup,rhodn,grho,gup,gdn,g2up,g2dn,g3rho,g3up,g3dn)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   is    : species number (in,integer)
!   rhoup : spin-up density in spherical coordinates (in,real(lmmaxvr,nrmtmax))
!   rhodn : spin-down density (in,real(lmmaxvr,nrmtmax))
!   grho  : |grad rho| (out,real(lmmaxvr,nrmtmax))
!   gup   : |grad rhoup| (out,real(lmmaxvr,nrmtmax))
!   gdn   : |grad rhodn| (out,real(lmmaxvr,nrmtmax))
!   g2up  : grad^2 rhoup (out,real(lmmaxvr,nrmtmax))
!   g2dn  : grad^2 rhodn (out,real(lmmaxvr,nrmtmax))
!   g3rho : (grad rho).(grad |grad rho|) (out,real(lmmaxvr,nrmtmax))
!   g3up  : (grad rhoup).(grad |grad rhoup|) (out,real(lmmaxvr,nrmtmax))
!   g3dn  : (grad rhodn).(grad |grad rhodn|) (out,real(lmmaxvr,nrmtmax))
! !DESCRIPTION:
!   Computes $|\nabla\rho|$, $|\nabla\rho^{\uparrow}|$,
!   $|\nabla\rho^{\downarrow}|$, $\nabla^2\rho^{\uparrow}$,
!   $\nabla^2\rho^{\downarrow}$, $\nabla\rho\cdot(\nabla|\nabla\rho|)$,
!   $\nabla\rho^{\uparrow}\cdot(\nabla|\nabla\rho^{\uparrow}|)$ and
!   $\nabla\rho^{\downarrow}\cdot(\nabla|\nabla\rho^{\downarrow}|)$
!   for a muffin-tin charge density, as required by the generalised gradient
!   approximation functionals of type 1 for spin-polarised densities. The input
!   densities and output gradients are in terms of spherical coordinates. See
!   routines {\tt potxc} and {\tt modxcifc}.
!
! !REVISION HISTORY:
!   Created April 2004 (JKD)
!   Simplified and improved, October 2009 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: is
real(8), intent(in) :: rhoup(lmmaxvr,nrmtmax)
real(8), intent(in) :: rhodn(lmmaxvr,nrmtmax)
real(8), intent(out) :: grho(lmmaxvr,nrmtmax)
real(8), intent(out) :: gup(lmmaxvr,nrmtmax)
real(8), intent(out) :: gdn(lmmaxvr,nrmtmax)
real(8), intent(out) :: g2up(lmmaxvr,nrmtmax)
real(8), intent(out) :: g2dn(lmmaxvr,nrmtmax)
real(8), intent(out) :: g3rho(lmmaxvr,nrmtmax)
real(8), intent(out) :: g3up(lmmaxvr,nrmtmax)
real(8), intent(out) :: g3dn(lmmaxvr,nrmtmax)
! local variables
integer nr,nri,i
! allocatable arrays
real(8), allocatable :: rfmt1(:,:),rfmt2(:,:),grfmt(:,:,:)
real(8), allocatable :: gvup(:,:,:),gvdn(:,:,:)
allocate(rfmt1(lmmaxvr,nrmtmax),rfmt2(lmmaxvr,nrmtmax))
allocate(grfmt(lmmaxvr,nrmtmax,3))
allocate(gvup(lmmaxvr,nrmtmax,3),gvdn(lmmaxvr,nrmtmax,3))
nr=nrmt(is)
nri=nrmtinr(is)
!----------------!
!     rho up     !
!----------------!
! convert rhoup to spherical harmonics
call rfsht(nr,nri,1,rhoup,1,rfmt1)
! grad rhoup in spherical coordinates
call gradrfmt(nr,nri,rsp(:,is),rfmt1,nrmtmax,grfmt)
do i=1,3
  call rbsht(nr,nri,1,grfmt(:,:,i),1,gvup(:,:,i))
end do
! |grad rhoup|
gup(:,1:nr)=sqrt(gvup(:,1:nr,1)**2+gvup(:,1:nr,2)**2+gvup(:,1:nr,3)**2)
! grad^2 rhoup in spherical coordinates
call grad2rfmt(nr,nri,rsp(:,is),rfmt1,rfmt2)
call rbsht(nr,nri,1,rfmt2,1,g2up)
! (grad rhoup).(grad |grad rhoup|)
call rfsht(nr,nri,1,gup,1,rfmt1)
call gradrfmt(nr,nri,rsp(:,is),rfmt1,nrmtmax,grfmt)
g3up(:,1:nr)=0.d0
do i=1,3
  call rbsht(nr,nri,1,grfmt(:,:,i),1,rfmt1)
  g3up(:,1:nr)=g3up(:,1:nr)+gvup(:,1:nr,i)*rfmt1(:,1:nr)
end do
!------------------!
!     rho down     !
!------------------!
! convert rhodn to spherical harmonics
call rfsht(nr,nri,1,rhodn,1,rfmt1)
! grad rhodn in spherical coordinates
call gradrfmt(nr,nri,rsp(:,is),rfmt1,nrmtmax,grfmt)
do i=1,3
  call rbsht(nr,nri,1,grfmt(:,:,i),1,gvdn(:,:,i))
end do
gdn(:,1:nr)=sqrt(gvdn(:,1:nr,1)**2+gvdn(:,1:nr,2)**2+gvdn(:,1:nr,3)**2)
! grad^2 rhodn in spherical coordinates
call grad2rfmt(nr,nri,rsp(:,is),rfmt1,rfmt2)
call rbsht(nr,nri,1,rfmt2,1,g2dn)
! (grad rhodn).(grad |grad rhodn|)
call rfsht(nr,nri,1,gdn,1,rfmt1)
call gradrfmt(nr,nri,rsp(:,is),rfmt1,nrmtmax,grfmt)
g3dn(:,1:nr)=0.d0
do i=1,3
  call rbsht(nr,nri,1,grfmt(:,:,i),1,rfmt1)
  g3dn(:,1:nr)=g3dn(:,1:nr)+gvdn(:,1:nr,i)*rfmt1(:,1:nr)
end do
!-------------!
!     rho     !
!-------------!
! |grad rho|
grho(:,1:nr)=sqrt((gvup(:,1:nr,1)+gvdn(:,1:nr,1))**2 &
                 +(gvup(:,1:nr,2)+gvdn(:,1:nr,2))**2 &
                 +(gvup(:,1:nr,3)+gvdn(:,1:nr,3))**2)
! (grad rho).(grad |grad rho|)
call rfsht(nr,nri,1,grho,1,rfmt1)
call gradrfmt(nr,nri,rsp(:,is),rfmt1,nrmtmax,grfmt)
g3rho(:,1:nr)=0.d0
do i=1,3
  call rbsht(nr,nri,1,grfmt(:,:,i),1,rfmt1)
  g3rho(:,1:nr)=g3rho(:,1:nr)+(gvup(:,1:nr,i)+gvdn(:,1:nr,i))*rfmt1(:,1:nr)
end do
deallocate(rfmt1,rfmt2,grfmt,gvup,gvdn)
return
end subroutine
!EOC
