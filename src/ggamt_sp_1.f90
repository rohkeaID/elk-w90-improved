
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
!   rhoup : spin-up density in spherical coordinates (in,real(npmtmax))
!   rhodn : spin-down density (in,real(npmtmax))
!   grho  : |grad rho| (out,real(npmtmax))
!   gup   : |grad rhoup| (out,real(npmtmax))
!   gdn   : |grad rhodn| (out,real(npmtmax))
!   g2up  : grad^2 rhoup (out,real(npmtmax))
!   g2dn  : grad^2 rhodn (out,real(npmtmax))
!   g3rho : (grad rho).(grad |grad rho|) (out,real(npmtmax))
!   g3up  : (grad rhoup).(grad |grad rhoup|) (out,real(npmtmax))
!   g3dn  : (grad rhodn).(grad |grad rhodn|) (out,real(npmtmax))
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
real(8), intent(in) :: rhoup(npmtmax),rhodn(npmtmax)
real(8), intent(out) :: grho(npmtmax),gup(npmtmax),gdn(npmtmax)
real(8), intent(out) :: g2up(npmtmax),g2dn(npmtmax)
real(8), intent(out) :: g3rho(npmtmax),g3up(npmtmax),g3dn(npmtmax)
! local variables
integer nr,nri,np,i
! allocatable arrays
real(8), allocatable :: rfmt1(:),rfmt2(:),grfmt(:,:)
real(8), allocatable :: gvup(:,:),gvdn(:,:)
allocate(rfmt1(npmtmax),rfmt2(npmtmax))
allocate(grfmt(npmtmax,3),gvup(npmtmax,3),gvdn(npmtmax,3))
nr=nrmt(is)
nri=nrmti(is)
np=npmt(is)
!----------------!
!     rho up     !
!----------------!
! convert rhoup to spherical harmonics
call rfsht(nr,nri,rhoup,rfmt1)
! grad rhoup in spherical coordinates
call gradrfmt(nr,nri,rsp(:,is),rfmt1,npmtmax,grfmt)
do i=1,3
  call rbsht(nr,nri,grfmt(:,i),gvup(:,i))
end do
! |grad rhoup|
gup(1:np)=sqrt(gvup(1:np,1)**2+gvup(1:np,2)**2+gvup(1:np,3)**2)
! grad^2 rhoup in spherical coordinates
call grad2rfmt(nr,nri,rsp(:,is),rfmt1,rfmt2)
call rbsht(nr,nri,rfmt2,g2up)
! (grad rhoup).(grad |grad rhoup|)
call rfsht(nr,nri,gup,rfmt1)
call gradrfmt(nr,nri,rsp(:,is),rfmt1,npmtmax,grfmt)
g3up(1:np)=0.d0
do i=1,3
  call rbsht(nr,nri,grfmt(:,i),rfmt1)
  g3up(1:np)=g3up(1:np)+gvup(1:np,i)*rfmt1(1:np)
end do
!------------------!
!     rho down     !
!------------------!
! convert rhodn to spherical harmonics
call rfsht(nr,nri,rhodn,rfmt1)
! grad rhodn in spherical coordinates
call gradrfmt(nr,nri,rsp(:,is),rfmt1,npmtmax,grfmt)
do i=1,3
  call rbsht(nr,nri,grfmt(:,i),gvdn(:,i))
end do
gdn(1:np)=sqrt(gvdn(1:np,1)**2+gvdn(1:np,2)**2+gvdn(1:np,3)**2)
! grad^2 rhodn in spherical coordinates
call grad2rfmt(nr,nri,rsp(:,is),rfmt1,rfmt2)
call rbsht(nr,nri,rfmt2,g2dn)
! (grad rhodn).(grad |grad rhodn|)
call rfsht(nr,nri,gdn,rfmt1)
call gradrfmt(nr,nri,rsp(:,is),rfmt1,npmtmax,grfmt)
g3dn(1:np)=0.d0
do i=1,3
  call rbsht(nr,nri,grfmt(:,i),rfmt1)
  g3dn(1:np)=g3dn(1:np)+gvdn(1:np,i)*rfmt1(1:np)
end do
!-------------!
!     rho     !
!-------------!
! |grad rho|
grho(1:np)=sqrt((gvup(1:np,1)+gvdn(1:np,1))**2 &
               +(gvup(1:np,2)+gvdn(1:np,2))**2 &
               +(gvup(1:np,3)+gvdn(1:np,3))**2)
! (grad rho).(grad |grad rho|)
call rfsht(nr,nri,grho,rfmt1)
call gradrfmt(nr,nri,rsp(:,is),rfmt1,npmtmax,grfmt)
g3rho(1:np)=0.d0
do i=1,3
  call rbsht(nr,nri,grfmt(:,i),rfmt1)
  g3rho(1:np)=g3rho(1:np)+(gvup(1:np,i)+gvdn(1:np,i))*rfmt1(1:np)
end do
deallocate(rfmt1,rfmt2,grfmt,gvup,gvdn)
return
end subroutine
!EOC

