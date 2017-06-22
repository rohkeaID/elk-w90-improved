
! Copyright (C) 2014 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rfmtaddc
! !INTERFACE:
subroutine rfmtaddc(nr,nri,c,rfmt1,rfmt2)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   nr    : number of radial mesh points (in,integer)
!   nri   : number of points on the inner part of the muffin-tin (in,integer)
!   c     : real constant (in,real)
!   rfmt1 : first real muffin-tin function (in,real(lmmaxvr,nr))
!   rfmt2 : second real muffin-tin function (inout,real(lmmaxvr,nr))
! !DESCRIPTION:
!   Adds a real muffin-tin function times a constant to another:
!   $f_2\rightarrow f_2+c f_1$.
!
! !REVISION HISTORY:
!   Created August 2014 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: nr,nri
real(8), intent(in) :: c
real(8), intent(in) :: rfmt1(lmmaxvr,nr)
real(8), intent(inout) :: rfmt2(lmmaxvr,nr)
! local variables
integer nro,iro,ir
! add on inner part of muffin-tin
do ir=1,nri
  call daxpy(lmmaxinr,c,rfmt1(:,ir),1,rfmt2(:,ir),1)
end do
! add on outer part of muffin-tin
nro=nr-nri
if (nro.eq.0) return
iro=nri+1
call daxpy(lmmaxvr*nro,c,rfmt1(:,iro),1,rfmt2(:,iro),1)
return
end subroutine
!EOC

