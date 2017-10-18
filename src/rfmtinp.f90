
! Copyright (C) 2003-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: rfmtinp
! !INTERFACE:
real(8) function rfmtinp(nr,nri,r,r2,rfmt1,rfmt2)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   nr    : number of radial mesh points (in,integer)
!   nri   : number of radial mesh points on the inner part of the muffin-tin
!           (in,integer)
!   r     : radial mesh (in,real(nr))
!   r2    : r^2 on radial mesh (in,real(nr))
!   rfmt1 : first real function inside muffin-tin (in,real(*))
!   rfmt2 : second real function inside muffin-tin (in,real(*))
! !DESCRIPTION:
!   Calculates the inner product of two real functions in the muffin-tin. So
!   given two real functions of the form
!   $$ f({\bf r})=\sum_{l=0}^{l_{\rm max}}\sum_{m=-l}^{l}f_{lm}(r)R_{lm}
!    (\hat{\bf r}) $$
!   where $R_{lm}$ are the real spherical harmonics, the function returns
!   $$ I=\int\sum_{l=0}^{l_{\rm max}}\sum_{m=-l}^{l}f_{lm}^1(r)f_{lm}^2(r)r^2
!    dr\;. $$
!   The radial integral is performed using a high accuracy cubic spline method.
!   See routines {\tt genrlm} and {\tt fintgt}.
!
! !REVISION HISTORY:
!   Created November 2003 (Sharma)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: nr,nri
real(8), intent(in) :: r(nr),r2(nr)
real(8), intent(in) :: rfmt1(*),rfmt2(*)
! local variables
integer ir,i0,i1
! automatic arrays
real(8) fr(nr)
! external functions
real(8) fintgt
external fintgt
! inner part of muffin-tin
i1=0
do ir=1,nri
  i0=i1+1
  i1=i1+lmmaxi
  fr(ir)=dot_product(rfmt1(i0:i1),rfmt2(i0:i1))*r2(ir)
end do
! outer part of muffin-tin
do ir=nri+1,nr
  i0=i1+1
  i1=i1+lmmaxo
  fr(ir)=dot_product(rfmt1(i0:i1),rfmt2(i0:i1))*r2(ir)
end do
! integrate
rfmtinp=fintgt(-1,nr,r,fr)
return
end function
!EOC

