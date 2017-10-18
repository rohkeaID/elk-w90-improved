
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: potks
! !INTERFACE:
subroutine potks
! !USES:
use modmain
! !DESCRIPTION:
!   Computes the Kohn-Sham effective potential by adding together the Coulomb
!   and exchange-correlation potentials. Also computes the effective magnetic
!   field. See routines {\tt potcoul} and {\tt potxc}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ias,np
real(8) ts0,ts1
call timesec(ts0)
! compute the Coulomb potential
call potcoul
! compute the exchange-correlation potential and fields
call potxc
! effective potential from sum of Coulomb and exchange-correlation potentials
do ias=1,natmtot
  is=idxis(ias)
  np=npmt(is)
  vsmt(1:np,ias)=vclmt(1:np,ias)+vxcmt(1:np,ias)
end do
vsir(:)=vclir(:)+vxcir(:)
! generate the Kohn-Sham effective magnetic fields
call bfieldks
call timesec(ts1)
timepot=timepot+ts1-ts0
return
end subroutine
!EOC
