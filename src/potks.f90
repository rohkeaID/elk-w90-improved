
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
integer is,ias,ir
real(8) ts0,ts1
call timesec(ts0)
! compute the Coulomb potential
call potcoul
! compute the exchange-correlation potential and fields
call potxc
! effective potential from sum of Coulomb and exchange-correlation potentials
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(is,ir)
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
! inner part of muffin-tin: l <= lmaxinr
  do ir=1,nrmtinr(is)
    vsmt(1:lmmaxinr,ir,ias)=vclmt(1:lmmaxinr,ir,ias)+vxcmt(1:lmmaxinr,ir,ias)
    vsmt(lmmaxinr+1:lmmaxvr,ir,ias)=0.d0
  end do
! outer part of muffin-tin: l <= lmaxvr
  do ir=nrmtinr(is)+1,nrmt(is)
    vsmt(:,ir,ias)=vclmt(:,ir,ias)+vxcmt(:,ir,ias)
  end do
end do
!$OMP END DO
!$OMP END PARALLEL
!$OMP PARALLEL WORKSHARE
vsir(:)=vclir(:)+vxcir(:)
!$OMP END PARALLEL WORKSHARE
! generate the Kohn-Sham effective magnetic fields
call bfieldks
call timesec(ts1)
timepot=timepot+ts1-ts0
return
end subroutine
!EOC
