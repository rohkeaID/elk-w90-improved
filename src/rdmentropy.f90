
! Copyright (C) 2008 T. Baldsiefen, S. Sharma, J. K. Dewhurst and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

!BOP
! !ROUTINE: rdmentropy
! !INTERFACE:
subroutine rdmentropy
! !USES:
use modmain
use modrdm
! !DESCRIPTION:
!  Calculates RDMFT entropy $S=-\sum_i n_i\log(n_i/n_{\rm max})
!  +(n_{\rm max}-n_i)\log(1-n_i/n_{\rm max})$, where $n_{\rm max}$ is the
!  maximum allowed occupancy (1 or 2).
!
! !REVISION HISTORY:
!   Created 2008 (Baldsiefen)
!EOP
!BOC
implicit none
! local variables
integer ik,ist
real(8) t1
rdmentrpy=0.d0
do ik=1,nkpt
  do ist=1,nstsv
    t1=max(occsv(ist,ik),epsocc)
    t1=min(t1,occmax-epsocc)
    rdmentrpy=rdmentrpy-wkpt(ik)*(t1*log(t1/occmax) &
     +(occmax-t1)*log(1.d0-t1/occmax))
  end do
end do
rdmentrpy=kboltz*rdmentrpy
return
end subroutine
!EOC
