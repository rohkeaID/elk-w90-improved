
! Copyright (C) 2014 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: eulerrot
! !INTERFACE:
subroutine eulerrot(ang,rot)
! !INPUT/OUTPUT PARAMETERS:
!   ang : Euler angles (alpha, beta, gamma) (in,real(3))
!   rot : rotation matrix (out,real(3,3))
! !DESCRIPTION:
!   Given a set of Euler angles, $(\alpha,\beta,\gamma)$, this routine
!   determines the corresponding $3\times 3$ rotation matrix. The so-called
!   `y-convention' is taken for the Euler angles. See the routine {\tt roteuler}
!   for details.
!
! !REVISION HISTORY:
!   Created January 2014 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: ang(3)
real(8), intent(out) :: rot(3,3)
! local variables
real(8) sa,sb,sg,ca,cb,cg
sa=sin(ang(1)); sb=sin(ang(2)); sg=sin(ang(3))
ca=cos(ang(1)); cb=cos(ang(2)); cg=cos(ang(3))
rot(1,1)=cg*cb*ca-sg*sa
rot(1,2)=cg*cb*sa+sg*ca
rot(1,3)=-cg*sb
rot(2,1)=-sg*cb*ca-cg*sa
rot(2,2)=-sg*cb*sa+cg*ca
rot(2,3)=sg*sb
rot(3,1)=sb*ca
rot(3,2)=sb*sa
rot(3,3)=cb
return
end subroutine
!EOC

