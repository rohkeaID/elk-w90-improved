
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: zfinp
! !INTERFACE:
complex(8) function zfinp(zfmt1,zfir1,zfmt2,zfir2)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   zfmt1 : first complex function in spherical harmonics/coordinates for all
!           muffin-tins (in,complex(lmmaxvr,nrcmtmax,natmtot))
!   zfir1 : first complex interstitial function in real-space
!           (in,complex(ngtot))
!   zfmt2 : second complex function in spherical harmonics/coordinates for all
!           muffin-tins (in,complex(lmmaxvr,nrcmtmax,natmtot))
!   zfir2 : second complex interstitial function in real-space
!           (in,complex(ngtot))
! !DESCRIPTION:
!   Calculates the inner product of two complex fuctions over the entire unit
!   cell. The muffin-tin functions should be stored on the coarse radial grid
!   and have angular momentum cut-off {\tt lmaxvr}. In the intersitial region,
!   the integrand is multiplied with the characteristic function, to remove the
!   contribution from the muffin-tin. See routines {\tt zfmtinp} and
!   {\tt gencfun}.
!
! !REVISION HISTORY:
!   Created July 2004 (Sharma)
!EOP
!BOC
implicit none
! arguments
complex(8), intent(in) :: zfmt1(lmmaxvr,nrcmtmax,natmtot),zfir1(ngtot)
complex(8), intent(in) :: zfmt2(lmmaxvr,nrcmtmax,natmtot),zfir2(ngtot)
! local variables
integer ias,is,ir
complex(8) zsum
! external functions
complex(8) zfmtinp
external zfmtinp
zsum=0.d0
! interstitial contribution
do ir=1,ngtot
  zsum=zsum+cfunir(ir)*conjg(zfir1(ir))*zfir2(ir)
end do
zsum=zsum*(omega/dble(ngtot))
! muffin-tin contribution
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(is) REDUCTION(+:zsum)
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  zsum=zsum+zfmtinp(nrcmt(is),nrcmtinr(is),rcmt(:,is),r2cmt(:,is), &
   zfmt1(:,:,ias),zfmt2(:,:,ias))
end do
!$OMP END DO
!$OMP END PARALLEL
zfinp=zsum
return
end function
!EOC

