
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: zfinp
! !INTERFACE:
complex(8) function zfinp(zfmt1,zfir1,zfmt2,zfir2)
! !USES:
use modmain
use modomp
! !INPUT/OUTPUT PARAMETERS:
!   zfmt1 : first complex function in spherical harmonics/coordinates for all
!           muffin-tins (in,complex(npcmtmax,natmtot))
!   zfir1 : first complex interstitial function in real-space
!           (in,complex(ngtc))
!   zfmt2 : second complex function in spherical harmonics/coordinates for all
!           muffin-tins (in,complex(npcmtmax,natmtot))
!   zfir2 : second complex interstitial function in real-space
!           (in,complex(ngtc))
! !DESCRIPTION:
!   Calculates the inner product of two complex fuctions over the entire unit
!   cell. The muffin-tin functions should be stored on the coarse radial grid.
!   In the interstitial region, the integrand is multiplied with the
!   characteristic function to remove the contribution from the muffin-tin. See
!   routines {\tt zfmtinp} and {\tt gencfun}.
!
! !REVISION HISTORY:
!   Created July 2004 (Sharma)
!EOP
!BOC
implicit none
! arguments
complex(8), intent(in) :: zfmt1(npcmtmax,natmtot),zfir1(ngtc)
complex(8), intent(in) :: zfmt2(npcmtmax,natmtot),zfir2(ngtc)
! local variables
integer is,ias,ir,nthd
complex(8) zsum
! external functions
complex(8) zfmtinp
external zfmtinp
! interstitial contribution
zsum=cfrc(1)*conjg(zfir1(1))*zfir2(1)
do ir=2,ngtc
  zsum=zsum+cfrc(ir)*conjg(zfir1(ir))*zfir2(ir)
end do
zsum=zsum*(omega/dble(ngtc))
! muffin-tin contribution
call omp_hold(natmtot,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(is) REDUCTION(+:zsum) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  zsum=zsum+zfmtinp(nrcmt(is),nrcmti(is),rcmt(:,is),r2cmt(:,is),zfmt1(:,ias), &
   zfmt2(:,ias))
end do
!$OMP END DO
!$OMP END PARALLEL
call omp_free(nthd)
zfinp=zsum
return
end function
!EOC

