
! Copyright (C) 2014 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

complex(8) function rzfinp(rfmt,rfir,zfmt,zfir)
use modmain
use modomp
implicit none
! arguments
real(8), intent(in) :: rfmt(npcmtmax,natmtot),rfir(ngtot)
complex(8), intent(in) :: zfmt(npcmtmax,natmtot),zfir(ngtot)
! local variables
integer is,ias,ir,nthd
complex(8) zsum
! external functions
complex(8) rzfmtinp
external rzfmtinp
zsum=0.d0
! interstitial contribution
do ir=1,ngtot
  zsum=zsum+(cfunir(ir)*rfir(ir))*zfir(ir)
end do
zsum=zsum*(omega/dble(ngtot))
! muffin-tin contribution
call omp_hold(natmtot,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(is) REDUCTION(+:zsum) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  zsum=zsum+rzfmtinp(nrcmt(is),nrcmti(is),rcmt(:,is),r2cmt(:,is),rfmt(:,ias), &
   zfmt(:,ias))
end do
!$OMP END DO
!$OMP END PARALLEL
call omp_free(nthd)
rzfinp=zsum
return
end function

