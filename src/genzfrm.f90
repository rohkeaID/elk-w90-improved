
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genzfrm(wfmt11,wfmt12,wfir11,wfir12,wfmt21,wfmt22,wfir21,wfir22, &
 zrhomt,zrhoir,zmagmt,zmagir)
use modmain
use modomp
implicit none
! arguments
complex(8), intent(in) ::  wfmt11(npcmtmax,natmtot),wfmt12(npcmtmax,natmtot)
complex(8), intent(in) ::  wfir11(ngtot),wfir12(ngtot)
complex(8), intent(in) ::  wfmt21(npcmtmax,natmtot),wfmt22(npcmtmax,natmtot)
complex(8), intent(in) ::  wfir21(ngtot),wfir22(ngtot)
complex(8), intent(out) :: zrhomt(npcmtmax,natmtot),zrhoir(ngtot)
complex(8), intent(out) :: zmagmt(npcmtmax,natmtot,ndmag),zmagir(ngtot,ndmag)
! local variables
integer ld,is,ias,nthd
!-------------------------!
!     muffin-tin part     !
!-------------------------!
ld=npcmtmax*natmtot
call omp_hold(natmtot,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(is) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  call genzrm(npcmt(is),wfmt11(:,ias),wfmt12(:,ias),wfmt21(:,ias), &
   wfmt22(:,ias),zrhomt(:,ias),ld,zmagmt(:,ias,1))
end do
!$OMP END DO
!$OMP END PARALLEL
call omp_free(nthd)
!---------------------------!
!     interstitial part     !
!---------------------------!
call genzrm(ngtot,wfir11,wfir12,wfir21,wfir22,zrhoir,ngtot,zmagir)
return
end subroutine

