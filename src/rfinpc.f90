
! Copyright (C) 2016 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

real(8) function rfinpc(ld,rfmt1,rfir1,rfmt2,rfir2)
use modmain
use modomp
implicit none
integer, intent(in) :: ld
real(8), intent(in) :: rfmt1(ld,natmtot),rfir1(ngtot)
real(8), intent(in) :: rfmt2(ld,natmtot),rfir2(ngtot)
! local variables
integer is,ias,ir,nthd
real(8) sum
! external functions
real(8) rfmtinp
external rfmtinp
sum=0.d0
! interstitial contribution
do ir=1,ngtot
  sum=sum+rfir1(ir)*rfir2(ir)*cfunir(ir)
end do
sum=sum*omega/dble(ngtot)
! muffin-tin contribution
call omp_hold(natmtot,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(is) REDUCTION(+:sum) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  sum=sum+rfmtinp(nrcmt(is),nrcmti(is),rcmt(:,is),r2cmt(:,is),rfmt1(:,ias), &
   rfmt2(:,ias))
end do
!$OMP END DO
!$OMP END PARALLEL
call omp_free(nthd)
rfinpc=sum
return
end function

