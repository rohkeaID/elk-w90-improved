
! Copyright (C) 2017 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine zftwfir(ngp,igpig,wfir)
use modmain
use modomp
implicit none
! arguments
integer, intent(in) :: ngp(nspnfv),igpig(ngkmax,nspnfv)
complex(8), intent(inout) :: wfir(ngtc,nspinor,nstsv)
! local variables
integer ist,ispn,jspn
integer igp,ifg,nthd
real(8) t0
! allocatable arrays
complex(8), allocatable :: z(:)
t0=1.d0/sqrt(omega)
call omp_hold(nstsv,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(z,ispn,jspn,igp,ifg) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ist=1,nstsv
  allocate(z(ngkmax))
  do ispn=1,nspinor
    jspn=jspnfv(ispn)
    call zcopy(ngp(jspn),wfir(:,ispn,ist),1,z,1)
    wfir(:,ispn,ist)=0.d0
    do igp=1,ngp(jspn)
      ifg=igfc(igpig(igp,jspn))
      wfir(ifg,ispn,ist)=t0*z(igp)
    end do
    call zfftifc(3,ngdc,1,wfir(:,ispn,ist))
  end do
  deallocate(z)
end do
!$OMP END DO
!$OMP END PARALLEL
call omp_free(nthd)
return
end subroutine

