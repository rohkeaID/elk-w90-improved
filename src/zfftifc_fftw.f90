
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine zfftifc(nd,n,sgn,z)
implicit none
! arguments
integer, intent(in) :: nd,n(nd),sgn
complex(8), intent(inout) :: z(*)
! local variables
integer, parameter :: FFTW_ESTIMATE=64
integer i,p
integer(8) plan
real(8) t1
! interface to FFTW version 3
call dfftw_plan_dft(plan,nd,n,z,z,sgn,FFTW_ESTIMATE)
call dfftw_execute(plan)
call dfftw_destroy_plan(plan)
if (sgn.eq.-1) then
  p=1
  do i=1,nd
    p=p*n(i)
  end do
  t1=1.d0/dble(p)
  call zdscal(p,t1,z,1)
end if
return
end subroutine

