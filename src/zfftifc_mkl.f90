
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine zfftifc(nd,n,sgn,z)
use mkl_dfti
implicit none
! arguments
integer, intent(in) :: nd
integer, intent(in) :: n(nd)
integer, intent(in) :: sgn
complex(8), intent(inout) :: z(*)
! local variables
integer status,i,p
real(8) t1
type(DFTI_DESCRIPTOR), pointer :: handle
! interface to the Intel MKL advanced Discreet Fourier Transform (DFT) routines
! (with thanks to Torbjorn Bjorkman)
p=1
do i=1,nd
  p=p*n(i)
end do
t1=1.d0/dble(p)
status=DftiCreateDescriptor(handle,DFTI_DOUBLE,DFTI_COMPLEX,nd,n)
status=DftiSetValue(handle,DFTI_FORWARD_SCALE,t1)
status=DftiCommitDescriptor(handle)
if (sgn.eq.-1) then
  status=DftiComputeForward(handle,z)
else
  status=DftiComputeBackward(handle,z)
end if
status=DftiFreeDescriptor(handle)
return
end subroutine

