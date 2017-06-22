
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: zfftifc
! !INTERFACE:
subroutine zfftifc(nd,n,sgn,z)
! !INPUT/OUTPUT PARAMETERS:
!   nd   : number of dimensions (in,integer)
!   n    : grid sizes (in,integer(nd))
!   sgn  : FFT direction, -1: forward; 1: backward (in,integer)
!   z    : array to transform (inout,complex(n(1)*n(2)*...*n(nd)))
! !DESCRIPTION:
!   Interface to the double-precision complex fast Fourier transform routine.
!   This is to allow machine-optimised routines to be used without affecting the
!   rest of the code. See routine {\tt nfftifc}.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: nd
integer, intent(in) :: n(nd)
integer, intent(in) :: sgn
complex(8), intent(inout) :: z(*)
! interface to modified FFTPACK5
call cfftnd(nd,n,sgn,z)
return
end subroutine
!EOC

