
! Copyright (C) 2009 J. K. Dewhurst, S. Sharma and E. K. U. Gross
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

real(8) function sdelta_lr(x)
implicit none
! arguments
real(8), intent(in) :: x
! local variables
real(8), parameter :: twopi=6.2831853071795864769d0
sdelta_lr=1.d0/(twopi*(x**2+0.25d0))
return
end function

