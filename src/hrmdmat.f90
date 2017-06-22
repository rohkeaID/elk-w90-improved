
! Copyright (C) 2014 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine hrmdmat(n,dmat)
implicit none
! arguments
integer, intent(in) :: n
complex(8), intent(inout) :: dmat(n,n)
! local variables
integer i,j
complex(8) z1
do i=1,n
  do j=1,i-1
    z1=dmat(i,j)+conjg(dmat(j,i))
    dmat(i,j)=z1
    dmat(j,i)=conjg(z1)
  end do
  dmat(i,i)=dble(dmat(i,i))
end do
return
end subroutine

