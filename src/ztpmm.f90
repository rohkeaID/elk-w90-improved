
! Copyright (C) 2016 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: ztpmm
! !INTERFACE:
subroutine ztpmm(m,n,a,b,c)
! !INPUT/OUTPUT PARAMETERS:
!   m : order of matrix A (in,integer)
!   n : order of matrix B (in,integer)
!   a : matrix A (in,complex(m,m))
!   b : matrix B (in,complex(n,n))
!   c : matrix C (inout,complex(m*n,m*n))
! !DESCRIPTION:
!   Performs an in-place multiplication of a matrix tensor product with another
!   matrix:
!   $$ C\rightarrow (A\otimes B)C, $$
!   where $A$, $B$ and $C$ are $m\times m$, $n\times n$ and $mn\times mn$
!   general complex matrices, respectively. This is done most efficiently by
!   considering each column of the matrix $C$ as an $m\times n$ matrix $D$,
!   performing the multiplication $BDA^t$ and replacing the column of $C$ with
!   the result.
!
! !REVISION HISTORY:
!   Created December 2016 (TM,JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: m,n
complex(8), intent(in) :: a(m,m),b(n,n)
complex(8), intent(inout) :: c(m*n,*)
! local variables
integer j
complex(8), parameter :: z0=(0.d0,0.d0),z1=(1.d0,0.d0)
! allocatable arrays
complex(8), allocatable :: d(:,:)
allocate(d(n,m))
do j=1,m*n
  call zgemm('N','N',n,m,n,z1,b,n,c(:,j),n,z0,d,n)
  call zgemm('N','T',n,m,m,z1,d,n,a,m,z0,c(:,j),n)
end do
deallocate(d)
return
end subroutine
!EOC

