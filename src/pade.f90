
! Copyright (C) A. Sanna and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: pade
! !INTERFACE:
subroutine pade(nin,zin,uin,nout,zout,uout)
! !INPUT/OUTPUT PARAMETERS:
!   nin  : number of input points (in,integer)
!   zin  : input points (in,complex(nin))
!   uin  : input function values (in,complex(nin))
!   nout : number of output points (in,integer)
!   zout : output points (in,complex(nout))
!   uout : output function values (out,complex(nout))
! !DESCRIPTION:
!   Calculates a Pad\'{e} approximant of a function, given the function
!   evaluated on a set of points in the complex plane. The function is returned
!   for a set of complex output points. The algorithm from H. J. Vidberg and
!   J. W. Serene {\it J. Low Temp. Phys.} {\bf } 29, 179 (1977) is used.
!
! !REVISION HISTORY:
!   Created December 2010 (Antonio Sanna)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: nin
complex(8), intent(in) :: zin(nin)
complex(8), intent(in) :: uin(nin)
integer, intent(in) :: nout
complex(8), intent(in) :: zout(nout)
complex(8), intent(out) :: uout(nout)
! local variables
integer i,j
real(8) t1
complex(8) a0,a1,b0,b1,z1,z2
! allocatable arrays
complex(8), allocatable :: g(:,:)
if ((nin.le.0).or.(nout.le.0)) then
  write(*,*)
  write(*,'("Error(pade): invalid number of input or output points : ",2I8)') &
   nin,nout
  write(*,*)
  stop
end if
allocate(g(nin,nin))
! define the g functions using Eq. (A2)
g(1,:)=uin(:)
do i=2,nin
  do j=i,nin
    g(i,j)=(g(i-1,i-1)-g(i-1,j))/((zin(j)-zin(i-1))*g(i-1,j))
  end do
end do
! loop over output points
do i=1,nout
! use recursive algorithm in Eq. (A3) to evaluate function
  a0=0.d0
  a1=g(1,1)
  b0=1.d0
  b1=1.d0
  do j=2,nin
    z1=(zout(i)-zin(j-1))*g(j,j)
    z2=a1+z1*a0
    a0=a1
    a1=z2
    z2=b1+z1*b0
    b0=b1
    b1=z2
! check for overflow and rescale
    if ((abs(dble(a1)).gt.1.d100).or.(abs(aimag(a1)).gt.1.d100)) then
      t1=1.d0/abs(a1)
      a0=a0*t1
      b0=b0*t1
      a1=a1*t1
      b1=b1*t1
    end if
    if ((abs(dble(b1)).gt.1.d100).or.(abs(aimag(b1)).gt.1.d100)) then
      t1=1.d0/abs(b1)
      a0=a0*t1
      b0=b0*t1
      a1=a1*t1
      b1=b1*t1
    end if
  end do
  t1=abs(dble(b1))+abs(aimag(b1))
  if (t1.ne.0.d0) then
    uout(i)=a1/b1
  else
    uout(i)=0.d0
  end if
end do
deallocate(g)
return
end subroutine
!EOC

