
! Copyright (C) 2015 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genmcph(w,ev,a)
use modmain
use modphonon
implicit none
! arguments
real(8), intent(in) :: w(nbph)
complex(8), intent(in) :: ev(nbph,nbph)
complex(8), intent(out) :: a(nbph,nbph)
! local variables
integer is,ia,ip,i,j
real(8) t1
do j=1,nbph
  i=0
  do is=1,nspecies
    t1=2.d0*spmass(is)*w(j)
    if (t1.gt.1.d-8) then
      t1=1.d0/sqrt(t1)
    else
      t1=0.d0
    end if
    do ia=1,natoms(is)
      do ip=1,3
        i=i+1
        a(i,j)=t1*ev(i,j)
      end do
    end do
  end do
end do
return
end subroutine

