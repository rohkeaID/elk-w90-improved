
! Copyright (C) 2014 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dynevs(ev,a,w)
use modmain
use modphonon
implicit none
! arguments
complex(8), intent(in) :: ev(nbph,nbph)
complex(8), intent(in) :: a(nbph,nbph)
real(8), intent(out) :: w(nbph)
! local variables
integer i,j,k
integer lwork,info
real(8) t1,t2
complex(8) z1
! automatic arrays
real(8) wt(nbph)
! allocatable arrays
real(8), allocatable :: rwork(:)
complex(8), allocatable :: work(:)
! external functions
complex(8) zdotc
external zdotc
allocate(rwork(3*nbph))
lwork=2*nbph
allocate(work(lwork))
call zheev('V','U',nbph,a,nbph,w,work,lwork,rwork,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(dynevs): diagonalisation failed")')
  write(*,'(" ZHEEV returned INFO = ",I8)') info
  write(*,*)
  stop
end if
! reorder eigenvalues so that the eigenvectors maximally overlap with ev
wt(:)=w(:)
do i=1,nbph
  j=1
  t1=0.d0
  do k=1,nbph
    z1=zdotc(nbph,ev(:,i),1,a(:,k),1)
    t2=dble(z1)**2+aimag(z1)**2
    if (t2.gt.t1) then
      j=k
      t1=t2
    end if
  end do
  w(i)=wt(j)
end do
deallocate(rwork,work)
return
end subroutine

