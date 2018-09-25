
! Copyright (C) 2017 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine pades(ns,r,ni,zi,ui,no,zo,uo)
implicit none
! arguments
integer, intent(in) :: ns
real(8), intent(in) :: r
integer, intent(in) :: ni
complex(8), intent(in) :: zi(ni)
complex(8), intent(in) :: ui(ni)
integer, intent(in) :: no
complex(8), intent(in) :: zo(no)
complex(8), intent(out) :: uo(no)
! local variables
integer i
real(8), parameter :: pi=3.1415926535897932385d0
real(8) t1,t2
complex(8) z1
! allocatable arrays
complex(8), allocatable :: u1(:),u2(:)
if (ns.le.0) then
  write(*,*)
  write(*,'("Error(pades): ns <= 0 : ",I8)') ns
  write(*,*)
  stop
end if
if (ns.eq.1) then
  call pade(ni,zi,ui,no,zo,uo)
  return
end if
allocate(u1(ni),u2(no))
uo(:)=0.d0
do i=1,ns
  t1=dble(i-1)/dble(ns)
  t2=6.d0*pi*t1
  z1=r*t1*cmplx(cos(t2),sin(t2),8)
  u1(:)=ui(:)+z1
  call pade(ni,zi,u1,no,zo,u2)
  uo(:)=uo(:)+u2(:)-z1
end do
t1=1.d0/dble(ns)
uo(:)=t1*uo(:)
deallocate(u1,u2)
return
end subroutine

