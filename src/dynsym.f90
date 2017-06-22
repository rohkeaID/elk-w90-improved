
! Copyright (C) 2006-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dynsym(vpl,dynp)
use modmain
use modphonon
implicit none
! arguments
real(8), intent(in) :: vpl(3)
complex(8), intent(inout) :: dynp(nbph,nbph)
! local variables
integer isym,lspl,i,j,n
real(8) v1(3),v2(3),s(3,3),t1
! automatic arrays
complex(8) dyns(nbph,nbph)
! map input vector to first Brillouin zone
v1(:)=vpl(:)
call vecfbz(epslat,bvec,v1)
n=0
dyns(:,:)=0.d0
! use the symmetries which leave vpl invariant
do isym=1,nsymcrys
  if (.not.tvzsymc(isym)) cycle
  lspl=lsplsymc(isym)
  s(:,:)=dble(symlat(:,:,lspl))
  call r3mtv(s,v1,v2)
  t1=abs(v1(1)-v2(1))+abs(v1(2)-v2(2))+abs(v1(3)-v2(3))
  if (t1.lt.epslat) then
    call dynsymapp(isym,v1,dynp,dyns)
    n=n+1
  end if
end do
if (n.eq.0) then
  write(*,*)
  write(*,'("Error(dynsym): no symmetries leave vpl invariant")')
  write(*,*)
  stop
end if
t1=1.d0/dble(n)
dynp(:,:)=t1*dyns(:,:)
! make the matrix Hermitian
do i=1,nbph
  do j=i,nbph
    dynp(i,j)=0.5d0*(dynp(i,j)+conjg(dynp(j,i)))
    dynp(j,i)=conjg(dynp(i,j))
  end do
end do
return
end subroutine

