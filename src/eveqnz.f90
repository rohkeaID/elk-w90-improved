
! Copyright (C) 2016 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine eveqnz(n,ld,a,w)
use modmain
implicit none
! arguments
integer, intent(in) :: n,ld
complex(8), intent(inout) :: a(ld,n)
real(8), intent(out) :: w(n)
! local variables
integer liwork,lrwork,lwork,info
! allocatable arrays
integer, allocatable :: iwork(:)
real(8), allocatable :: rwork(:)
complex(8), allocatable :: work(:)
select case(evtype)
case(0)
! use the LAPACK routine zheev
  lwork=2*n
  allocate(rwork(3*n),work(lwork))
  call zheev('V','U',n,a,ld,w,work,lwork,rwork,info)
  if (info.ne.0) then
    write(*,*)
    write(*,'("Error(eveqnz): diagonalisation failed")')
    write(*,'(" ZHEEV returned INFO = ",I8)') info
    write(*,*)
    stop
  end if
case(1)
! use the divide-and-conquer LAPACK routine zheevd
  liwork=5*n+3
  lrwork=2*n**2+5*n+1
  lwork=n**2+2*n
  allocate(iwork(liwork),rwork(lrwork),work(lwork))
  call zheevd('V','U',n,a,ld,w,work,lwork,rwork,lrwork,iwork,liwork,info)
  if (info.ne.0) then
    write(*,*)
    write(*,'("Error(eveqnz): diagonalisation failed")')
    write(*,'(" ZHEEVD returned INFO = ",I8)') info
    write(*,*)
    stop
  end if
  deallocate(iwork,rwork,work)
case default
  write(*,*)
  write(*,'("Error(eveqnz): evtype not defined : ",I8)') evtype
  write(*,*)
  stop
end select
return
end subroutine

