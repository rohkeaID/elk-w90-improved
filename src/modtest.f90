
! Copyright (C) 2009 J. K. Dewhurst, S. Sharma and E. K. U. Gross
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modtest

use modmpi

! if test is .true. then the test variables are written to file
logical test

contains

subroutine writetest(id,descr,nv,iv,iva,tol,rv,rva,zv,zva)
implicit none
! arguments
integer, intent(in) :: id
character(*), intent(in) :: descr
integer, optional, intent(in) :: nv
integer, optional, intent(in) :: iv
integer, optional, intent(in) :: iva(*)
real(8), optional, intent(in) :: tol
real(8), optional, intent(in) :: rv
real(8), optional, intent(in) :: rva(*)
complex(8), optional, intent(in) :: zv
complex(8), optional, intent(in) :: zva(*)
! local variables
integer j
character(256) fname
if (.not.test) return
if (.not.mp_mpi) return
if ((id.lt.0).or.(id.gt.999)) then
  write(*,*)
  write(*,'("Error(writetest): id out of range : ",I8)') id
  write(*,*)
  stop
end if
if ((present(iva)).or.(present(rva)).or.(present(zva))) then
  if (.not.present(nv)) then
    write(*,*)
    write(*,'("Error(writetest): missing argument nv")')
    write(*,*)
    stop
  else
    if (nv.le.0) then
      write(*,*)
      write(*,'("Error(writetest): nv <= 0 : ",I8)') nv
      write(*,*)
      stop
    end if
  end if
end if
if ((present(rv)).or.(present(rva)).or.(present(zv)).or.(present(zva))) then
  if (.not.present(tol)) then
    write(*,*)
    write(*,'("Error(writetest): missing argument tol")')
    write(*,*)
    stop
  end if
end if
write(fname,'("TEST",I3.3,".OUT")') id
open(90,file=trim(fname),action='WRITE',form='FORMATTED')
write(90,'("''",A,"''")') trim(descr)
if (present(iv)) then
  write(90,'(2I8)') 1,1
  write(90,'(2I8)') 1,iv
else if (present(rv)) then
  write(90,'(2I8)') 2,1
  write(90,'(G22.12)') tol
  write(90,'(I8,G22.12)') 1,rv
else if (present(zv)) then
  write(90,'(2I8)') 3,1
  write(90,'(G22.12)') tol
  write(90,'(I8,2G22.12)') 1,dble(zv),aimag(zv)
else if (present(iva)) then
  write(90,'(2I8)') 1,nv
  do j=1,nv
    write(90,'(2I8)') j,iva(j)
  end do
else if (present(rva)) then
  write(90,'(2I8)') 2,nv
  write(90,'(G22.12)') tol
  do j=1,nv
    write(90,'(I8,G22.12)') j,rva(j)
  end do
else if (present(zva)) then
  write(90,'(2I8)') 3,nv
  write(90,'(G22.12)') tol
  do j=1,nv
    write(90,'(I8,2G22.12)') j,dble(zva(j)),aimag(zva(j))
  end do
end if
close(90)
return
end subroutine

end module

