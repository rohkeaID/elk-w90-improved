
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modvars

use modmain
use modmpi

contains

subroutine delvars
implicit none
if (.not.mp_mpi) return
! delete existing variables file
open(90,file='VARIABLES.OUT')
close(90,status='DELETE')
return
end subroutine

subroutine writevars(vname,l,m,nv,iv,iva,rv,rva,zv,zva,sv,sva)
implicit none
! arguments
character(*), intent(in) :: vname
integer, optional, intent(in) :: l,m
integer, optional, intent(in) :: nv
integer, optional, intent(in) :: iv
integer, optional, intent(in) :: iva(*)
real(8), optional, intent(in) :: rv
real(8), optional, intent(in) :: rva(*)
complex(8), optional, intent(in) :: zv
complex(8), optional, intent(in) :: zva(*)
character(*), optional, intent(in) :: sv
character(*), optional, intent(in) :: sva(*)
! local variables
integer i
if (.not.wrtvars) return
if (.not.mp_mpi) return
if ((present(iva)).or.(present(rva)).or.(present(zva)).or.(present(sva))) then
  if (.not.present(nv)) then
    write(*,*)
    write(*,'("Error(writevars): missing argument nv")')
    write(*,*)
    stop
  else
    if (nv.lt.0) then
      write(*,*)
      write(*,'("Error(writevars): nv < 0 : ",I8)') nv
      write(*,*)
      stop
    end if
  end if
end if
open(90,file='VARIABLES.OUT',position='APPEND',form='FORMATTED')
write(90,*)
write(90,'(A)',advance='NO') trim(vname)
if (present(l)) write(90,'(I8)',advance='NO') l
if (present(m)) write(90,'(I8)',advance='NO') m
write(90,*)
if (present(iv)) then
  write(90,'(2I8)') 1,1
  write(90,'(I8)') iv
else if (present(rv)) then
  write(90,'(2I8)') 2,1
  write(90,'(G22.12)') rv
else if (present(zv)) then
  write(90,'(2I8)') 3,1
  write(90,'(2G22.12)') dble(zv),aimag(zv)
else if (present(sv)) then
  write(90,'(2I8)') 4,1
  write(90,'(A)') trim(sv)
else if (present(iva)) then
  write(90,'(2I8)') 1,nv
  do i=1,nv
    write(90,'(I8)') iva(i)
  end do
else if (present(rva)) then
  write(90,'(2I8)') 2,nv
  do i=1,nv
    write(90,'(G22.12)') rva(i)
  end do
else if (present(zva)) then
  write(90,'(2I8)') 3,nv
  do i=1,nv
    write(90,'(2G22.12)') dble(zva(i)),aimag(zva(i))
  end do
else if (present(sva)) then
  write(90,'(2I8)') 4,nv
  do i=1,nv
    write(90,'(A)') trim(sva(i))
  end do
end if
close(90)
return
end subroutine

end module

