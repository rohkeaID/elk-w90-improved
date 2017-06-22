
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dyntask(fnum,fext)
use modmain
use modphonon
use modmpi
implicit none
! arguments
integer, intent(in) :: fnum
character(*), intent(out) :: fext
! local variables
logical exist
! only master process should search for file
if (.not.mp_mpi) goto 10
do ipph=1,3
  do isph=1,nspecies
    do iaph=1,natoms(isph)
      do iqph=1,nqpt
! construct the phonon file extension
        call phfext(iqph,isph,iaph,ipph,fext)
! determine if the DYN file with this extension exists
        inquire(file='DYN'//trim(fext),exist=exist)
        if (.not.exist) then
          open(fnum,file='DYN'//trim(fext),action='WRITE',form='FORMATTED')
          iasph=idxas(iaph,isph)
          goto 10
        end if
      end do
    end do
  end do
end do
iqph=0; isph=0; iaph=0; iasph=0; ipph=0
write(*,*)
write(*,'("Info(dyntask): nothing more to do")')
10 continue
! broadcast to all other processes
call mpi_bcast(iqph,1,mpi_integer,0,mpi_comm_kpt,ierror)
call mpi_bcast(isph,1,mpi_integer,0,mpi_comm_kpt,ierror)
call mpi_bcast(iaph,1,mpi_integer,0,mpi_comm_kpt,ierror)
call mpi_bcast(iasph,1,mpi_integer,0,mpi_comm_kpt,ierror)
call mpi_bcast(ipph,1,mpi_integer,0,mpi_comm_kpt,ierror)
if (iqph.eq.0) then
  fext='.OUT'
else
  call phfext(iqph,isph,iaph,ipph,fext)
end if
! set the q=0 flag
if (iqph.eq.iq0) then
  tphiq0=.true.
else
  tphiq0=.false.
end if
return
end subroutine

