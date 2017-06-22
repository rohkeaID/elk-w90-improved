
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

! Stub routines for MPI.

module mpi

integer mpi_comm_world
integer mpi_in_place
integer mpi_sum
integer mpi_logical
integer mpi_integer
integer mpi_double_precision
integer mpi_double_complex

end module

subroutine mpi_init(ierror)
implicit none
integer, intent(out) :: ierror
ierror=0
return
end subroutine

subroutine mpi_finalize(ierror)
implicit none
integer, intent(out) :: ierror
ierror=0
return
end subroutine

subroutine mpi_comm_dup(comm,newcomm,ierror)
implicit none
integer, intent(in) :: comm
integer, intent(out) :: newcomm
integer, intent(out) :: ierror
newcomm=comm
ierror=0
return
end subroutine

subroutine mpi_comm_size(comm,size,ierror)
implicit none
integer, intent(in) :: comm
integer, intent(out) :: size
integer, intent(out) :: ierror
size=1
ierror=0
return
end subroutine

subroutine mpi_comm_rank(comm,rank,ierror)
implicit none
integer, intent(in) :: comm
integer, intent(out) :: rank
integer, intent(out) :: ierror
rank=0
ierror=0
return
end subroutine

subroutine mpi_barrier(comm,ierror)
implicit none
integer, intent(in) :: comm
integer, intent(out) :: ierror
ierror=0
return
end subroutine

subroutine mpi_bcast(buffer,count,datatype,root,comm,ierror)
implicit none
real(8), intent(in) :: buffer(*)
integer, intent(in) :: count
integer, intent(in) :: datatype
integer, intent(in) :: root
integer, intent(in) :: comm
integer, intent(out) :: ierror
ierror=0
return
end subroutine

subroutine mpi_allreduce(sendbuf,recvbuf,count,datatype,op,comm,ierror)
implicit none
real(8), intent(in) :: sendbuf(*)
real(8), intent(in) :: recvbuf(*)
integer, intent(in) :: count
integer, intent(in) :: datatype
integer, intent(in) :: op
integer, intent(in) :: comm
integer, intent(out) :: ierror
ierror=0
return
end subroutine

