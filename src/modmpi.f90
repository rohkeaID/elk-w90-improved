
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modmpi

use mpi

! MPI communicator for main code
integer mpicom
! number of MPI processes
integer np_mpi
! local MPI process number
integer lp_mpi
! mp_mpi is .true. if the local MPI process is the master (0)
logical mp_mpi
! commonly used error variable
integer ierror

end module

