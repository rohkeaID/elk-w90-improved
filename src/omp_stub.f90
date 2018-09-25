
! Copyright (C) 2015 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

! Stub routines for OpenMP.

integer function omp_get_num_procs()
implicit none
omp_get_num_procs=1
return
end function

integer function omp_get_max_threads()
implicit none
omp_get_max_threads=1
return
end function

integer function omp_get_num_threads()
implicit none
omp_get_num_threads=1
return
end function

integer function omp_get_thread_num()
implicit none
omp_get_thread_num=0
return
end function

logical function omp_get_nested()
implicit none
omp_get_nested=.false.
return
end function

integer function omp_get_max_active_levels()
implicit none
omp_get_max_active_levels=1
return
end function

logical function omp_get_dynamic()
implicit none
omp_get_dynamic=.false.
return
end function

integer function omp_get_level()
implicit none
omp_get_level=0
return
end function

subroutine omp_set_num_threads(num_threads)
implicit none
integer, intent(in) :: num_threads
return
end subroutine

subroutine omp_set_nested(nested)
implicit none
logical, intent(in) :: nested
return
end subroutine

subroutine omp_set_max_active_levels(max_levels)
implicit none
integer, intent(in) :: max_levels
return
end subroutine

subroutine omp_set_dynamic(dynamic)
implicit none
logical, intent(in) :: dynamic
return
end subroutine

subroutine omp_init_lock(lock)
implicit none
integer(8), intent(in) :: lock
return
end subroutine

subroutine omp_destroy_lock(lock)
implicit none
integer(8), intent(in) :: lock
return
end subroutine

subroutine omp_set_lock(lock)
implicit none
integer(8), intent(in) :: lock
return
end subroutine

subroutine omp_unset_lock(lock)
implicit none
integer(8), intent(in) :: lock
return
end subroutine
