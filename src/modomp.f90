
! Copyright (C) 2015 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modomp

! number of processors available for OpenMP
integer ncpu
! maximum number of OpenMP threads available
integer maxthd
! maximum OpenMP nesting level
integer maxlvl
! number of active OpenMP threads for each nesting level
integer, allocatable :: nathd(:)
! number of idle OpenMP threads
integer nithd

interface

integer function omp_get_num_procs()
end function

integer function omp_get_max_threads()
end function

integer function omp_get_num_threads()
end function

integer function omp_get_thread_num()
end function

logical function omp_get_nested()
end function

integer function omp_get_max_active_levels()
end function

logical function omp_get_dynamic()
end function

integer function omp_get_level()
end function

subroutine omp_set_num_threads(num_threads)
integer, intent(in) :: num_threads
end subroutine

subroutine omp_set_nested(nested)
logical, intent(in) :: nested
end subroutine

subroutine omp_set_max_active_levels(max_levels)
integer, intent(in) :: max_levels
end subroutine

subroutine omp_set_dynamic(dynamic_threads)
logical, intent(in) :: dynamic_threads
end subroutine

end interface

contains

subroutine omp_init
implicit none
! get number of processors
ncpu=omp_get_num_procs()
if (maxthd.lt.0) then
! set the number of threads equal to the number of processors
  maxthd=ncpu
  call omp_set_num_threads(maxthd)
else if (maxthd.eq.0) then
! use the system default number of threads
  maxthd=omp_get_max_threads()
else
! use the number of threads specified in the input file
  call omp_set_num_threads(maxthd)
end if
! switch off dynamic allocation of threads
call omp_set_dynamic(.false.)
! allow nested parallelism
call omp_set_nested(.true.)
! set the maximum nesting level
if (maxlvl.le.0) then
  maxlvl=omp_get_max_active_levels()
end if
call omp_set_max_active_levels(maxlvl)
if (allocated(nathd)) deallocate(nathd)
allocate(nathd(0:maxlvl))
! initialise the number of active and idle threads
call omp_reset
return
end subroutine

subroutine omp_reset
implicit none
! number of active threads at each nesting level
nathd(0)=1
nathd(1:)=0
! number of idle threads
nithd=maxthd-1
return
end subroutine

subroutine omp_hold(nloop,nthd)
implicit none
! arguments
integer, intent(in) :: nloop
integer, intent(out) :: nthd
! local variables
integer lvl,nl,na,n,i
! nesting level
lvl=omp_get_level()+1
if ((lvl.le.0).or.(lvl.gt.maxlvl)) then
  nthd=1
  return
end if
! the number of loops should be greater than zero
nl=max(nloop,1)
!$OMP CRITICAL(omp_hold_)
! determine number of active threads at the next highest nesting level
na=1
do i=lvl-1,0,-1
  if (nathd(i).gt.0) then
    na=nathd(i)
    exit
  end if
end do
! number of threads allowed for this loop
n=maxthd/na
n=min(max(n,1),nl,nithd+1,maxthd)
nthd=n
! add to number of active threads in this nesting level
n=nathd(lvl)+nthd
n=min(max(n,0),maxthd)
nathd(lvl)=n
! subtract from the number of idle threads
n=nithd-nthd+1
n=min(max(n,0),maxthd)
nithd=n
!$OMP END CRITICAL(omp_hold_)
return
end subroutine

subroutine omp_free(nthd)
implicit none
! arguments
integer, intent(in) :: nthd
! local variables
integer lvl,n
! nesting level
lvl=omp_get_level()+1
if ((lvl.le.0).or.(lvl.gt.maxlvl)) return
!$OMP CRITICAL(omp_free_)
! subtract from the number of active threads in this nesting level
n=nathd(lvl)-nthd
n=min(max(n,0),maxthd)
nathd(lvl)=n
! add to the number of idle threads
n=nithd+nthd-1
n=min(max(n,0),maxthd)
nithd=n
!$OMP END CRITICAL(omp_free_)
return
end subroutine

end module

