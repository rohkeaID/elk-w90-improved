
! Copyright (C) 2015 Manh Duc Le, 2017-18 Arsenii Gerasimov, Yaroslav Kvashnin
! and Lars Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: getw90bands
! !INTERFACE:
subroutine getw90bands(str)
! !USES:
use modmain
use modw90
! !INPUT/OUTPUT PARAMETERS:
!   str    : string defining which bands to pass to the Wannier90
!
! !DESCRIPTION:
!   Parses the string from {\tt elk.in} denoting which bands to pass to the Wannier90.
!   This should be a comma separated list of indices (e.g. 1-4 or 1,2,4 or
!   1-3,5-6).
!
! !REVISION HISTORY:
!   Created January 2015 (Manh Duc Le)
!   Modified October 2018 (Yaroslav Kvashnin, Arsenii Gerasimov)
!EOP
!BOC
implicit none
! arguments
character(256), intent(in) :: str
! parameters
integer, parameter :: max_bands = 1024  ! Maximum num bands to pass to Wannier90
! local variables
integer :: i0,i1,idf,idb,ib,ibs,ibn,ibi,ios
integer bands(max_bands)
logical end_next
!-------------------------------------------------------------------------------

! Band index initialisation
ib = 1
! Split str into comma delimited sections and parses them.
i0 = 1
i1 = index(str,',')
if( i1 .eq. 0 ) then
  if( index(str,'-') .eq. 0 ) then
    wann_nband = 1
    allocate(wann_bands(wann_nband))
    read(str,*,iostat=ios,err=100) wann_bands(1)
    return
  else
    i1 = len(str)
  end if
end if

end_next = .false.
do
  idf = index(str(i0:i1),'-')
  idb = index(str(i0:i1),'-',.true.)
  if( idf .ne. idb ) then
    write(*,*)
    write(*,'("Error(readinput): more than one dash given for a band range&
                      & in wann_projections for band definition: ",A)') str
    write(*,*)
    stop
  end if
  if( idf .eq. 0 ) then
    read(str(i0:i1),*,iostat=ios,err=100) bands(ib)
    ib = ib + 1
  else
    read(str(i0:(idf+i0-2)),*,iostat=ios,err=100) ibs
    read(str((idf+i0):i1),*,iostat=ios,err=100)   ibn
    do ibi = ibs,ibn
      bands(ib) = ibi
      ib = ib + 1
    end do
  end if
  if ( end_next ) exit
  i0 = i1+1
  if ( i0 .gt. len(str) ) exit
  i1 = index(str(i0:len(str)),',') + i0 - 1
  if ( i1 .eq. i0 - 1 ) then
    end_next = .true.
    i1 = len(str)
  end if
end do
wann_nband = ib - 1
allocate(wann_bands(wann_nband))
wann_bands(1:wann_nband) = bands(1:wann_nband)

return

100 write(*,*)
    write(*,'("Error(readinput): error reading the band definition: ''",A,"''&
                                                  & in wann_projections")') str
    write(*,*)
    stop

end subroutine
!EOC
