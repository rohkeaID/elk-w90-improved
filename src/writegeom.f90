
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writegeom
! !INTERFACE:
subroutine writegeom(fnum)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   fnum : file number for writing output (in,integer)
! !DESCRIPTION:
!   Outputs the lattice vectors and atomic positions to file, in a format
!   which may be then used directly in {\tt elk.in}.
!
! !REVISION HISTORY:
!   Created January 2004 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: fnum
! local variables
integer is,ia,i
real(8) v1(3),v2(3)
write(fnum,*)
write(fnum,'("scale")')
write(fnum,'(" 1.0")')
write(fnum,*)
write(fnum,'("scale1")')
write(fnum,'(" 1.0")')
write(fnum,*)
write(fnum,'("scale2")')
write(fnum,'(" 1.0")')
write(fnum,*)
write(fnum,'("scale3")')
write(fnum,'(" 1.0")')
write(fnum,*)
write(fnum,'("avec")')
write(fnum,'(3G18.10)') avec(:,1)
write(fnum,'(3G18.10)') avec(:,2)
write(fnum,'(3G18.10)') avec(:,3)
if (molecule) then
  write(fnum,*)
  write(fnum,'("molecule")')
  write(fnum,'(" ",L1)') molecule
end if
write(fnum,*)
write(fnum,'("atoms")')
write(fnum,'(I4,T40," : nspecies")') nspecies
do is=1,nspecies
  write(fnum,'("''",A,"''",T40," : spfname")') trim(spfname(is))
  write(fnum,'(I4,T40," : natoms; atpos, bfcmt below")') natoms(is)
  do ia=1,natoms(is)
    if (molecule) then
! map lattice coordinates to [-0.5,0.5) and write as Cartesian coordinates
      v1(:)=atposl(:,ia,is)
      do i=1,3
        if (v1(i).gt.0.5d0) v1(i)=v1(i)-1.d0
      end do
      call r3mv(avec,v1,v2)
    else
! otherwise write lattice coordinates
      v2(:)=atposl(:,ia,is)
    end if
    write(fnum,'(3F14.8,"  ",3F12.8)') v2(:),bfcmt(:,ia,is)
  end do
end do
return
end subroutine
!EOC

