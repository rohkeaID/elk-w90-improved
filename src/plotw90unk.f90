
! Copyright (C) 2017-18 Arsenii Gerasimov, Yaroslav Kvashnin and Lars Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: plotw90unk
! !INTERFACE:
subroutine plotw90unk(fnum,zfmt,zfir)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   fnum : plot file ID (in,integer)
!   zfmt : complex muffin-tin  function (in,complex(npcmtmax,natmtot))
!   zfir : complex intersitial function (in,complex(ngtot))
! !DESCRIPTION:
!   Produces a 3D plot of the complex functions contained in {\tt zfmt} and
!   {\tt zfir} arrays in the parallelepiped defined by the corner vertices in
!   the global array {\tt vclp3d}. Similar routine to the {\tt plot3d}.
!
! !REVISION HISTORY:
!   Created September 2018 (Arsenii Gerasimov)
!EOP
!BOC
implicit none
! arguments
integer,    intent(in) :: fnum
complex(8), intent(in) :: zfmt(npcmtmax,natmtot),zfir(ngtot)
! local variables
integer np,ip,i,j,k
! allocatable arrays
real(8),    allocatable :: vpl(:,:)
complex(8), allocatable :: fp(:)
!-------------------------------------------------------------------------------

! Total number of plot points
np = np3d(1)*np3d(2)*np3d(3)
! Allocate local arrays
allocate(vpl(3,np),fp(np))
! Generate the 3D plotting points
call plotpt3d(vpl)
! Evaluate the functions at the grid points
call zfplotw90unk(np,vpl,zfmt(:,:),zfir(:),fp(:))

! Write functions to file
do i = 1,np3d(3)
  do j = 1,np3d(2)
    do k = 1,np3d(1)
      ip = ( j - 1 + ( i - 1 )*np3d(2) )*np3d(3) + k
      write(fnum,*) dreal(fp(ip)),dimag(fp(ip))
    end do
  end do
end do

deallocate(vpl,fp)

return

end subroutine
!EOC
