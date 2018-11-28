
! Copyright (C) 2017-18 Arsenii Gerasimov, Yaroslav Kvashnin and Lars Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: plotUNK
! !INTERFACE:
subroutine plotUNK(fnum,zfmt,zfir)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   fnum : plot file number (in,integer)
!   nf   : number of functions (in,integer)
!   zfmt : complex muffin-tin  function (in,complex(npmtmax,natmtot,nf))
!   zfir : complex intersitial function (in,complex(ngtot,nf))
! !DESCRIPTION:
!
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
real(8),    allocatable :: vpl(:,:),vpc(:,:),fp_real(:),fp_imag(:)
complex(8), allocatable :: fp(:)
!-------------------------------------------------------------------------------
! Total number of plot points
np = np3d(1)*np3d(2)*np3d(3)
! Allocate local arrays
allocate(vpl(3,np),vpc(3,np),fp_real(np),fp_imag(np),fp(np))
! Generate the 3D plotting points
call plotpt3d(vpl)
! Evaluate the functions at the grid points
call zfplotUNK(np,vpl,zfmt(:,:),zfir(:),fp(:))

! Write functions to file
do i = 1,np3d(3)
  do j = 1,np3d(2)
    do k = 1,np3d(1)
      ip = ( j - 1 + (i-1)*np3d(2) )*np3d(3) + k
      write(fnum,*) dreal(fp(ip)),dimag(fp(ip))
    end do
  end do
end do

deallocate(vpl,fp_real,fp_imag,fp)
return
end subroutine
!EOC
