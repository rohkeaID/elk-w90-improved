
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: plot3d
! !INTERFACE:
subroutine plotUNK(fnum,zfmt,zfir)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   fnum : plot file number (in,integer)
!   nf   : number of functions (in,integer)
!   zfmt : real muffin-tin function (in,real(npmtmax,natmtot,nf))
!   zfir : real intersitial function (in,real(ngtot,nf))
! !DESCRIPTION:
!   Produces a 3D plot of the real functions contained in arrays {\tt zfmt} and
!   {\tt zfir} in the parallelepiped defined by the corner vertices in the
!   global array {\tt vclp3d}. See routine {\tt rfarray}.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!   Modified, October 2008 (F. Bultmark, F. Cricchio, L. Nordstrom)
!EOP
!BOC
implicit none
! arguments
integer,    intent(in) :: fnum
complex(8), intent(in) :: zfmt(npcmtmax,natmtot),zfir(ngtot)
! local variables
integer np,ip,i,j,k
! allocatable arrays
real(8), allocatable :: vpl(:,:),vpc(:,:),fp_real(:),fp_imag(:)
complex(8), allocatable :: fp(:)
! total number of plot points
np=np3d(1)*np3d(2)*np3d(3)
! allocate local arrays
allocate(vpl(3,np),vpc(3,np),fp_real(np),fp_imag(np),fp(np))
! generate the 3D plotting points
call plotpt3d(vpl)
! evaluate the functions at the grid points
call zfplotUNK(np,vpl,zfmt(:,:),zfir(:),fp(:))

! write functions to file
do i=1,np3d(3)
  do j=1,np3d(2)
    do k=1,np3d(1)
      ip = (j-1+(i-1)*np3d(2))*np3d(3) + k
      write(fnum,*) dreal(fp(ip)),dimag(fp(ip))
      !write(fnum,'(2F7.3)') dreal(fp(ip)),dimag(fp(ip))
    end do
  end do
end do

deallocate(vpl,fp_real,fp_imag,fp)
return
end subroutine
!EOC
