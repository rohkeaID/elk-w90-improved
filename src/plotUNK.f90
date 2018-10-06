
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: plot3d
! !INTERFACE:
subroutine plotUNK(fnum,rfmt,rfir)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   fnum : plot file number (in,integer)
!   nf   : number of functions (in,integer)
!   rfmt : real muffin-tin function (in,real(npmtmax,natmtot,nf))
!   rfir : real intersitial function (in,real(ngtot,nf))
! !DESCRIPTION:
!   Produces a 3D plot of the real functions contained in arrays {\tt rfmt} and
!   {\tt rfir} in the parallelepiped defined by the corner vertices in the
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
complex(8), intent(in) :: rfmt(npmtmax,natmtot),rfir(ngtot)
! local variables
integer np,ip,i,j,k
! allocatable arrays
real(8), allocatable :: vpl(:,:),fp_real(:),fp_imag(:)
complex(8), allocatable :: fp(:)
! total number of plot points
np=np3d(1)*np3d(2)*np3d(3)
! allocate local arrays
allocate(vpl(3,np),fp_real(np),fp_imag(np),fp(np))
! generate the 3D plotting points
call plotpt3d(vpl)
! evaluate the functions at the grid points
call rfplot(np,vpl,dreal(rfmt(:,:)),dreal(rfir(:)),fp_real(:))
call rfplot(np,vpl,dimag(rfmt(:,:)),dimag(rfir(:)),fp_imag(:))
!call rfplotUNK(np,vpl,rfmt(:,:),rfir(:),fp_imag(:))

! write functions to file
!write(fnum,'(3I6," : grid size")') np3d(:)
!do ip=1,np
!  write(fnum,'(2F7.3)') fp_real(ip),fp_imag(ip)
!end do
!np3d(1)=3
!np3d(2)=4
!np3d(3)=5
do k=1,np3d(3)
  do j=1,np3d(2)
    do i=1,np3d(1)
      ip = (j-1+(i-1)*np3d(2))*np3d(3) + k
      !write(fnum,'(2F7.3)') fp_real(ip),fp_imag(ip)
      write(fnum,*) fp_real(ip),fp_imag(ip)
    end do
  end do
end do
deallocate(vpl,fp_real,fp_imag,fp)
return
end subroutine
!EOC
