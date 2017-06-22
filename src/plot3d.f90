
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: plot3d
! !INTERFACE:
subroutine plot3d(fnum,nf,rfmt,rfir)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   fnum : plot file number (in,integer)
!   nf   : number of functions (in,integer)
!   rfmt : real muffin-tin function (in,real(lmmaxvr,nrmtmax,natmtot,nf))
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
integer, intent(in) :: fnum,nf
real(8), intent(in) :: rfmt(lmmaxvr,nrmtmax,natmtot,nf),rfir(ngtot,nf)
! local variables
integer np,ip,i
real(8) v1(3)
! allocatable arrays
real(8), allocatable :: vpl(:,:),fp(:,:)
if ((nf.lt.1).or.(nf.gt.4)) then
  write(*,*)
  write(*,'("Error(plot3d): invalid number of functions : ",I8)') nf
  write(*,*)
  stop
end if
! total number of plot points
np=np3d(1)*np3d(2)*np3d(3)
! allocate local arrays
allocate(vpl(3,np),fp(np,nf))
! generate the 3D plotting points
call plotpt3d(vpl)
! evaluate the functions at the grid points
do i=1,nf
  call rfplot(np,vpl,rfmt(:,:,:,i),rfir(:,i),fp(:,i))
end do
! write functions to file
write(fnum,'(3I6," : grid size")') np3d(:)
do ip=1,np
  call r3mv(avec,vpl(:,ip),v1)
  write(fnum,'(7G18.10)') v1(:),(fp(ip,i),i=1,nf)
end do
deallocate(vpl,fp)
return
end subroutine
!EOC
