
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: plot2d
! !INTERFACE:
subroutine plot2d(fnum,nf,rfmt,rfir)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   fnum : plot file number (in,integer)
!   nf   : number of functions (in,integer)
!   rfmt : real muffin-tin function (in,real(lmmaxvr,nrmtmax,natmtot,nf))
!   rfir : real intersitial function (in,real(ngtot,nf))
! !DESCRIPTION:
!   Produces a 2D plot of the real functions contained in arrays {\tt rfmt} and
!   {\tt rfir} on the parallelogram defined by the corner vertices in the global
!   array {\tt vclp2d}. See routine {\tt rfplot}.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: fnum,nf
real(8), intent(in) :: rfmt(lmmaxvr,nrmtmax,natmtot,nf),rfir(ngtot,nf)
! local variables
integer np,ip,i
real(8) vpnc(3)
! allocatable arrays
real(8), allocatable :: vpl(:,:),vppc(:,:),fp(:,:)
if ((nf.lt.1).or.(nf.gt.4)) then
  write(*,*)
  write(*,'("Error(plot2d): invalid number of functions : ",I8)') nf
  write(*,*)
  stop
end if
! allocate local arrays
np=np2d(1)*np2d(2)
allocate(vpl(3,np),vppc(2,np),fp(np,nf))
! generate the 2D plotting points
call plotpt2d(avec,ainv,vpnc,vpl,vppc)
! evaluate the functions at the grid points
do i=1,nf
  call rfplot(np,vpl,rfmt(:,:,:,i),rfir(:,i),fp(:,i))
end do
! write the functions to file
write(fnum,'(2I6," : grid size")') np2d(:)
do ip=1,np
  write(fnum,'(6G18.10)') vppc(1,ip),vppc(2,ip),(fp(ip,i),i=1,nf)
end do
deallocate(vpl,vppc,fp)
return
end subroutine
!EOC
