
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: gridsize
! !INTERFACE:
subroutine gridsize(avec,gmaxvr,ngridg,ngtot,intgv)
! !INPUT/OUTPUT PARAMETERS:
!   avec   : lattice vectors (in,real(3,3))
!   gmaxvr : G-vector cut-off (in,real)
!   ngridg : G-vector grid sizes (out,integer(3))
!   ngtot  : total number of G-vectors (out,integer)
!   intgv  : integer grid intervals for each direction (out,integer(2,3))
! !DESCRIPTION:
!   Finds the ${\bf G}$-vector grid which completely contains the vectors with
!   $G<G_{\rm max}$ and is compatible with the FFT routine. The optimal sizes
!   are given by
!   $$ n_i=\frac{G_{\rm max}|{\bf a}_i|}{\pi}+1, $$
!   where ${\bf a}_i$ is the $i$th lattice vector.
!
! !REVISION HISTORY:
!   Created July 2003 (JKD)
!   Removed modmain and added arguments, September 2012 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: avec(3,3),gmaxvr
integer, intent(out) :: ngridg(3),ngtot,intgv(2,3)
! local variables
real(8), parameter :: pi=3.1415926535897932385d0
! find optimal grid size for potential and density
ngridg(:)=int(gmaxvr*sqrt(avec(1,:)**2+avec(2,:)**2+avec(3,:)**2)/pi)+1
! find next largest FFT-compatible grid size
call nfftifc(ngridg(1))
call nfftifc(ngridg(2))
call nfftifc(ngridg(3))
if ((ngridg(1).le.0).or.(ngridg(2).le.0).or.(ngridg(3).le.0)) then
  write(*,*)
  write(*,'("Error(gridsize): invalid ngridg : ",3I8)') ngridg
  write(*,*)
  stop
end if
! total number of points in grid
ngtot=ngridg(1)*ngridg(2)*ngridg(3)
! determine integer ranges for grid
intgv(1,:)=ngridg(:)/2-ngridg(:)+1
intgv(2,:)=ngridg(:)/2
return
end subroutine
!EOC

