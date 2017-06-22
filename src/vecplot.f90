
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: vecplot
! !INTERFACE:
subroutine vecplot
! !DESCRIPTION:
!   Outputs a 2D or 3D vector field for plotting. The vector field can be the
!   magnetisation vector field, ${\bf m}$; the exchange-correlation magnetic
!   field, ${\bf B}_{\rm xc}$; or the electric field
!   ${\bf E}\equiv-\nabla V_{\rm C}$. The magnetisation is obtained from the
!   spin density matrix, $\rho_{\alpha\beta}$, by solving
!   $$ \rho_{\alpha\beta}({\bf r})=\frac{1}{2}\left(n({\bf r})
!    \delta_{\alpha\beta}+\sigma\cdot {\bf m}({\bf r})\right), $$
!   where $n\equiv\tr\rho_{\alpha\beta}$ is the total density. In the case of 2D
!   plots, the magnetisation vectors are still 3D, but are in the coordinate
!   system of the plane.
!
! !REVISION HISTORY:
!   Created August 2004 (JKD)
!   Included electric field plots, August 2006 (JKD)
!EOP
!BOC
use modmain
implicit none
! allocatable arrays
real(8), allocatable :: rvfmt(:,:,:,:),rvfir(:,:)
! initialise universal variables
call init0
if ((task.eq.72).or.(task.eq.73).or.(task.eq.82).or.(task.eq.83)) then
  if (.not.spinpol) then
    write(*,*)
    write(*,'("Error(vecplot): spin-unpolarised magnetisation/field is zero")')
    write(*,*)
    stop
  end if
end if
! read magnetisation from file
call readstate
allocate(rvfmt(lmmaxvr,nrmtmax,natmtot,3),rvfir(ngtot,3))
select case(task)
case(72,73)
! magnetisation
  if (ncmag) then
! non-collinear
    rvfmt(:,:,:,:)=magmt(:,:,:,:)
    rvfir(:,:)=magir(:,:)
  else
! collinear
    rvfmt(:,:,:,1:2)=0.d0
    rvfir(:,1:2)=0.d0
    rvfmt(:,:,:,3)=magmt(:,:,:,1)
    rvfir(:,3)=magir(:,1)
  end if
case(82,83)
! exchange-correlation magnetic field
  if (ncmag) then
! non-collinear
    rvfmt(:,:,:,:)=bxcmt(:,:,:,:)
    rvfir(:,:)=bxcir(:,:)
  else
! collinear
    rvfmt(:,:,:,1:2)=0.d0
    rvfir(:,1:2)=0.d0
    rvfmt(:,:,:,3)=bxcmt(:,:,:,1)
    rvfir(:,3)=bxcir(:,1)
  end if
case(142,143)
! electric field
  call gradrf(vclmt,vclir,rvfmt,rvfir)
! use the negative of the gradient
  rvfmt(:,:,:,:)=-rvfmt(:,:,:,:)
  rvfir(:,:)=-rvfir(:,:)
case(152,153)
  if (ndmag.lt.3) then
    write(*,*)
    write(*,'("Error(vecplot): collinear m(r)xB(r) is zero")')
    write(*,*)
    stop
  end if
  call rvfcross(magmt,magir,bxcmt,bxcir,rvfmt,rvfir)
end select
select case(task)
case(72,82,142,152)
! project the magnetisation/field onto the plotting plane
  call proj2d(rvfmt,rvfir)
  if (task.eq.72) then
    open(50,file='MAG2D.OUT',action='WRITE',form='FORMATTED')
  else if (task.eq.82) then
    open(50,file='BXC2D.OUT',action='WRITE',form='FORMATTED')
  else if (task.eq.142) then
    open(50,file='EF2D.OUT',action='WRITE',form='FORMATTED')
  else
    open(50,file='MCBXC2D.OUT',action='WRITE',form='FORMATTED')
  end if
  call plot2d(50,3,rvfmt,rvfir)
  close(50)
  write(*,*)
  write(*,'("Info(vecplot):")')
  if (task.eq.72) then
    write(*,'(" 2D magnetisation density written to MAG2D.OUT")')
  else if (task.eq.82) then
    write(*,'(" 2D exchange-correlation field written to BXC2D.OUT")')
  else if (task.eq.142) then
    write(*,'(" 2D electric field written to EF2D.OUT")')
  else
    write(*,'(" 2D m(r) x B_xc(r) written to MCBXC2D.OUT")')
  end if
case(73,83,143,153)
  if (task.eq.73) then
    open(50,file='MAG3D.OUT',action='WRITE',form='FORMATTED')
  else if (task.eq.83) then
    open(50,file='BXC3D.OUT',action='WRITE',form='FORMATTED')
  else if (task.eq.143) then
    open(50,file='EF3D.OUT',action='WRITE',form='FORMATTED')
  else
    open(50,file='MCBXC3D.OUT',action='WRITE',form='FORMATTED')
  end if
  call plot3d(50,3,rvfmt,rvfir)
  close(50)
  write(*,*)
  write(*,'("Info(vecplot):")')
  if (task.eq.73) then
    write(*,'(" 3D magnetisation density written to MAG3D.OUT")')
  else if (task.eq.83) then
    write(*,'(" 3D exchange-correlation field written to BXC3D.OUT")')
  else if (task.eq.143) then
    write(*,'(" 3D electric field written to EF3D.OUT")')
  else
    write(*,'(" 3D m(r) x B_xc(r) written to MCBXC3D.OUT")')
  end if
end select
deallocate(rvfmt,rvfir)
return
end subroutine
!EOC

