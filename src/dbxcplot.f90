
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dbxcplot
use modmain
implicit none
! local variables
integer idm,is,ias,np
! allocatable arrays
real(8), allocatable :: rvfmt(:,:,:),rvfir(:,:)
real(8), allocatable :: rfmt(:,:),rfir(:)
real(8), allocatable :: grfmt(:,:,:),grfir(:,:)
! initialise universal variables
call init0
if (.not.spinpol) then
  write(*,*)
  write(*,'("Error(dbxcplot): spin-unpolarised magnetic field is zero")')
  write(*,*)
  stop
end if
! read magnetisation from file
call readstate
allocate(rvfmt(npmtmax,natmtot,3),rvfir(ngtot,3))
allocate(rfmt(npmtmax,natmtot),rfir(ngtot))
allocate(grfmt(npmtmax,natmtot,3),grfir(ngtot,3))
if (ncmag) then
! non-collinear
  rvfmt(:,:,:)=bxcmt(:,:,:)
  rvfir(:,:)=bxcir(:,:)
else
! collinear
  rvfmt(:,:,1:2)=0.d0
  rvfir(:,1:2)=0.d0
  rvfmt(:,:,3)=bxcmt(:,:,1)
  rvfir(:,3)=bxcir(:,1)
end if
rfmt(:,:)=0.d0
rfir(:)=0.d0
do idm=1,3
  call gradrf(rvfmt(:,:,idm),rvfir(:,idm),grfmt,grfir)
  do ias=1,natmtot
    is=idxis(ias)
    np=npmt(is)
    rfmt(1:np,ias)=rfmt(1:np,ias)+grfmt(1:np,ias,idm)
  end do
  rfir(:)=rfir(:)+grfir(:,idm)
end do
select case(task)
case(91)
  open(50,file='DBXC1D.OUT',form='FORMATTED')
  open(51,file='DBXCLINES.OUT',form='FORMATTED')
  call plot1d(50,51,1,rfmt,rfir)
  close(50)
  close(51)
  write(*,*)
  write(*,'("Info(dbxcplot):")')
  write(*,'(" 1D divergence of exchange-correlation field written to &
   &DBXC1D.OUT")')
  write(*,'(" vertex location lines written to DBXCLINES.OUT")')
case(92)
  open(50,file='DBXC2D.OUT',form='FORMATTED')
  call plot2d(50,1,rfmt,rfir)
  close(50)
  write(*,'("Info(dbxcplot):")')
  write(*,'(" 2D divergence of exchange-correlation field written to &
   &DBXC2D.OUT")')
case(93)
  open(50,file='DBXC3D.OUT',form='FORMATTED')
  call plot3d(50,1,rfmt,rfir)
  close(50)
  write(*,'("Info(dbxcplot):")')
  write(*,'(" 3D divergence of exchange-correlation field written to &
   &DBXC3D.OUT")')
end select
deallocate(rvfmt,rvfir,rfmt,rfir,grfmt,grfir)
return
end subroutine

