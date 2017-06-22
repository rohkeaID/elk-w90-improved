
! Copyright (C) 2014 K. Krieger, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writetddft
use modmain
use modtddft
use moddftu
implicit none
! local variables
integer ias
real(8) t1
character(256) fext
! allocatable arrays
real(8), allocatable :: rvfmt(:,:,:,:),rvfir(:,:)
! file extension
write(fext,'("_TS",I8.8,".OUT")') itimes
! delete all files at first time-step
if (itimes.le.1) then
  open(50,file='CHARGEMT_TD.OUT')
  close(50,status='DELETE')
  open(50,file='CHARGEIR_TD.OUT')
  close(50,status='DELETE')
  open(50,file='MOMENT_TD.OUT')
  close(50,status='DELETE')
  open(50,file='MOMENTM_TD.OUT')
  close(50,status='DELETE')
  open(50,file='CURRENT_TD.OUT')
  close(50,status='DELETE')
  open(50,file='CURRENTM_TD.OUT')
  close(50,status='DELETE')
end if
! muffin-tin charges
open(50,file='CHARGEMT_TD.OUT',action='WRITE',form='FORMATTED', &
 position='APPEND')
write(50,'(G18.10)',advance='NO') times(itimes)
do ias=1,natmtot
  write(50,'(G18.10)',advance='NO') chgmt(ias)
end do
write(50,*)
close(50)
! interstitial charge
open(50,file='CHARGEIR_TD.OUT',action='WRITE',form='FORMATTED', &
 position='APPEND')
write(50,'(2G18.10)') times(itimes),chgir
close(50)
! spin moment
open(50,file='MOMENT_TD.OUT',action='WRITE',form='FORMATTED',position='APPEND')
write(50,'(4G18.10)') times(itimes),momtot(1:ndmag)
close(50)
open(50,file='MOMENTM_TD.OUT',action='WRITE',form='FORMATTED',position='APPEND')
write(50,'(2G18.10)') times(itimes),momtotm
close(50)
! total current
open(50,file='CURRENT_TD.OUT',action='WRITE',form='FORMATTED',position='APPEND')
write(50,'(4G18.10)') times(itimes),curtot(:)
close(50)
! total current magnitude
open(50,file='CURRENTM_TD.OUT',action='WRITE',form='FORMATTED', &
 position='APPEND')
t1=sqrt(sum(curtot(:)**2))
write(50,'(4G18.10)') times(itimes),t1
close(50)
! write optional quantities
if (ntswrite.le.0) goto 10
if (mod(itimes-1,ntswrite).ne.0) goto 10
! charge density in 1D
if (tdrho1d) then
  open(50,file='RHO1D'//trim(fext),action='WRITE',form='FORMATTED')
  open(51,file='RHOLINES'//trim(fext),action='WRITE',form='FORMATTED')
  call plot1d(50,51,1,rhomt,rhoir)
  close(50)
  close(51)
end if
! charge density in 2D
if (tdrho2d) then
  open(50,file='RHO2D'//trim(fext),action='WRITE',form='FORMATTED')
  call plot2d(50,1,rhomt,rhoir)
  close(50)
end if
! charge density in 3D
if (tdrho3d) then
  open(50,file='RHO3D'//trim(fext),action='WRITE',form='FORMATTED')
  call plot3d(50,1,rhomt,rhoir)
  close(50)
end if
! magnetisation in 2D or 3D
if (tdmag2d.or.tdmag3d) then
  allocate(rvfmt(lmmaxvr,nrmtmax,natmtot,3),rvfir(ngtot,3))
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
  if (tdmag2d) then
! project the magnetisation onto the plotting plane
    call proj2d(rvfmt,rvfir)
    open(50,file='MAG2D'//trim(fext),action='WRITE',form='FORMATTED')
    call plot2d(50,3,rvfmt,rvfir)
    close(50)
  else
    open(50,file='MAG3D'//trim(fext),action='WRITE',form='FORMATTED')
    call plot3d(50,3,rvfmt,rvfir)
    close(50)
  end if
  deallocate(rvfmt,rvfir)
end if
! calculate and write tensor moments
if (dftu.ne.0) then
  if (tmwrite) then
    open(50,file='TMDFTU'//trim(fext),action='WRITE',form='FORMATTED')
    if (spinorb) then
      call writetm3du(50)
    else
      call writetm2du(50)
    end if
    close(50)
  end if
end if
10 continue
return
end subroutine

