
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine init3
use modmain
use modvars
implicit none
! local variables
integer ig,iw
real(8) w1,w2,t1,t2

!-------------------------------------------------------------!
!     response function and perturbation theory variables     !
!-------------------------------------------------------------!
! G-vectors for response functions
ngrf=1
do ig=ngvec,1,-1
  if (gc(ig).lt.gmaxrf) then
    ngrf=ig
    exit
  end if
end do
! frequencies for reponse functions
nwrf=1
if (allocated(wrf)) deallocate(wrf)
if ((task.eq.188).or.(task.eq.320).or.(task.eq.330).or.(task.eq.331)) then
  nwrf=nwplot
  allocate(wrf(nwrf))
  w1=wplot(1)
  w2=max(wplot(2),w1)
  t1=(w2-w1)/dble(nwplot)
  do iw=1,nwplot
    t2=w1+t1*dble(iw-1)
    wrf(iw)=cmplx(t2,swidth,8)
  end do
else
  nwrf=1
  allocate(wrf(nwrf))
  wrf(1)=cmplx(0.d0,swidth,8)
end if
! write to VARIABLES.OUT
call writevars('gmaxrf',rv=gmaxrf)
call writevars('ngrf',iv=ngrf)
call writevars('nwrf',iv=nwrf)
call writevars('wrf',nv=nwrf,zva=wrf)

return
end subroutine

