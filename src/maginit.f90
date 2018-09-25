
! Copyright (C) 2016 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine maginit
use modmain
implicit none
! local variables
integer idm,is,ia,ias,np
! magnetisation as fraction of density
real(8), parameter :: fmr=0.15d0
real(8) v(3),t1
if (.not.spinpol) return
! initialise muffin-tin magnetisation
do ias=1,natmtot
  is=idxis(ias)
  ia=idxia(ias)
  np=npmt(is)
  v(:)=bfcmt(:,ia,is)+bfieldc(:)
  t1=sqrt(v(1)**2+v(2)**2+v(3)**2)
  if (t1.gt.1.d-8) then
    t1=-fmr/t1
    v(:)=t1*v(:)
    if (.not.ncmag) v(1)=v(3)
    do idm=1,ndmag
      t1=v(idm)
      magmt(1:np,ias,idm)=t1*rhomt(1:np,ias)
    end do
  else
    magmt(1:np,ias,:)=0.d0
  end if
end do
! initialise interstitial magnetisation
v(:)=bfieldc(:)
t1=sqrt(v(1)**2+v(2)**2+v(3)**2)
if (t1.gt.1.d-8) then
  t1=-fmr/t1
  v(:)=t1*v(:)
  if (.not.ncmag) v(1)=v(3)
  do idm=1,ndmag
    t1=v(idm)
    magir(:,idm)=t1*rhoir(:)
  end do
else
  magir(:,:)=0.d0
end if
return
end subroutine

