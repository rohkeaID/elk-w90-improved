
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rhomagsm
use modmain
implicit none
! local variables
integer is,ias,idm
if (msmooth.eq.0) return
! smooth the muffin-tin density
do ias=1,natmtot
  is=idxis(ias)
  call rfmtsm(msmooth,lmmaxvr,nrmt(is),lmmaxvr,rhomt(:,:,ias))
end do
! smooth the muffin-tin magnetisation
do idm=1,ndmag
  do ias=1,natmtot
    is=idxis(ias)
    call rfmtsm(msmooth,lmmaxvr,nrmt(is),lmmaxvr,magmt(:,:,ias,idm))
  end do
end do
return
end subroutine

