
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dpotxc
use modmain
use modphonon
implicit none
! local variables
integer idm,is,ias,ir
integer nr,nri,nrc,nrci
! allocatable arrays
real(8), allocatable :: fxcmt(:,:,:,:,:),fxcir(:,:,:)
complex(8), allocatable :: dvmt(:,:),dbmt(:,:,:)
! compute the exchange-correlation kernel
if (spinpol) then
  allocate(fxcmt(lmmaxvr,nrmtmax,natmtot,4,4))
  allocate(fxcir(ngtot,4,4))
  call genspfxcr(.false.,fxcmt,fxcir)
else
  allocate(fxcmt(lmmaxvr,nrmtmax,natmtot,1,1))
  allocate(fxcir(ngtot,1,1))
  call genfxcr(.false.,fxcmt,fxcir)
end if
allocate(dvmt(lmmaxvr,nrmtmax))
if (spinpol) allocate(dbmt(lmmaxvr,nrmtmax,3))
!---------------------------------------!
!     muffin-tin potential and field    !
!---------------------------------------!
! note: muffin-tin functions are in spherical coordinates
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmtinr(is)
  nrc=nrcmt(is)
  nrci=nrcmtinr(is)
! charge-charge contribution to potential derivative
  do ir=1,nr
    dvmt(:,ir)=fxcmt(:,ir,ias,1,1)*drhomt(:,ir,ias)
  end do
! spin-polarised case
  if (spinpol) then
    if (ncmag) then
! non-collinear
      do ir=1,nr
! add charge-spin contribution to potential derivative
        dvmt(:,ir)=dvmt(:,ir) &
         +fxcmt(:,ir,ias,1,2)*dmagmt(:,ir,ias,1) &
         +fxcmt(:,ir,ias,1,3)*dmagmt(:,ir,ias,2) &
         +fxcmt(:,ir,ias,1,4)*dmagmt(:,ir,ias,3)
! spin-charge contribution to B-field derivative
        dbmt(:,ir,1)=fxcmt(:,ir,ias,1,2)*drhomt(:,ir,ias)
        dbmt(:,ir,2)=fxcmt(:,ir,ias,1,3)*drhomt(:,ir,ias)
        dbmt(:,ir,3)=fxcmt(:,ir,ias,1,4)*drhomt(:,ir,ias)
! add spin-spin contribution to B-field derivative
! (note: fxc is stored as an upper triangular matrix)
        dbmt(:,ir,1)=dbmt(:,ir,1) &
         +fxcmt(:,ir,ias,2,2)*dmagmt(:,ir,ias,1) &
         +fxcmt(:,ir,ias,2,3)*dmagmt(:,ir,ias,2) &
         +fxcmt(:,ir,ias,2,4)*dmagmt(:,ir,ias,3)
        dbmt(:,ir,2)=dbmt(:,ir,2) &
         +fxcmt(:,ir,ias,2,3)*dmagmt(:,ir,ias,1) &
         +fxcmt(:,ir,ias,3,3)*dmagmt(:,ir,ias,2) &
         +fxcmt(:,ir,ias,3,4)*dmagmt(:,ir,ias,3)
        dbmt(:,ir,3)=dbmt(:,ir,3) &
         +fxcmt(:,ir,ias,2,4)*dmagmt(:,ir,ias,1) &
         +fxcmt(:,ir,ias,3,4)*dmagmt(:,ir,ias,2) &
         +fxcmt(:,ir,ias,4,4)*dmagmt(:,ir,ias,3)
      end do
    else
! collinear
      do ir=1,nr
! add charge-spin contribution to potential derivative
        dvmt(:,ir)=dvmt(:,ir)+fxcmt(:,ir,ias,1,4)*dmagmt(:,ir,ias,1)
! spin-charge contribution to B-field derivative
        dbmt(:,ir,1)=fxcmt(:,ir,ias,1,4)*drhomt(:,ir,ias)
! add spin-spin contribution to B-field derivative
        dbmt(:,ir,1)=dbmt(:,ir,1)+fxcmt(:,ir,ias,4,4)*dmagmt(:,ir,ias,1)
      end do
    end if
  end if
! convert potential derivative to spherical harmonics
  call zfsht(nr,nri,dvmt,dvsmt(:,:,ias))
! convert magnetic field derivative to spherical harmonics
  do idm=1,ndmag
    call zfsht(nrc,nrci,dbmt(:,:,idm),dbsmt(:,:,ias,idm))
  end do
end do
!------------------------------------------!
!     interstitial potential and field     !
!------------------------------------------!
! charge-charge contribution to potential derivative
do ir=1,ngtot
  dvsir(ir)=fxcir(ir,1,1)*drhoir(ir)
end do
! spin-polarised case
if (spinpol) then
  if (ncmag) then
! non-collinear
    do ir=1,ngtot
! add charge-spin contribution to potential derivative
      dvsir(ir)=dvsir(ir) &
       +fxcir(ir,1,2)*dmagir(ir,1) &
       +fxcir(ir,1,3)*dmagir(ir,2) &
       +fxcir(ir,1,4)*dmagir(ir,3)
! spin-charge contribution to B-field derivative
      dbsir(ir,1)=fxcir(ir,1,2)*drhoir(ir)
      dbsir(ir,2)=fxcir(ir,1,3)*drhoir(ir)
      dbsir(ir,3)=fxcir(ir,1,4)*drhoir(ir)
! add spin-spin contribution to B-field derivative
      dbsir(ir,1)=dbsir(ir,1) &
       +fxcir(ir,2,2)*dmagir(ir,1) &
       +fxcir(ir,2,3)*dmagir(ir,2) &
       +fxcir(ir,2,4)*dmagir(ir,3)
      dbsir(ir,2)=dbsir(ir,2) &
       +fxcir(ir,2,3)*dmagir(ir,1) &
       +fxcir(ir,3,3)*dmagir(ir,2) &
       +fxcir(ir,3,4)*dmagir(ir,3)
      dbsir(ir,3)=dbsir(ir,3) &
       +fxcir(ir,2,4)*dmagir(ir,1) &
       +fxcir(ir,3,4)*dmagir(ir,2) &
       +fxcir(ir,4,4)*dmagir(ir,3)
    end do
  else
! collinear
    do ir=1,ngtot
! add charge-spin contribution to potential derivative
      dvsir(ir)=dvsir(ir)+fxcir(ir,1,4)*dmagir(ir,1)
! spin-charge contribution to B-field derivative
      dbsir(ir,1)=fxcir(ir,1,4)*drhoir(ir)
! add spin-spin contribution to B-field derivative
      dbsir(ir,1)=dbsir(ir,1)+fxcir(ir,4,4)*dmagir(ir,1)
    end do
  end if
end if
deallocate(fxcmt,fxcir,dvmt)
if (spinpol) deallocate(dbmt)
return
end subroutine

