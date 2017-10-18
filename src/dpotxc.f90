
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dpotxc
use modmain
use modphonon
implicit none
! local variables
integer idm,is,ias
integer nr,nri,nrc,nrci
integer ir,np
! allocatable arrays
real(8), allocatable :: fxcmt(:,:,:,:),fxcir(:,:,:)
complex(8), allocatable :: dvmt(:),dbmt(:,:),zfmt(:)
! compute the exchange-correlation kernel
if (spinpol) then
  allocate(fxcmt(npmtmax,natmtot,4,4),fxcir(ngtot,4,4))
  call genspfxcr(.false.,fxcmt,fxcir)
else
  allocate(fxcmt(npmtmax,natmtot,1,1),fxcir(ngtot,1,1))
  call genfxcr(.false.,fxcmt,fxcir)
end if
allocate(dvmt(npmtmax))
if (spinpol) allocate(dbmt(npmtmax,3),zfmt(npcmtmax))
!---------------------------------------!
!     muffin-tin potential and field    !
!---------------------------------------!
! note: muffin-tin functions are in spherical coordinates
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  np=npmt(is)
! charge-charge contribution to potential derivative
  dvmt(1:np)=fxcmt(1:np,ias,1,1)*drhomt(1:np,ias)
! spin-polarised case
  if (spinpol) then
    if (ncmag) then
! non-collinear
! add charge-spin contribution to potential derivative
      dvmt(1:np)=dvmt(1:np) &
       +fxcmt(1:np,ias,1,2)*dmagmt(1:np,ias,1) &
       +fxcmt(1:np,ias,1,3)*dmagmt(1:np,ias,2) &
       +fxcmt(1:np,ias,1,4)*dmagmt(1:np,ias,3)
! spin-charge contribution to B-field derivative
      dbmt(1:np,1)=fxcmt(1:np,ias,1,2)*drhomt(1:np,ias)
      dbmt(1:np,2)=fxcmt(1:np,ias,1,3)*drhomt(1:np,ias)
      dbmt(1:np,3)=fxcmt(1:np,ias,1,4)*drhomt(1:np,ias)
! add spin-spin contribution to B-field derivative
! (note: fxc is stored as an upper triangular matrix)
      dbmt(1:np,1)=dbmt(1:np,1) &
       +fxcmt(1:np,ias,2,2)*dmagmt(1:np,ias,1) &
       +fxcmt(1:np,ias,2,3)*dmagmt(1:np,ias,2) &
       +fxcmt(1:np,ias,2,4)*dmagmt(1:np,ias,3)
      dbmt(1:np,2)=dbmt(1:np,2) &
       +fxcmt(1:np,ias,2,3)*dmagmt(1:np,ias,1) &
       +fxcmt(1:np,ias,3,3)*dmagmt(1:np,ias,2) &
       +fxcmt(1:np,ias,3,4)*dmagmt(1:np,ias,3)
      dbmt(1:np,3)=dbmt(1:np,3) &
       +fxcmt(1:np,ias,2,4)*dmagmt(1:np,ias,1) &
       +fxcmt(1:np,ias,3,4)*dmagmt(1:np,ias,2) &
       +fxcmt(1:np,ias,4,4)*dmagmt(1:np,ias,3)
    else
! collinear
! add charge-spin contribution to potential derivative
      dvmt(1:np)=dvmt(1:np)+fxcmt(1:np,ias,1,4)*dmagmt(1:np,ias,1)
! spin-charge contribution to B-field derivative
      dbmt(1:np,1)=fxcmt(1:np,ias,1,4)*drhomt(1:np,ias)
! add spin-spin contribution to B-field derivative
      dbmt(1:np,1)=dbmt(1:np,1)+fxcmt(1:np,ias,4,4)*dmagmt(1:np,ias,1)
    end if
  end if
! convert potential derivative to spherical harmonics
  call zfsht(nr,nri,dvmt,dvsmt(:,ias))
! convert magnetic field derivative to spherical harmonics on coarse mesh
  do idm=1,ndmag
    call zfmtftoc(nr,nri,dbmt(:,idm),zfmt)
    call zfsht(nrc,nrci,zfmt,dbsmt(:,ias,idm))
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
if (spinpol) deallocate(dbmt,zfmt)
return
end subroutine

