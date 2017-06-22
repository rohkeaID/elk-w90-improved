
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dpotks
use modmain
use modphonon
implicit none
! local variables
integer is,ias,nr,ir
! allocatable arrays
complex(8), allocatable :: zfmt(:,:)
! compute the exchange-correlation potential derivative
! (at this stage the density derivative is in spherical coordinates)
call dpotxc
! convert density derivative to spherical harmonics
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(zfmt,is,nr)
!$OMP DO
do ias=1,natmtot
  allocate(zfmt(lmmaxvr,nrmtmax))
  is=idxis(ias)
  nr=nrmt(is)
  call zcopy(lmmaxvr*nr,drhomt(:,:,ias),1,zfmt,1)
  call zfsht(nr,nrmtinr(is),zfmt,drhomt(:,:,ias))
  deallocate(zfmt)
end do
!$OMP END DO
!$OMP END PARALLEL
! generate the Coulomb potential derivative
call dpotcoul
! add to the Kohn-Sham potential derivative
do ias=1,natmtot
  is=idxis(ias)
  call zfmtadd(nrmt(is),nrmtinr(is),zone,dvclmt(:,:,ias),dvsmt(:,:,ias))
end do
call zaxpy(ngtot,zone,dvclir,1,dvsir,1)
! remove the gradient part of the potential derivative for displaced muffin-tin
call zfmtadd(nrmt(isph),nrmtinr(isph),zone,gvsmt,dvsmt(:,:,iasph))
! zero high (l,m) components on inner part of muffin-tin
do ias=1,natmtot
  is=idxis(ias)
  do ir=1,nrmtinr(is)
    dvsmt(lmmaxinr+1:lmmaxvr,ir,ias)=0.d0
  end do
end do
return
end subroutine

