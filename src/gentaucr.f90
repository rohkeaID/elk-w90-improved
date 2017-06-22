
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gentaucr(taumt)
use modmain
implicit none
! arguments
real(8), intent(inout) :: taumt(lmmaxvr,nrmtmax,natmtot,nspinor)
! local variables
integer ist,ispn,jspn
integer is,ia,ias
integer nr,nri,ir,i,m
! allocatable arrays
complex(8), allocatable :: wfcr(:,:,:)
complex(8), allocatable :: gzfmt(:,:,:),zfmt(:,:)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(wfcr,gzfmt,zfmt) &
!$OMP PRIVATE(is,ia,nr,nri,ist,m) &
!$OMP PRIVATE(ispn,jspn,i,ir)
!$OMP DO
do ias=1,natmtot
  allocate(wfcr(lmmaxvr,nrmtmax,2))
  allocate(gzfmt(lmmaxvr,nrmtmax,3),zfmt(lmmaxvr,nrmtmax))
  is=idxis(ias)
  ia=idxia(ias)
  nr=nrmt(is)
  nri=nrmtinr(is)
  do ist=1,nstsp(is)
    if (spcore(ist,is)) then
      do m=-ksp(ist,is),ksp(ist,is)-1
! generate the core wavefunction in spherical harmonics (pass in m-1/2)
        call wavefcr(.true.,1,is,ia,ist,m,nrmtmax,wfcr)
        do ispn=1,2
          if (spinpol) then
            jspn=ispn
          else
            jspn=1
          end if
! compute the gradient
          call gradzfmt(nr,nri,rsp(:,is),wfcr(:,:,ispn),nrmtmax,gzfmt)
          do i=1,3
! convert gradient to spherical coordinates
            call zbsht(nr,nri,gzfmt(:,:,i),zfmt)
! add on inner part of muffin-tin
            do ir=1,nri
              taumt(1:lmmaxinr,ir,ias,jspn)=taumt(1:lmmaxinr,ir,ias,jspn) &
               +0.5d0*(dble(zfmt(1:lmmaxinr,ir))**2 &
               +aimag(zfmt(1:lmmaxinr,ir))**2)
            end do
! add on outer part of muffin-tin
            do ir=nri+1,nr
              taumt(:,ir,ias,jspn)=taumt(:,ir,ias,jspn) &
               +0.5d0*(dble(zfmt(:,ir))**2+aimag(zfmt(:,ir))**2)
            end do
          end do
        end do
      end do
    end if
  end do
  deallocate(wfcr,gzfmt,zfmt)
end do
!$OMP END DO
!$OMP END PARALLEL
return
end subroutine

