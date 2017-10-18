
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gentaucr(taumt)
use modmain
implicit none
! arguments
real(8), intent(inout) :: taumt(npmtmax,natmtot,nspinor)
! local variables
integer ist,ispn,jspn
integer is,ia,ias
integer nr,nri,np,i,m
! allocatable arrays
complex(8), allocatable :: wfcr(:,:)
complex(8), allocatable :: gzfmt(:,:),zfmt(:)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(wfcr,gzfmt,zfmt) &
!$OMP PRIVATE(is,ia,nr,nri,np) &
!$OMP PRIVATE(ist,m,ispn,jspn,i)
!$OMP DO
do ias=1,natmtot
  allocate(wfcr(npmtmax,2))
  allocate(gzfmt(npmtmax,3),zfmt(npmtmax))
  is=idxis(ias)
  ia=idxia(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  np=npmt(is)
  do ist=1,nstsp(is)
    if (spcore(ist,is)) then
      do m=-ksp(ist,is),ksp(ist,is)-1
! generate the core wavefunction in spherical harmonics (pass in m-1/2)
        call wavefcr(.true.,1,is,ia,ist,m,npmtmax,wfcr)
        do ispn=1,2
          if (spinpol) then
            jspn=ispn
          else
            jspn=1
          end if
! compute the gradient
          call gradzfmt(nr,nri,rsp(:,is),wfcr(:,ispn),npmtmax,gzfmt)
          do i=1,3
! convert gradient to spherical coordinates
            call zbsht(nr,nri,gzfmt(:,i),zfmt)
! add to total in muffin-tin
            taumt(1:np,ias,jspn)=taumt(1:np,ias,jspn) &
             +0.5d0*(dble(zfmt(1:np))**2+aimag(zfmt(1:np))**2)
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

