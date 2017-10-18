
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine xc_c_tb09
use modmain
implicit none
! local variables
integer is,ias,i
integer nr,nri,ir
real(8), parameter :: alpha=-0.012d0, beta=1.023d0
real(8) t1
! allocatable arrays
real(8), allocatable :: grfmt(:,:,:),grfir(:,:)
real(8), allocatable :: rfmt(:,:),rfir(:)
real(8), allocatable :: rfmt1(:),rfmt2(:,:)
! external functions
real(8) rfint
external rfint
! if Tran-Blaha constant has been read in return
if (tc_tb09) return
! compute the gradient of the density
allocate(grfmt(npmtmax,natmtot,3),grfir(ngtot,3))
call gradrf(rhomt,rhoir,grfmt,grfir)
allocate(rfmt(npmtmax,natmtot),rfmt1(npmtmax),rfmt2(npmtmax,3))
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
! convert muffin-tin density to spherical coordinates
  call rbsht(nr,nri,rhomt(:,ias),rfmt1)
! convert muffin-tin gradient to spherical coordinates
  do i=1,3
    call rbsht(nr,nri,grfmt(:,ias,i),rfmt2(:,i))
  end do
! integrand in muffin-tin
  do i=1,npmt(is)
    t1=sqrt(rfmt2(i,1)**2+rfmt2(i,2)**2+rfmt2(i,3)**2)
    rfmt1(i)=t1/rfmt1(i)
  end do
! convert to spherical harmonics
  call rfsht(nr,nri,rfmt1,rfmt(:,ias))
end do
deallocate(grfmt,rfmt1,rfmt2)
! integrand in interstitial
allocate(rfir(ngtot))
do ir=1,ngtot
  t1=sqrt(grfir(ir,1)**2+grfir(ir,2)**2+grfir(ir,3)**2)
  rfir(ir)=t1/rhoir(ir)
end do
! integrate over the unit cell
t1=rfint(rfmt,rfir)
! set the constant
c_tb09=alpha+beta*sqrt(abs(t1)/omega)
deallocate(grfir,rfmt,rfir)
return
end subroutine

