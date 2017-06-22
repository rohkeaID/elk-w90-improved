
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gengqrf(vqpc,igq0,vgqc,gqc,ylmgq,sfacgq)
use modmain
implicit none
! arguments
real(8), intent(in) :: vqpc(3)
integer, intent(out) :: igq0
real(8), intent(out) :: vgqc(3,ngrf),gqc(ngrf)
complex(8), intent(out) :: ylmgq(lmmaxvr,ngrf)
complex(8), intent(out) :: sfacgq(ngrf,natmtot)
! local variables
integer ig
real(8) tp(2)
do ig=1,ngrf
! determine the G+q-vectors
  vgqc(:,ig)=vgc(:,ig)+vqpc(:)
! G+q-vector length and (theta, phi) coordinates
  call sphcrd(vgqc(:,ig),gqc(ig),tp)
! spherical harmonics for G+q-vectors
  call genylm(lmaxvr,tp,ylmgq(:,ig))
end do
! structure factors for G+q
call gensfacgp(ngrf,vgqc,ngrf,sfacgq)
! find the shortest G+q-vector
call findigp0(ngrf,gqc,igq0)
return
end subroutine

