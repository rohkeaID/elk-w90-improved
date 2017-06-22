
! Copyright (C) 2014 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gengqvec(iq)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: iq
! local variables
integer is,ig
real(8) tp(2)
! loop over G-vectors
do ig=1,ngtot
! G+q-vector in Cartesian coordinates
  vgqc(:,ig)=vgc(:,ig)+vqc(:,iq)
! G+q-vector length and (theta, phi) coordinates
  call sphcrd(vgqc(:,ig),gqc(ig),tp)
! spherical harmonics for G+q-vectors
  call genylm(lmaxvr,tp,ylmgq(:,ig))
end do
! compute the spherical Bessel functions j_l(|G+q|R_mt)
call genjlgpr(lnpsd,gqc,jlgqr)
! structure factors for G+q
call gensfacgp(ngvec,vgqc,ngvec,sfacgq)
! generate the smooth step function form factors for G+q
do is=1,nspecies
  call genffacgp(is,gqc,ffacgq(:,is))
end do
return
end subroutine

