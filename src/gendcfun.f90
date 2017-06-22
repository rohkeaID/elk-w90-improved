
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gendcfun
use modmain
use modphonon
implicit none
! local variables
integer ig
real(8) t1,t2
complex(8) z1,z2
do ig=1,ngtot
  t1=-dot_product(vgqc(:,ig),atposc(:,iaph,isph))
  t2=ffacgq(ig,isph)*vgqc(ipph,ig)
  z1=t2*cmplx(cos(t1),sin(t1),8)
  z2=cmplx(-aimag(z1),dble(z1),8)
  dcfunig(ig)=z2
  dcfunir(igfft(ig))=z2
end do
call zfftifc(3,ngridg,1,dcfunir)
return
end subroutine

