
! Copyright (C) 2016 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine zfmtftoc(nr,nri,zfmt1,zfmt2)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri
complex(8), intent(in) :: zfmt1(*)
complex(8), intent(out) :: zfmt2(*)
! local variables
integer ir,i1,i2,j1
i1=1
i2=1
j1=lmmaxi*lradstp
do ir=1,nri,lradstp
  call zcopy(lmmaxi,zfmt1(i1),1,zfmt2(i2),1)
  i1=i1+j1
  i2=i2+lmmaxi
end do
i1=i1+(lradstp-1)*(lmmaxo-lmmaxi)
j1=lmmaxo*lradstp
do ir=nri+lradstp,nr,lradstp
  call zcopy(lmmaxo,zfmt1(i1),1,zfmt2(i2),1)
  i1=i1+j1
  i2=i2+lmmaxo
end do
return
end subroutine

