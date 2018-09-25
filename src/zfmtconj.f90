
! Copyright (C) 2016 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine zfmtconj(nr,nri,np,zfmt)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri,np
complex(8), intent(inout) :: zfmt(np)
! local variables
integer ir,i
! automatic arrays
complex(8) zfmt1(np)
call zcopy(np,zfmt,1,zfmt1,1)
i=1
do ir=1,nri
  call zflmconj(lmaxi,zfmt1(i),zfmt(i))
  i=i+lmmaxi
end do
do ir=nri+1,nr
  call zflmconj(lmaxo,zfmt1(i),zfmt(i))
  i=i+lmmaxo
end do
return
end subroutine

