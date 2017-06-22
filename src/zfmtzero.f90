
! Copyright (C) 2015 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine zfmtzero(nr,nri,zfmt)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri
complex(8), intent(out) :: zfmt(lmmaxvr,nr)
zfmt(1:lmmaxinr,1:nri)=0.d0
zfmt(:,nri+1:)=0.d0
return
end subroutine

