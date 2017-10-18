
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine zfpack(tpack,n,np,ld,zfmt,zfir,v)
use modmain
implicit none
! arguments
logical, intent(in) :: tpack
integer, intent(inout) :: n
integer, intent(in) :: np(nspecies)
integer, intent(in) :: ld
complex(8), intent(inout) :: zfmt(ld,natmtot),zfir(ngtot)
real(8), intent(out) :: v(*)
! local variables
integer is,ias,k
if (tpack) then
! pack the function
  do ias=1,natmtot
    is=idxis(ias)
    k=2*np(is)
    call dcopy(k,zfmt(:,ias),1,v(n+1),1)
    n=n+k
  end do
  k=2*ngtot
  call dcopy(k,zfir,1,v(n+1),1)
  n=n+k
else
! unpack the function
  do ias=1,natmtot
    is=idxis(ias)
    k=2*np(is)
    call dcopy(k,v(n+1),1,zfmt(:,ias),1)
    n=n+k
  end do
  k=2*ngtot
  call dcopy(k,v(n+1),1,zfir,1)
  n=n+k
end if
return
end subroutine

