
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rfpack(tpack,n,np,ld,rfmt,rfir,v)
use modmain
implicit none
! arguments
logical, intent(in) :: tpack
integer, intent(inout) :: n
integer, intent(in) :: np(nspecies)
integer, intent(in) :: ld
real(8), intent(inout) :: rfmt(ld,natmtot),rfir(ngtot)
real(8), intent(out) :: v(*)
! local variables
integer is,ias
if (tpack) then
! pack the function
  do ias=1,natmtot
    is=idxis(ias)
    call dcopy(np(is),rfmt(:,ias),1,v(n+1),1)
    n=n+np(is)
  end do
  call dcopy(ngtot,rfir,1,v(n+1),1)
  n=n+ngtot
else
! unpack the function
  do ias=1,natmtot
    is=idxis(ias)
    call dcopy(np(is),v(n+1),1,rfmt(:,ias),1)
    n=n+np(is)
  end do
  call dcopy(ngtot,v(n+1),1,rfir,1)
  n=n+ngtot
end if
return
end subroutine

