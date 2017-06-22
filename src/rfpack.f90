
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rfpack(tpack,n,nr,nri,ld,rfmt,rfir,v)
use modmain
implicit none
! arguments
logical, intent(in) :: tpack
integer, intent(inout) :: n
integer, intent(in) :: nr(nspecies),nri(nspecies)
integer, intent(in) :: ld
real(8), intent(inout) :: rfmt(lmmaxvr,ld,natmtot),rfir(ngtot)
real(8), intent(out) :: v(*)
! local variables
integer is,ias,ir,k
if (tpack) then
! pack the function
  do ias=1,natmtot
    is=idxis(ias)
! inner part of muffin-tin
    do ir=1,nri(is)
      call dcopy(lmmaxinr,rfmt(:,ir,ias),1,v(n+1),1)
      n=n+lmmaxinr
    end do
! outer part of muffin-tin
    k=lmmaxvr*(nr(is)-nri(is))
    call dcopy(k,rfmt(:,nri(is)+1,ias),1,v(n+1),1)
    n=n+k
  end do
  call dcopy(ngtot,rfir,1,v(n+1),1)
  n=n+ngtot
else
! unpack the function
  do ias=1,natmtot
    is=idxis(ias)
    do ir=1,nri(is)
      call dcopy(lmmaxinr,v(n+1),1,rfmt(:,ir,ias),1)
      n=n+lmmaxinr
    end do
    k=lmmaxvr*(nr(is)-nri(is))
    call dcopy(k,v(n+1),1,rfmt(:,nri(is)+1,ias),1)
    n=n+k
  end do
  call dcopy(ngtot,v(n+1),1,rfir,1)
  n=n+ngtot
end if
return
end subroutine

