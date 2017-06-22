
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine zfpack(tpack,n,nr,nri,ld,zfmt,zfir,v)
use modmain
implicit none
! arguments
logical, intent(in) :: tpack
integer, intent(inout) :: n
integer, intent(in) :: nr(nspecies),nri(nspecies)
integer, intent(in) :: ld
complex(8), intent(inout) :: zfmt(lmmaxvr,ld,natmtot)
complex(8), intent(inout) :: zfir(ngtot)
real(8), intent(out) :: v(*)
! local variables
integer is,ias,ir,k
if (tpack) then
! pack the function
  do ias=1,natmtot
    is=idxis(ias)
! inner part of muffin-tin
    k=2*lmmaxinr
    do ir=1,nri(is)
      call dcopy(k,zfmt(:,ir,ias),1,v(n+1),1)
      n=n+k
    end do
! outer part of muffin-tin
    k=2*lmmaxvr*(nr(is)-nri(is))
    call dcopy(k,zfmt(:,nri(is)+1,ias),1,v(n+1),1)
    n=n+k
  end do
  k=2*ngtot
  call dcopy(k,zfir,1,v(n+1),1)
  n=n+k
else
! unpack the function
  do ias=1,natmtot
    is=idxis(ias)
    k=2*lmmaxinr
    do ir=1,nri(is)
      call dcopy(k,v(n+1),1,zfmt(:,ir,ias),1)
      n=n+k
    end do
    k=2*lmmaxvr*(nr(is)-nri(is))
    call dcopy(k,v(n+1),1,zfmt(:,nri(is)+1,ias),1)
    n=n+k
  end do
  k=2*ngtot
  call dcopy(k,v(n+1),1,zfir,1)
  n=n+k
end if
return
end subroutine

