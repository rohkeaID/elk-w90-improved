
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gradwfcr2(gwf2mt)
use modmain
implicit none
! arguments
real(8), intent(inout) :: gwf2mt(lmmaxvr,nrmtmax,natmtot)
! local variables
integer ist,is,ias
integer nr,nri,ir
integer l,m,lm,i
! allocatable arrays
complex(8), allocatable :: wfmt(:,:),gwfmt(:,:,:),zfmt(:,:)
allocate(wfmt(lmmaxvr,nrmtmax))
allocate(gwfmt(lmmaxvr,nrmtmax,3))
allocate(zfmt(lmmaxvr,nrmtmax))
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmtinr(is)
  do ist=1,nstsp(is)
    if (spcore(ist,is).and.(ksp(ist,is).eq.lsp(ist,is)+1)) then
      l=lsp(ist,is)
      do m=-l,l
        lm=idxlm(l,m)
        wfmt(:,1:nr)=0.d0
        do ir=1,nr
          wfmt(lm,ir)=rwfcr(ir,1,ist,ias)/rsp(ir,is)
        end do
        call gradzfmt(nr,nri,rsp(:,is),wfmt,nrmtmax,gwfmt)
        do i=1,3
          call zbsht(nr,nri,gwfmt(:,:,i),zfmt)
! inner part of muffin-tin
          do ir=1,nri
! factor of 2 from spin
            gwf2mt(1:lmmaxinr,ir,ias)=gwf2mt(1:lmmaxinr,ir,ias) &
             +2.d0*(dble(zfmt(1:lmmaxinr,ir))**2+aimag(zfmt(1:lmmaxinr,ir))**2)
          end do
! outer part of muffin tin
          do ir=nri+1,nr
            gwf2mt(:,ir,ias)=gwf2mt(:,ir,ias) &
             +2.d0*(dble(zfmt(:,ir))**2+aimag(zfmt(:,ir))**2)
          end do
        end do
      end do
    end if
  end do
! end loops over atoms
end do
deallocate(wfmt,gwfmt,zfmt)
return
end subroutine

