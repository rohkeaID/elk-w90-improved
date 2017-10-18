
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gradwfcr2(gwf2mt)
use modmain
implicit none
! arguments
real(8), intent(inout) :: gwf2mt(npmtmax,natmtot)
! local variables
integer ist,is,ias
integer nr,nri,iro,ir
integer np,l,m,lm,i
! allocatable arrays
complex(8), allocatable :: wfmt(:),gwfmt(:,:),zfmt(:)
allocate(wfmt(npmtmax),gwfmt(npmtmax,3),zfmt(npmtmax))
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  iro=nri+1
  np=npmt(is)
  do ist=1,nstsp(is)
    if (spcore(ist,is).and.(ksp(ist,is).eq.lsp(ist,is)+1)) then
      l=lsp(ist,is)
      do m=-l,l
        lm=idxlm(l,m)
        wfmt(1:np)=0.d0
        i=lm
        do ir=1,nri
          wfmt(i)=rwfcr(ir,1,ist,ias)/rsp(ir,is)
          i=i+lmmaxi
        end do
        do ir=iro,nr
          wfmt(i)=rwfcr(ir,1,ist,ias)/rsp(ir,is)
          i=i+lmmaxo
        end do
        call gradzfmt(nr,nri,rsp(:,is),wfmt,npmtmax,gwfmt)
        do i=1,3
          call zbsht(nr,nri,gwfmt(:,i),zfmt)
! factor of 2 from spin
          gwf2mt(1:np,ias)=gwf2mt(1:np,ias) &
           +2.d0*(dble(zfmt(1:np))**2+aimag(zfmt(1:np))**2)
        end do
      end do
    end if
  end do
! end loops over atoms
end do
deallocate(wfmt,gwfmt,zfmt)
return
end subroutine

