
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dpotcoul
use modmain
use modphonon
implicit none
! local variables
integer nr,nri,ir
integer np,i,igq0
complex(8) zgq0
! automatic arrays
real(8) vn(nrmtmax)
! allocatable arrays
complex(8), allocatable :: zfmt(:),gzfmt(:,:)
nr=nrmt(isph)
nri=nrmti(isph)
np=npmt(isph)
! solve the complex Poisson's equation in the muffin-tins
call genzvclmt(nrmt,nrmti,nrspmax,rsp,npmtmax,drhomt,dvclmt)
! compute the monopole potential for effective nuclear charge
call potnucl(ptnucl,nr,rsp(:,isph),spzn(isph),vn)
! calculate the gradient of the monopole potential
allocate(zfmt(npmtmax),gzfmt(npmtmax,3))
zfmt(1:np)=0.d0
i=1
do ir=1,nri
  zfmt(i)=vn(ir)/y00
  i=i+lmmaxi
end do
do ir=nri+1,nr
  zfmt(i)=vn(ir)/y00
  i=i+lmmaxo
end do
call gradzfmt(nr,nri,rsp(:,isph),zfmt,npmtmax,gzfmt)
! subtract gradient component corresponding to the phonon polarisation
dvclmt(1:np,iasph)=dvclmt(1:np,iasph)-gzfmt(1:np,ipph)
deallocate(zfmt,gzfmt)
! solve Poisson's equation in the entire unit cell
if (tphq0) then
  igq0=1
else
  igq0=0
end if
call zpotcoul(nrmt,nrmti,npmt,npmti,nrspmax,rsp,ngvec,igq0,gqc,ngvec,jlgqrmt, &
 ylmgq,sfacgq,drhoir,npmtmax,dvclmt,dvclir,zgq0)
return
end subroutine

