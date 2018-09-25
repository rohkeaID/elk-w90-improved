
! Copyright (C) 2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine energynn
use modmain
implicit none
! local variables
integer is,ia,ias
integer nr,nri,ir,np,i
real(8) t1
! automatic arrays
real(8) vn(nrmtmax),vn0(nspecies)
! allocatable arrays
complex(8), allocatable :: zvclmt(:,:),zvclir(:),zrhoir(:)
allocate(zvclmt(npmtmax,natmtot),zvclir(ngtot))
! generate the nuclear monopole potentials
t1=1.d0/y00
do is=1,nspecies
  nr=nrmt(is)
  nri=nrmti(is)
  np=npmt(is)
  call potnucl(ptnucl,nr,rsp(:,is),spzn(is),vn)
  vn0(is)=vn(1)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    zvclmt(1:np,ias)=0.d0
    i=1
    do ir=1,nri
      zvclmt(i,ias)=t1*vn(ir)
      i=i+lmmaxi
    end do
    do ir=nri+1,nr
      zvclmt(i,ias)=t1*vn(ir)
      i=i+lmmaxo
    end do
  end do
end do
! set the interstitial density to zero
allocate(zrhoir(ngtot))
zrhoir(:)=0.d0
! solve the complex Poisson's equation
call zpotcoul(nrmt,nrmti,npmt,npmti,nrspmax,rsp,ngridg,igfft,ngvec,gc,gclg, &
 ngvec,jlgrmt,ylmg,sfacg,zrhoir,npmtmax,zvclmt,zvclir)
! compute the nuclear-nuclear energy
engynn=0.d0
do ias=1,natmtot
  is=idxis(ias)
  t1=dble(zvclmt(1,ias))*y00-vn0(is)
  engynn=engynn+spzn(is)*t1
end do
engynn=0.5d0*engynn
deallocate(zvclmt,zvclir,zrhoir)
return
end subroutine

