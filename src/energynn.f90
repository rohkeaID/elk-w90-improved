
! Copyright (C) 2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine energynn
use modmain
implicit none
! local variables
integer is,ia,ias,ir
real(8) t1
complex(8) zrho0
! automatic arrays
real(8) vn(nrmtmax),vn0(nspecies)
! allocatable arrays
complex(8), allocatable :: zrhoir(:),zvclmt(:,:,:),zvclir(:)
allocate(zrhoir(ngtot))
allocate(zvclmt(lmmaxvr,nrmtmax,natmtot))
allocate(zvclir(ngtot))
! set the density to zero
zrhoir(:)=0.d0
! generate the nuclear monopole potentials
t1=1.d0/y00
do is=1,nspecies
  call potnucl(ptnucl,nrmt(is),rsp(:,is),spzn(is),vn)
  vn0(is)=vn(1)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ir=1,nrmt(is)
      zvclmt(1,ir,ias)=t1*vn(ir)
      zvclmt(2:,ir,ias)=0.d0
    end do
  end do
end do
! solve the complex Poisson's equation
call zpotcoul(nrmt,nrmtinr,nrspmax,rsp,1,gc,jlgr,ylmg,sfacg,zrhoir,nrmtmax, &
 zvclmt,zvclir,zrho0)
! compute the nuclear-nuclear energy
engynn=0.d0
do ias=1,natmtot
  is=idxis(ias)
  t1=dble(zvclmt(1,1,ias))*y00-vn0(is)
  engynn=engynn+spzn(is)*t1
end do
engynn=0.5d0*engynn
deallocate(zrhoir,zvclmt,zvclir)
return
end subroutine

