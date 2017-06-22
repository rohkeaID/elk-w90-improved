
! Copyright (C) 2014 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genfspecies(zn,symb)
use modmain
use modmpi
implicit none
! arguments
real(8), intent(in) :: zn
character(*), intent(in) :: symb
! local variables
integer, parameter :: nit=4
integer nst,ist,jst
integer nmax,in,il,ik
integer nrm,nr,ir,it
integer n(maxstsp),l(maxstsp),k(maxstsp)
integer idx(maxstsp),iv(maxstsp)
real(8) rm,rmin,rmax
real(8) mass,t1,t2,t3
real(8) occ(maxstsp),eval(maxstsp),rv(maxstsp)
character(256) name
! allocatable arrays
real(8), allocatable :: r(:),rho(:),vr(:),rwf(:,:,:)
! external functions
real(8) massnucl
external massnucl
name='Fractional species'
! set up the initial occupancies
occ(:)=0.d0
t1=abs(zn)
nmax=1
ist=0
do in=1,maxstsp
  do il=0,in-1
    do ik=max(il,1),il+1
      t2=dble(2*ik)
      t2=min(t2,t1)
      ist=ist+1
      n(ist)=in
      l(ist)=il
      k(ist)=ik
      occ(ist)=t2
      if (t2.gt.epsocc) nmax=in
      t1=t1-t2
      if (ist.eq.maxstsp) then
        if (t1.gt.epsocc) then
          write(*,*)
          write(*,'("Error(genfspecies): too many states for fractional &
           &species ",A)') trim(symb)
          write(*,*)
          stop
        else
          goto 10
        end if
      end if
    end do
  end do
end do
10 continue
! minimum radius
rmin=2.d-6/sqrt(abs(zn))
! initial maximum radius
rmax=100.d0
! initial muffin-tin radius
rm=2.d0
! number of points to muffin-tin radius
nrm=100*(nmax+1)
! iterate the solution but not to self-consistency
do it=1,nit
! number of points to effective infinity
  t1=log(rm/rmin)
  t2=log(rmax/rmin)
  t3=dble(nrm)*t2/t1
  nr=int(t3)
  allocate(r(nr),rho(nr),vr(nr),rwf(nr,2,maxstsp))
! generate logarithmic radial mesh
  t2=t1/dble(nrm-1)
  do ir=1,nr
    r(ir)=rmin*exp(dble(ir-1)*t2)
  end do
! solve the Kohn-Sham-Dirac equation for the atom
  call atom(sol,.true.,zn,maxstsp,n,l,k,occ,3,0,nr,r,eval,rho,vr,rwf)
! check for spurious eigenvalues
  do ist=2,maxstsp
    if (eval(ist).lt.eval(1)) eval(ist)=1.d6
  end do
! recompute the effective infinity
  do ir=nr,1,-1
    if (rho(ir).gt.1.d-20) then
      rmax=1.75d0*r(ir)
      exit
    end if
  end do
! estimate the muffin-tin radius
  do ir=nr,1,-1
    if (rho(ir).gt.2.d-2) then
      rm=r(ir)
      exit
    end if
  end do
  if (rm.lt.1.d0) rm=1.d0
  if (rm.gt.3.2d0) rm=3.2d0
! sort the eigenvalues
  call sortidx(maxstsp,eval,idx)
! recompute the occupancies
  occ(:)=0.d0
  t1=abs(zn)
  do ist=1,maxstsp
    jst=idx(ist)
    ik=k(jst)
    t2=dble(2*ik)
    t2=min(t2,t1)
    occ(jst)=t2
    t1=t1-t2
  end do
  deallocate(r,rho,vr,rwf)
end do
! rearrange the arrays
iv(:)=n(:)
n(:)=iv(idx(:))
iv(:)=l(:)
l(:)=iv(idx(:))
iv(:)=k(:)
k(:)=iv(idx(:))
rv(:)=occ(:)
occ(:)=rv(idx(:))
rv(:)=eval(:)
eval(:)=rv(idx(:))
! find the number of occupied states
nst=0
do ist=1,maxstsp
  if (occ(ist).lt.epsocc) then
    nst=ist
    exit
  end if
end do
! estimate the nuclear mass
mass=massnucl(zn)
! convert from 'atomic mass units' to atomic units
mass=mass*amu
! write the species file
call writespecies(symb,name,zn,mass,rmin,rm,rmax,nrm,nst,n,l,k,occ,eval)
if (mp_mpi) then
  write(*,*)
  write(*,'("Info(genfspecies): wrote fractional species file ",A,".in")') &
   trim(symb)
end if
return
end subroutine

