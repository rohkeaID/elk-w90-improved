
! Copyright (C) 2009 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genspecies(fnum)
use modmain
use modmpi
implicit none
! arguments
integer, intent(in) :: fnum
! local variables
integer nz,nmax,nst,ist
integer ne,nrm,nr,ir,it,i
integer n(maxstsp),l(maxstsp),k(maxstsp)
real(8) mass,zn,t1,t2,t3
real(8) rm,rmin,rmax
real(8) occ(maxstsp),eval(maxstsp)
character(256) symb,name
! allocatable arrays
real(8), allocatable :: r(:),rho(:),vr(:),rwf(:,:,:)
read(fnum,*,err=20) nz
if (nz.le.0) then
  write(*,*)
  write(*,'("Error(genspecies): atomic number negative : ",I8)') nz
  write(*,*)
  stop
end if
read(fnum,*,err=20) symb,name
read(fnum,*,err=20) mass
! convert from 'atomic mass units' to atomic units
mass=mass*amu
read(fnum,*,err=20) rm
read(fnum,*,err=20) nst
if ((nst.le.0).or.(nst.gt.maxstsp)) then
  write(*,*)
  write(*,'("Error(genspecies): nst out of range : ",I8)') nst
  write(*,'(" for species ",A)') trim(name)
  write(*,*)
  stop
end if
ne=0
nmax=1
do ist=1,nst
  read(fnum,*,err=20) n(ist),l(ist),k(ist),i
  ne=ne+i
  occ(ist)=i
  nmax=max(nmax,n(ist))
end do
if (mp_mpi) then
  write(*,'("Info(genspecies): running Z = ",I4,", (",A,")")') nz,trim(name)
  if (ne.ne.nz) then
    write(*,*)
    write(*,'("Warning(genspecies): atom not neutral, electron number : ",&
     &I4)') ne
  end if
end if
! nuclear charge in units of e
zn=-dble(nz)
! minimum radial mesh point proportional to 1/sqrt(Z)
rmin=2.d-6/sqrt(dble(nz))
! default effective infinity
rmax=100.d0
! set the number of radial mesh points proportional to number of nodes
nrm=100*(nmax+1)
do it=1,2
! number of points to effective infinity
  t1=log(rm/rmin)
  t2=log(rmax/rmin)
  t3=dble(nrm)*t2/t1
  nr=int(t3)
  allocate(r(nr),rho(nr),vr(nr),rwf(nr,2,nst))
! generate logarithmic radial mesh
  t2=t1/dble(nrm-1)
  do ir=1,nr
    r(ir)=rmin*exp(dble(ir-1)*t2)
  end do
! solve the Kohn-Sham-Dirac equation for the atom
  call atom(sol,.true.,zn,nst,n,l,k,occ,3,0,nr,r,eval,rho,vr,rwf)
! recompute the effective infinity
  do ir=nr,1,-1
    if (rho(ir).gt.1.d-20) then
      rmax=1.75d0*r(ir)
      exit
    end if
  end do
  deallocate(r,rho,vr,rwf)
end do
! write the species file
call writespecies(symb,name,zn,mass,rmin,rm,rmax,nrm,nst,n,l,k,occ,eval)
return
20 continue
write(*,*)
write(*,'("Error(genspecies): error reading species data")')
write(*,*)
stop
end subroutine

