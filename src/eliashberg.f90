
! Copyright (C) 2011 A. Sanna and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: eliashberg
! !INTERFACE:
subroutine eliashberg
! !USES:
use modmain
use modphonon
! !DESCRIPTION:
!   Calculates the superconducting gap within Eliashberg theory. This
!   implementation is isotropic and assumes a flat density of states. The
!   Eliashberg function $\alpha^2F$ is required as input for this calculation.
!
! !REVISION HISTORY:
!   Created December 2010 (Antonio Sanna)
!   Modified, June 2011 (JKD)
!EOP
!BOC
implicit none
! local variables
! maximum allowed number of Matsubara frequencies
integer, parameter :: maxwf=40000
! maximum number of iterations
integer, parameter :: maxit=1000
integer nwf,nwfcl,nin,nout
integer itemp,it,i,m,n
! convergence tolerance
real(8), parameter :: eps=1.d-12
! mixing paramter
real(8), parameter :: beta=0.5d0
real(8) lambda,wlog,wrms,tc
real(8) wfmax,tmin,tmax,dtemp,temp
real(8) dw,dmu,sum,a,b,t0,t1
! allocatable arrays
real(8), allocatable :: w(:),a2f(:),wf(:),l(:)
real(8), allocatable :: d0(:),d(:),z0(:),z(:),r(:)
complex(8), allocatable :: zin(:),uin(:),zout(:),uout(:)
! initialise universal variables
call init0
call init1
allocate(w(nwplot),a2f(nwplot))
! read in the Eliashberg function
call readalpha2f(w,a2f)
dw=(w(nwplot)-w(1))/dble(nwplot)
! compute the McMillan parameters
call mcmillan(w,a2f,lambda,wlog,wrms,tc)
! Matsubara frequency cut-off
wfmax=20.d0*wrms
! minumum temperature
tmin=tc/6.d0
if (tmin.lt.1.d-2) tmin=0.1d0
! maximum temperature
tmax=3.d0*tc
if (tmax.lt.1.d0) tmax=1.d0
! temperature step size
dtemp=(tmax-tmin)/dble(ntemp)
! maximum number of fermionic Matsubara frequencies
nwf=nint(wfmax/(twopi*kboltz*dtemp))
if (nwf.lt.1) nwf=1
if (nwf.gt.maxwf) nwf=maxwf
allocate(wf(-nwf:nwf))
allocate(l(-2*nwf:2*nwf))
allocate(d0(0:nwf),d(0:nwf))
allocate(z0(0:nwf),z(0:nwf))
allocate(r(0:nwf))
allocate(zin(0:nwf),uin(0:nwf))
! generate output points for continuation on the real axis
nout=4*nwplot
allocate(zout(nout),uout(nout))
do i=1,nout
  zout(i)=cmplx(2.d0*dble(i-1)*dw,0.d0,8)
end do
! open files for writing
open(62,file='ELIASHBERG.OUT',action='WRITE',form='FORMATTED')
open(63,file='ELIASHBERG_IA.OUT',action='WRITE',form='FORMATTED')
open(64,file='ELIASHBERG_GAP_T.OUT',action='WRITE',form='FORMATTED')
open(65,file='ELIASHBERG_GAP_RA.OUT',action='WRITE',form='FORMATTED')
open(66,file='ELIASHBERG_Z_RA.OUT',action='WRITE',form='FORMATTED')
write(62,'("+------------------------------+")')
write(62,'("|     Eliashberg equations     |")')
write(62,'("+------------------------------+")')
write(62,*)
write(62,'("Temperature range : ",2G18.10)') tmin,tmax
write(62,'("Number of temperature steps : ",I6)') ntemp
write(62,'("Number of output frequencies : ",I8)') nout
write(62,'("Fermionic Matsubara frequency cut-off")')
write(62,'(" phonons : ",G18.10)') wfmax
write(62,'(" Coulomb : ",G18.10)') wrms
call flushifc(62)
d0(:)=1.d-4
z0(:)=1.d0
! main loop over temperature
do itemp=1,ntemp
  write(*,'("Info(eliashberg): temperature step ",I6," of ",I6)') itemp,ntemp
  temp=dble(itemp)*dtemp+tmin
  write(62,*)
  write(62,'("Temperature (kelvin) : ",G18.10)') temp
  t0=pi*kboltz*temp
! number of Matsubara frequencies
  nwf=nint(wfmax/(2.d0*t0))
  if (nwf.gt.maxwf) nwf=maxwf
  nwfcl=nint(wrms/(2.d0*t0))
  if (nwfcl.lt.1) nwfcl=1
  if (nwfcl.gt.nwf) nwfcl=nwf
  write(62,'("Number of Matsubara frequencies")')
  write(62,'(" phonons : ",I8)') nwf
  write(62,'(" Coulomb : ",I8)') nwfcl
! make Pade approximant input points same as Matsubara frequencies
  nin=nwf
! generate fermionic Matsubara frequencies
  do m=-nwf,nwf
    wf(m)=t0*dble(2*m+1)
  end do
! compute lambda
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(t1,sum,i)
!$OMP DO
  do m=-2*nwf,2*nwf
    t1=(t0*dble(2*m))**2
    sum=0.d0
    do i=1,nwplot
      sum=sum+w(i)*a2f(i)/(w(i)**2+t1)
    end do
    l(m)=2.d0*sum*dw
  end do
!$OMP END DO
!$OMP END PARALLEL
! begin iteration loop
  do it=1,maxit
    do m=0,nwf
      r(m)=sqrt((wf(m)**2+d0(m)**2)*z0(m)**2)
    end do
    do n=0,nwf
      sum=0.d0
      do m=0,nwf-1
        sum=sum+(l(n-m)-l(n+m+1))*z0(m)*wf(m)/r(m)
      end do
      z(n)=t0*sum/wf(n)
    end do
    z(0:nwf)=z(0:nwf)+1.d0
    z0(0:nwf)=z(0:nwf)
! Coulomb part of summation
    dmu=0.d0
    do n=0,nwfcl
      dmu=dmu+mustar*d0(n)*z(n)/r(n)
    end do
    dmu=dmu*2.d0
! Gap
    do n=0,nwf
      sum=0.d0
      do m=0,nwf-1
        sum=sum+(l(n-m)+l(n+m+1))*d0(m)*z(m)/r(m)
      end do
      d(n)=t0*(sum-dmu)/z(n)
    end do
! mix old and new gap functions
    d(0:nwf)=beta*d(0:nwf)+(1.d0-beta)*d0(0:nwf)
    sum=0.d0
    do m=0,nwf
      sum=sum+abs(d0(m)-d(m))
    end do
    sum=sum/dble(2*nwf)
    d0(0:nwf)=d(0:nwf)
    if (sum.le.eps) then
      write(62,'("Eliashberg equations converged in ",I6," iterations")') it
      goto 10
    end if
! end iteration loop
  end do
  write(*,*)
  write(*,'("Warning(eliashberg): failed to converge: possibly close to T_c")')
  write(62,'("Failed to converge: possibly close to T_c")')
10 continue
  call flushifc(62)
  do n=-nwf,nwf
    if (n.ge.0) then
      m=n
    else
      m=-n-1
    end if
    write(63,'(3G18.10)') wf(n),d(m),z(m)
  end do
  write(63,*)
  call flushifc(63)
  write(64,'(3G18.10)') temp,d(0),z(0)
  call flushifc(64)
! analytic continuation to real axis
  do m=0,nin
    zin(m)=cmplx(0.d0,wf(m),8)
    uin(m)=cmplx(d(m),0.d0,8)
  end do
  call pade(nin,zin,uin,nout,zout,uout)
  do i=1,nout
    a=dble(uout(i))
    b=aimag(uout(i))
    write(65,'(3G18.10)') dble(zout(i)),a,b
  end do
  write(65,*)
  call flushifc(65)
  do m=0,nin
    uin(m)=cmplx(z(m),0.d0,8)
  end do
  call pade(nin,zin,uin,nout,zout,uout)
  do i=1,nout
    a=dble(uout(i))
    b=aimag(uout(i))
    write(66,'(3G18.10)') dble(zout(i)),a,b
  end do
  write(66,*)
  call flushifc(66)
! end loop over temperatures
end do
close(62); close(63); close(64); close(65); close(66)
write(*,*)
write(*,'("Info(eliashberg):")')
write(*,'(" calculation information written to ELIASHBERG.OUT")')
write(*,'(" gap and Z functions on the imaginary axis written to &
 &ELIASHBERG_IA.OUT")')
write(*,'(" gap vs. temperature written to ELIASHBERG_GAP_T.OUT")')
write(*,'(" gap function on the real axis written to ELIASHBERG_GAP_RA.OUT")')
write(*,'(" Z function on the real axis written to ELIASHBERG_Z_RA.OUT")')
deallocate(w,a2f,wf,l)
deallocate(d0,d,z0,z,r)
deallocate(zin,uin,zout,uout)
return
end subroutine
!EOC

