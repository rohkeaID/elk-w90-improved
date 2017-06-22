
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: findband
! !INTERFACE:
subroutine findband(sol,l,nr,r,vr,eps,demax,e,fnd)
! !INPUT/OUTPUT PARAMETERS:
!   sol   : speed of light in atomic units (in,real)
!   l     : angular momentum quantum number (in,integer)
!   nr    : number of radial mesh points (in,integer)
!   r     : radial mesh (in,real(nr))
!   vr    : potential on radial mesh (in,real(nr))
!   eps   : energy search tolerance (in,real)
!   demax : maximum allowed change from the input energy; enforced only if e < 0
!           (in,real)
!   e     : input energy and returned band energy (inout,real)
!   fnd   : set to .true. if the band energy is found (out,logical)
! !DESCRIPTION:
!   Finds the band energies for a given radial potential and angular momentum.
!   This is done by first searching upwards in energy, starting from the input
!   energy plus the offset energy, until the radial wavefunction at the
!   muffin-tin radius is zero. This is the energy at the top of the band,
!   denoted $E_{\rm t}$. A downward search is now performed from $E_{\rm t}$
!   until the slope of the radial wavefunction at the muffin-tin radius is zero.
!   This energy, $E_{\rm b}$, is at the bottom of the band. The band energy is
!   taken as $(E_{\rm t}+E_{\rm b})/2$. If either $E_{\rm t}$ or $E_{\rm b}$
!   cannot be found then the band energy is set to the input value.
!
! !REVISION HISTORY:
!   Created September 2004 (JKD)
!   Added two-pass loop, October 2013 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: sol
integer, intent(in) :: l,nr
real(8), intent(in) :: r(nr),vr(nr)
real(8), intent(in) :: eps,demax
real(8), intent(inout) :: e
logical, intent(out) :: fnd
! local variables
logical ft,fb
! maximum number of steps
integer, parameter :: maxstp=250
integer ip,ie,nn
! initial step size
real(8), parameter :: de0=0.001d0
real(8) de,et,eb,t,tp
! automatic arrays
real(8) p0(nr),p1(nr),q0(nr),q1(nr)
ft=.false.
fb=.false.
fnd=.false.
et=e
eb=e
! two-pass loop
do ip=1,2
! find the top of the band
  tp=0.d0
  de=de0
  do ie=1,maxstp
    et=et+de
    if (e.lt.0.d0) then
      if (et.gt.e+demax) exit
    end if
    call rschrodint(sol,l,et,nr,r,vr,nn,p0,p1,q0,q1)
    t=p0(nr)
    if (ie.gt.1) then
      if (t*tp.le.0.d0) then
        if (abs(de).lt.eps) then
          if (fb) goto 10
          ft=.true.
          eb=et+5.d0*abs(de0)
          exit
        end if
        de=-0.5d0*de
      else
        de=1.5d0*de
      end if
    end if
    tp=t
  end do
  if (fb) return
! find the bottom of the band
  tp=0.d0
  de=-de0
  do ie=1,maxstp
    eb=eb+de
    if (eb.lt.e-demax) return
    call rschrodint(sol,l,eb,nr,r,vr,nn,p0,p1,q0,q1)
    t=p1(nr)
    if (ie.gt.1) then
      if (t*tp.le.0.d0) then
        if (abs(de).lt.eps) then
          if (ft) goto 10
          fb=.true.
          et=eb-5.d0*abs(de0)
          exit
        end if
        de=-0.5d0*de
      else
        de=1.5d0*de
      end if
    end if
    tp=t
  end do
end do
return
10 continue
! set the band energy halfway between top and bottom
e=(et+eb)/2.d0
fnd=.true.
return
end subroutine
!EOC

