
! Copyright (C) 2014 K. Krieger, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genafieldt
! !INTERFACE:
subroutine genafieldt
! !USES:
use modmain
use modtddft
! !DESCRIPTION:
!   Generates a time-dependent vector potential, ${\bf A}(t)$, representing a
!   laser pulse and stores it in {\tt AFIELDT.OUT}. The vector potential is
!   constructed from a sum of sinusoidal waves, each modulated with a Gaussian
!   envelope function:
!   $$ {\bf A}(t)={\bf A}_0
!    \frac{e^{-(t-t_0)^2/2\sigma^2}}{\sigma\sqrt{2\pi}}
!    \sin(\omega(t-t_0)+\phi). $$
!   Seven real numbers have to be specified for each pulse, namely the vector
!   amplitude ${\bf A}_0$, peak time $t_0$, full-width at half-maximum
!   $d=2\sqrt{2\ln 2}\sigma$, frequency $\omega$ and phase $\phi$.
!
! !REVISION HISTORY:
!   Created May 2012 (K. Krieger)
!   Modified, January 2014 (S. Sharma)
!   Modified, February 2014 (JKD)
!EOP
!BOC
implicit none
! local variables
integer its,i
real(8) a0(3),t0,d,w,phi,rc
real(8) gs,ppd,s,t,t1,t2,t3
! conversion factor of power density to W/cm^2
real(8), parameter :: cpd=ha_si/(t_si*(100.d0*br_si)**2)
! generate the time-step grid
call gentimes
open(50,file='TD_INFO.OUT',form='FORMATTED')
write(50,*)
write(50,'("(All units are atomic unless otherwise specified)")')
write(50,*)
write(50,'("1 atomic unit of time is ",G18.10," attoseconds")') t_si*1.d18
write(50,*)
write(50,'("Total simulation time : ",G18.10)') tstime
write(50,'(" in attoseconds       : ",G18.10)') tstime*t_si*1.d18
write(50,*)
write(50,'("Time step length : ",G18.10)') dtimes
write(50,'(" in attoseconds  : ",G18.10)') dtimes*t_si*1.d18
write(50,*)
write(50,'("Number of time steps : ",I8)') ntimes
write(50,*)
write(50,'("Number of laser pulses : ",I6)') npulse
write(50,'("Number of ramps : ",I6)') nramp
! allocate and zero time-dependent A-field array
if (allocated(afieldt)) deallocate(afieldt)
allocate(afieldt(3,ntimes))
afieldt(:,:)=0.d0
! loop over laser pulses
do i=1,npulse
! vector amplitude
  a0(1:3)=pulse(1:3,i)
! frequency
  w=pulse(4,i)
! phase
  phi=pulse(5,i)
! chirp rate
  rc=pulse(6,i)
! peak time
  t0=pulse(7,i)
! full-width at half-maximum
  d=pulse(8,i)
! sigma
  s=d/(2.d0*sqrt(2.d0*log(2.d0)))
! write information to TD_INFO.OUT
  write(50,*)
  write(50,'("Pulse : ",I6)') i
  write(50,'(" vector amplitude : ",3G18.10)') a0(:)
  write(50,'(" laser frequency : ",G18.10)') w
  write(50,'("  in eV          : ",G18.10)') w*ha_ev
  write(50,'(" laser wavelength (Angstroms) : ",G18.10)') 1.d10/(w*ha_im)
  write(50,'(" phase (degrees) : ",G18.10)') phi
  write(50,'(" chirp rate : ",G18.10)') rc
  write(50,'(" peak time : ",G18.10)') t0
  write(50,'(" full-width at half-maximum : ",G18.10)') d
  write(50,'(" Gaussian sigma = FWHM/2*sqrt(2*ln 2) : ",G18.10)') s
  t1=a0(1)**2+a0(2)**2+a0(3)**2
  ppd=t1*(w**2)/(8.d0*pi*solsc)
  write(50,'(" peak laser power density : ",G18.10)') ppd
  write(50,'("  in W/cm^2               : ",G18.10)') ppd*cpd
! loop over time steps
  do its=1,ntimes
    t=times(its)
    t1=t-t0
    t2=-0.5d0*(t1/s)**2
    t3=w*t1+phi*pi/180.d0+0.5d0*rc*t**2
    gs=exp(t2)*sin(t3)
    if (abs(gs).lt.1.d-20) gs=0.d0
    afieldt(:,its)=afieldt(:,its)+a0(:)*gs
  end do
end do
! loop over ramps
do i=1,nramp
! vector amplitude
  a0(1:3)=ramp(1:3,i)
! ramp start time
  t0=ramp(4,i)
! linear coefficient
  t1=ramp(5,i)
! quadratic coefficient
  t2=ramp(6,i)
! write information to TD_INFO.OUT
  write(50,*)
  write(50,'("Ramp : ",I6)') i
  write(50,'(" vector amplitude : ",3G18.10)') a0(:)
  write(50,'(" ramp start time : ",G18.10)') t0
! loop over time steps
  do its=1,ntimes
    t=times(its)
    if (t.gt.t0) then
      t3=t1*(t-t0)+t2*(t-t0)**2
      afieldt(:,its)=afieldt(:,its)+a0(:)*t3
    end if
  end do
end do
close(50)
! write the vector potential to AFIELDT.OUT
open(50,file='AFIELDT.OUT',form='FORMATTED')
write(50,'(I8," : number of time steps")') ntimes
do its=1,ntimes
  write(50,'(I8,4G18.10)') its,times(its),afieldt(:,its)
end do
close(50)
write(*,*)
write(*,'("Info(genafieldt):")')
write(*,'(" Time-dependent A-field written to AFIELDT.OUT")')
write(*,'(" Laser pulse and ramp parameters written to TD_INFO.OUT")')
write(*,*)
write(*,'(" 1 atomic unit of time is ",G18.10," attoseconds")') t_si*1.d18
write(*,'(" Total simulation time : ",G18.10)') tstime
write(*,'("  in attoseconds       : ",G18.10)') tstime*t_si*1.d18
deallocate(times,afieldt)
return
end subroutine
!EOC

