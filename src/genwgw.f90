
! Copyright (C) 2016 A. Davydov, A. Sanna, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine genwgw
use modmain
use modgw
implicit none
! local variables
integer ik,ist,n
integer iw,jw,it
real(8) de,t0,t1
t0=kboltz*tempk
if (wmaxgw.le.0.d0) then
! read Fermi energy from file
  call readfermi
! find the maximum eigenvalue range over all k-points
  de=0.d0
  do ik=1,nkpt
    call getevalsv(filext,ik,vkl(:,ik),evalsv(:,ik))
    do ist=1,nstsv
      de=max(de,abs(evalsv(ist,ik)-efermi))
    end do
  end do
  wmaxgw=abs(wmaxgw)*de
end if
! number of Matsubara frequencies
t1=pi*t0
nwgw=2*nint(wmaxgw/t1)
nwgw=max(nwgw,2)
call nfftifc(nwgw)
! determine integer ranges for grid
intwgw(1)=nwgw/2-nwgw+1
intwgw(2)=nwgw/2
if (allocated(iwfft)) deallocate(iwfft)
allocate(iwfft(intwgw(1):intwgw(2)))
if (allocated(wgw)) deallocate(wgw)
allocate(wgw(intwgw(1):intwgw(2)))
do iw=intwgw(1),intwgw(2)
  if (iw.ge.0) then
    jw=iw
  else
    jw=nwgw+iw
  end if
  iwfft(iw)=jw+1
  wgw(iw)=dble(iw)*t1
end do
n=minval(abs(intwgw(:)))
if (n.eq.0) then
  write(*,*)
  write(*,'("Error(genwgw): not enough Matsubara frequencies")')
  write(*,'("Increase wmaxgw")')
  write(*,*)
  stop
end if
if (mod(n,2).eq.0) then
  nwbs=n
  nwfm=n-1
else
  nwfm=n
  nwbs=n-1
end if
! store the complex fermionic frequencies
if (allocated(wfm)) deallocate(wfm)
allocate(wfm(0:nwfm))
do iw=-nwfm,nwfm,2
  jw=(iw+nwfm)/2
  wfm(jw)=cmplx(0.d0,wgw(iw),8)
end do
! store the complex tau-points
if (allocated(taugw)) deallocate(taugw)
allocate(taugw(nwgw))
t1=2.d0/(t0*dble(nwgw))
do it=1,nwgw
  taugw(it)=t1*dble(it-1)
end do
! store the complex response function frequencies
nwrf=nwbs+1
if (allocated(wrf)) deallocate(wrf)
allocate(wrf(nwrf))
do iw=-nwbs,nwbs,2
  jw=(iw+nwbs)/2+1
  wrf(jw)=cmplx(0.d0,wgw(iw),8)
end do
return
end subroutine

