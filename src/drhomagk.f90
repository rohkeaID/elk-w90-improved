
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine drhomagk(ngp,ngpq,igpig,igpqig,occsvp,doccsvp,apwalm,apwalmq, &
 dapwalm,evecfv,devecfv,evecsv,devecsv)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: ngp(nspnfv),ngpq(nspnfv)
integer, intent(in) :: igpig(ngkmax,nspnfv),igpqig(ngkmax,nspnfv)
real(8), intent(in) :: occsvp(nstsv),doccsvp(nstsv)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv)
complex(8), intent(in) :: apwalmq(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv)
complex(8), intent(in) :: dapwalm(ngkmax,apwordmax,lmmaxapw,nspnfv)
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv)
complex(8), intent(in) :: devecfv(nmatmax,nstfv,nspnfv)
complex(8), intent(in) :: evecsv(nstsv,nstsv),devecsv(nstsv,nstsv)
! local variables
integer nst,ist,jst,is,ias
integer nr,nrci,ir,irc
integer lmmax
real(8) wo,dwo
! automatic arrays
integer idx(nstsv)
! allocatable arrays
complex(8), allocatable :: wfmt(:,:,:,:,:),wfir(:,:,:)
complex(8), allocatable :: dwfmt(:,:,:,:,:),dwfir(:,:,:)
! count and index the occupied states
nst=0
do ist=1,nstsv
  if (abs(occsvp(ist)).gt.epsocc) then
    nst=nst+1
    idx(nst)=ist
  end if
end do
! generate the wavefunctions
allocate(wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nst))
allocate(wfir(ngtot,nspinor,nst))
call genwfsv(.false.,.false.,nst,idx,ngp,igpig,apwalm,evecfv,evecsv,wfmt, &
 ngtot,wfir)
! generate the wavefunction derivatives
allocate(dwfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nst))
allocate(dwfir(ngtot,nspinor,nst))
call gendwfsv(.false.,.false.,nst,idx,ngp,ngpq,igpqig,apwalmq,dapwalm,evecfv, &
 devecfv,evecsv,devecsv,dwfmt,ngtot,dwfir)
! loop over occupied states
do ist=1,nst
  jst=idx(ist)
  wo=2.d0*wkptnr*occsvp(jst)
  dwo=wkptnr*doccsvp(jst)
!----------------------------------------------!
!     muffin-tin density and magnetisation     !
!----------------------------------------------!
  do ias=1,natmtot
    is=idxis(ias)
    nr=nrmt(is)
    nrci=nrcmtinr(is)
!$OMP CRITICAL
    if (spinpol) then
! spin-polarised
      if (ncmag) then
! non-collinear
        lmmax=lmmaxinr
        irc=0
        do ir=1,nr,lradstp
          irc=irc+1
          call drmk1(lmmax,wo,wfmt(:,irc,ias,1,jst),wfmt(:,irc,ias,2,jst), &
           dwfmt(:,irc,ias,1,jst),dwfmt(:,irc,ias,2,jst),drhomt(:,ir,ias), &
           dmagmt(:,ir,ias,1),dmagmt(:,ir,ias,2),dmagmt(:,ir,ias,3))
          if (tphiq0) then
            call drmk01(lmmax,dwo,wfmt(:,irc,ias,1,jst),wfmt(:,irc,ias,2,jst), &
             drhomt(:,ir,ias),dmagmt(:,ir,ias,1),dmagmt(:,ir,ias,2), &
             dmagmt(:,ir,ias,3))
          end if
          if (irc.eq.nrci) lmmax=lmmaxvr
        end do
      else
! collinear
        lmmax=lmmaxinr
        irc=0
        do ir=1,nr,lradstp
          irc=irc+1
          call drmk2(lmmax,wo,wfmt(:,irc,ias,1,jst),wfmt(:,irc,ias,2,jst), &
           dwfmt(:,irc,ias,1,jst),dwfmt(:,irc,ias,2,jst),drhomt(:,ir,ias), &
           dmagmt(:,ir,ias,1))
          if (tphiq0) then
            call drmk02(lmmax,dwo,wfmt(:,irc,ias,1,jst),wfmt(:,irc,ias,2,jst), &
             drhomt(:,ir,ias),dmagmt(:,ir,ias,1))
          end if
          if (irc.eq.nrci) lmmax=lmmaxvr
        end do
      end if
    else
! spin-unpolarised
      lmmax=lmmaxinr
      irc=0
      do ir=1,nr,lradstp
        irc=irc+1
        call drmk3(lmmax,wo,wfmt(:,irc,ias,1,jst),dwfmt(:,irc,ias,1,jst), &
         drhomt(:,ir,ias))
        if (tphiq0) then
          call drmk03(lmmax,dwo,wfmt(:,irc,ias,1,jst),drhomt(:,ir,ias))
        end if
        if (irc.eq.nrci) lmmax=lmmaxvr
      end do
    end if
!$OMP END CRITICAL
  end do
!------------------------------------------------!
!     interstitial density and magnetisation     !
!------------------------------------------------!
!$OMP CRITICAL
  if (spinpol) then
! spin-polarised
    if (ncmag) then
      call drmk1(ngtot,wo,wfir(:,1,jst),wfir(:,2,jst),dwfir(:,1,jst), &
       dwfir(:,2,jst),drhoir,dmagir(:,1),dmagir(:,2),dmagir(:,3))
      if (tphiq0) then
        call drmk01(ngtot,dwo,wfir(:,1,jst),wfir(:,2,jst),drhoir,dmagir(:,1), &
         dmagir(:,2),dmagir(:,3))
      end if
    else
! collinear
      call drmk2(ngtot,wo,wfir(:,1,jst),wfir(:,2,jst),dwfir(:,1,jst), &
       dwfir(:,2,jst),drhoir,dmagir)
      if (tphiq0) then
        call drmk02(ngtot,dwo,wfir(:,1,jst),wfir(:,2,jst),drhoir,dmagir)
      end if
    end if
  else
! spin-unpolarised
    call drmk3(ngtot,wo,wfir(:,1,jst),dwfir(:,1,jst),drhoir)
    if (tphiq0) then
      call drmk03(ngtot,dwo,wfir(:,1,jst),drhoir)
    end if
  end if
!$OMP END CRITICAL
! end loop over states
end do
deallocate(wfmt,wfir,dwfmt,dwfir)
return

contains

subroutine drmk1(n,wo,wf1,wf2,dwf1,dwf2,drho,dmag1,dmag2,dmag3)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: wo
complex(8), intent(in) :: wf1(n),wf2(n)
complex(8), intent(in) :: dwf1(n),dwf2(n)
complex(8), intent(inout) :: drho(n)
complex(8), intent(inout) :: dmag1(n),dmag2(n),dmag3(n)
! local variables
integer i
complex(8) z1,z2,z3,z4,z5,z6
do i=1,n
  z1=conjg(wf1(i))
  z2=conjg(wf2(i))
  z3=dwf1(i)
  z4=dwf2(i)
  z5=z1*z3
  z6=z2*z4
  drho(i)=drho(i)+wo*(z5+z6)
  dmag3(i)=dmag3(i)+wo*(z5-z6)
  z5=z1*z4
  z6=z2*z3
  dmag1(i)=dmag1(i)+wo*(z5+z6)
  z5=z5-z6
  dmag2(i)=dmag2(i)+wo*cmplx(aimag(z5),-dble(z5),8)
end do
return
end subroutine

subroutine drmk01(n,dwo,wf1,wf2,drho,dmag1,dmag2,dmag3)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: dwo
complex(8), intent(in) :: wf1(n),wf2(n)
complex(8), intent(inout) :: drho(n)
complex(8), intent(inout) :: dmag1(n),dmag2(n),dmag3(n)
! local variables
integer i
real(8) t1,t2
complex(8) z1,z2
do i=1,n
  z1=wf1(i)
  z2=wf2(i)
  t1=dble(z1)**2+aimag(z1)**2
  t2=dble(z2)**2+aimag(z2)**2
  z1=conjg(z1)*z2
  drho(i)=drho(i)+dwo*(t1+t2)
  dmag1(i)=dmag1(i)+dwo*2.d0*dble(z1)
  dmag2(i)=dmag2(i)+dwo*2.d0*aimag(z1)
  dmag3(i)=dmag3(i)+dwo*(t1-t2)
end do
return
end subroutine

subroutine drmk2(n,wo,wf1,wf2,dwf1,dwf2,drho,dmag)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: wo
complex(8), intent(in) :: wf1(n),wf2(n)
complex(8), intent(in) :: dwf1(n),dwf2(n)
complex(8), intent(inout) :: drho(n),dmag(n)
! local variables
integer i
complex(8) z1,z2
do i=1,n
  z1=conjg(wf1(i))*dwf1(i)
  z2=conjg(wf2(i))*dwf2(i)
  drho(i)=drho(i)+wo*(z1+z2)
  dmag(i)=dmag(i)+wo*(z1-z2)
end do
return
end subroutine

subroutine drmk02(n,dwo,wf1,wf2,drho,dmag)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: dwo
complex(8), intent(in) :: wf1(n),wf2(n)
complex(8), intent(inout) :: drho(n),dmag(n)
! local variables
integer i
real(8) t1,t2
do i=1,n
  t1=dble(wf1(i))**2+aimag(wf1(i))**2
  t2=dble(wf2(i))**2+aimag(wf2(i))**2
  drho(i)=drho(i)+dwo*(t1+t2)
  dmag(i)=dmag(i)+dwo*(t1-t2)
end do
return
end subroutine

subroutine drmk3(n,wo,wf,dwf,drho)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: wo
complex(8), intent(in) :: wf(n),dwf(n)
complex(8), intent(inout) :: drho(n)
drho(:)=drho(:)+wo*conjg(wf(:))*dwf(:)
return
end subroutine

subroutine drmk03(n,dwo,wf,drho)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: dwo
complex(8), intent(in) :: wf(n)
complex(8), intent(inout) :: drho(n)
drho(:)=drho(:)+dwo*(dble(wf(:))**2+aimag(wf(:))**2)
return
end subroutine

end subroutine

