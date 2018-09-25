
! Copyright (C) 2017 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rhomaguk(ik0,evecu)
use modmain
use modulr
implicit none
! arguments
integer, intent(in) :: ik0
complex(8), intent(in) :: evecu(nstulr,nstulr)
! local variables
integer ik,ist,ikpa
integer ngk0,igk,ifg
integer ispn,is,ias
integer npc,ir,i,j
real(8) wo
! automatic arrays
integer idx(nstsv)
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:),evecfv(:,:),evecsv(:,:)
complex(8), allocatable :: wfmt(:,:,:,:),wfgk(:,:,:)
complex(8), allocatable :: wfmtu(:,:,:,:),wfgku(:,:,:)
complex(8), allocatable :: wfir(:,:),zfft(:)
! central k-point
ik=(ik0-1)*nkpa+1
ngk0=ngk(1,ik)
! get the eigenvectors from file
allocate(evecfv(nmatmax,nstfv),evecsv(nstsv,nstsv))
call getevecfv(filext,ik,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
call getevecsv(filext,ik,vkl(:,ik),evecsv)
! find the matching coefficients
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
call match(ngk0,gkc(:,1,ik),tpgkc(:,:,1,ik),sfacgk(:,:,1,ik),apwalm)
! index to all states
do ist=1,nstsv
  idx(ist)=ist
end do
allocate(wfmt(npcmtmax,natmtot,nspinor,nstsv),wfgk(ngkmax,nspinor,nstsv))
call genwfsv(.false.,.true.,nstsv,idx,ngridg,igfft,ngk(1,ik),igkig(:,1,ik), &
 apwalm,evecfv,evecsv,wfmt,ngkmax,wfgk)
deallocate(apwalm,evecfv,evecsv)
allocate(wfmtu(npcmtmax,natmtot,nspinor,nqpt),wfgku(ngkmax,nspinor,nqpt))
allocate(wfir(ngtot,nspinor),zfft(nqpt))
! loop over ultra long-range states
do j=1,nstulr
  wo=wkpt(ik0)*occulr(j,ik0)
  if (abs(wo).lt.epsocc) cycle
! zero the ultra long-range wavefunctions
  do ir=1,nqpt
    do ispn=1,nspinor
      call wfmt0(wfmtu(:,:,ispn,ir))
    end do
  end do
  wfgku(:,:,:)=0.d0
  i=0
! loop over second-variational states
  do ist=1,nstsv
    zfft(1:nqpt)=0.d0
! loop over kappa-points
    do ikpa=1,nkpa
      i=i+1
! store the wavefunction in Q-space
      zfft(iqfft(ikpa))=evecu(i,j)
    end do
! Fourier transform to R-space
    call zfftifc(3,ngridq,1,zfft)
! loop over R-points
    do ir=1,nqpt
      call wfadd(ngk0,zfft(ir),wfmt(:,:,:,ist),wfgk(:,:,ist),wfmtu(:,:,:,ir), &
       wfgku(:,:,ir))
    end do
  end do
! loop over R-points
  do ir=1,nqpt
    do ispn=1,nspinor
! Fourier transform the interstitial part to real-space
      wfir(1:ngtot,ispn)=0.d0
      do igk=1,ngk0
        ifg=igfft(igkig(igk,1,ik))
        wfir(ifg,ispn)=wfgku(igk,ispn,ir)
      end do
      call zfftifc(3,ngridg,1,wfir(:,ispn))
    end do
! add to the density and magnetisation
    do ias=1,natmtot
      is=idxis(ias)
      npc=npcmt(is)
      if (spinpol) then
        if (ncmag) then
          call rmk1(npc,wo,wfmtu(:,ias,1,ir),wfmtu(:,ias,2,ir), &
           rhormt(:,ias,ir),magrmt(:,ias,1,ir),magrmt(:,ias,2,ir), &
           magrmt(:,ias,3,ir))
        else
          call rmk2(npc,wo,wfmtu(:,ias,1,ir),wfmtu(:,ias,2,ir), &
           rhormt(:,ias,ir),magrmt(:,ias,1,ir))
        end if
      else
        call rmk3(npc,wo,wfmtu(:,ias,1,ir),rhormt(:,ias,ir))
      end if
    end do
    if (spinpol) then
      if (ncmag) then
        call rmk1(ngtot,wo,wfir(:,1),wfir(:,2),rhorir(:,ir),magrir(:,1,ir), &
         magrir(:,2,ir),magrir(:,3,ir))
      else
        call rmk2(ngtot,wo,wfir(:,1),wfir(:,2),rhorir(:,ir),magrir(:,1,ir))
      end if
    else
      call rmk3(ngtot,wo,wfir(:,1),rhorir(:,ir))
    end if
! end loop over R-points
  end do
! end loop over long-range states
end do
deallocate(wfmt,wfgk,wfmtu,wfgku,wfir,zfft)
return

contains

subroutine wfmt0(wfmt)
implicit none
! arguments
complex(8), intent(out) :: wfmt(npcmtmax,natmtot)
! local variables
integer is,ias
do ias=1,natmtot
  is=idxis(ias)
  wfmt(1:npcmt(is),ias)=0.d0
end do
return
end subroutine

subroutine wfadd(ngp,za,wfmt1,wfir1,wfmt2,wfir2)
implicit none
! arguments
integer, intent(in) :: ngp
complex(8), intent(in) :: za
complex(8), intent(in) :: wfmt1(npcmtmax,natmtot,nspinor)
complex(8), intent(in) :: wfir1(ngkmax,nspinor)
complex(8), intent(inout) :: wfmt2(npcmtmax,natmtot,nspinor)
complex(8), intent(inout) :: wfir2(ngkmax,nspinor)
! local variables
integer ispn,is,ias
do ispn=1,nspinor
  do ias=1,natmtot
    is=idxis(ias)
    call zaxpy(npcmt(is),za,wfmt1(:,ias,ispn),1,wfmt2(:,ias,ispn),1)
  end do
end do
do ispn=1,nspinor
  call zaxpy(ngp,za,wfir1(:,ispn),1,wfir2(:,ispn),1)
end do
return
end subroutine

subroutine rmk1(n,wo,wf1,wf2,rho,mag1,mag2,mag3)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: wo
complex(8), intent(in) :: wf1(n),wf2(n)
real(8), intent(inout) :: rho(n),mag1(n),mag2(n),mag3(n)
! local variables
integer i
real(8) wo2,t1,t2
complex(8) z1,z2
wo2=2.d0*wo
do i=1,n
  z1=wf1(i)
  z2=wf2(i)
  t1=dble(z1)**2+aimag(z1)**2
  t2=dble(z2)**2+aimag(z2)**2
  z1=conjg(z1)*z2
  rho(i)=rho(i)+wo*(t1+t2)
  mag1(i)=mag1(i)+wo2*dble(z1)
  mag2(i)=mag2(i)+wo2*aimag(z1)
  mag3(i)=mag3(i)+wo*(t1-t2)
end do
return
end subroutine

subroutine rmk2(n,wo,wf1,wf2,rho,mag)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: wo
complex(8), intent(in) :: wf1(n),wf2(n)
real(8), intent(inout) :: rho(n),mag(n)
! local variables
integer i
real(8) t1,t2
do i=1,n
  t1=dble(wf1(i))**2+aimag(wf1(i))**2
  t2=dble(wf2(i))**2+aimag(wf2(i))**2
  rho(i)=rho(i)+wo*(t1+t2)
  mag(i)=mag(i)+wo*(t1-t2)
end do
return
end subroutine

subroutine rmk3(n,wo,wf,rho)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: wo
complex(8), intent(in) :: wf(n)
real(8), intent(inout) :: rho(n)
rho(:)=rho(:)+wo*(dble(wf(:))**2+aimag(wf(:))**2)
return
end subroutine

end subroutine

