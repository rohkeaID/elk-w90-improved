
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gradwf2(ik,evecfv,evecsv,gwf2mt,gwf2ir)
use modmain
implicit none
! arguments
integer, intent(in) :: ik
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv)
real(8), intent(inout) :: gwf2mt(npmtmax,natmtot),gwf2ir(ngtot)
! local variables
integer ist,ispn,jspn
integer is,ia,ias
integer nr,nri,ir,np
integer igk,ifg,i,j
real(8) wo,t1
complex(8) zq(2),z1
! automatic arrays
logical done(nstfv,nspnfv)
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:,:)
complex(8), allocatable :: wfmt1(:,:,:),wfmt2(:,:)
complex(8), allocatable :: gwfmt(:,:),zfmt(:)
complex(8), allocatable :: zfft1(:,:),zfft2(:)
!-------------------------!
!     muffin-tin part     !
!-------------------------!
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
allocate(wfmt1(npmtmax,nstfv,nspnfv),wfmt2(npmtmax,nspinor))
allocate(gwfmt(npmtmax,3),zfmt(npmtmax))
! find the matching coefficients
do ispn=1,nspnfv
  call match(ngk(ispn,ik),gkc(:,ispn,ik),tpgkc(:,:,ispn,ik), &
   sfacgk(:,:,ispn,ik),apwalm(:,:,:,:,ispn))
end do
do is=1,nspecies
  nr=nrmt(is)
  nri=nrmti(is)
  np=npmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! de-phasing factor for spin-spirals
    if (spinsprl.and.ssdph) then
      t1=-0.5d0*dot_product(vqcss(:),atposc(:,ia,is))
      zq(1)=cmplx(cos(t1),sin(t1),8)
      zq(2)=conjg(zq(1))
    end if
    done(:,:)=.false.
    do j=1,nstsv
      if (abs(occsv(j,ik)).lt.epsocc) cycle
      wo=wkpt(ik)*occsv(j,ik)
      if (tevecsv) then
! generate spinor wavefunction from second-variational eigenvectors
        wfmt2(1:np,:)=0.d0
        i=0
        do ispn=1,nspinor
          jspn=jspnfv(ispn)
          do ist=1,nstfv
            i=i+1
            z1=evecsv(i,j)
            if (spinsprl.and.ssdph) z1=z1*zq(ispn)
            if (abs(dble(z1))+abs(aimag(z1)).gt.epsocc) then
              if (.not.done(ist,jspn)) then
                call wavefmt(1,ias,ngk(jspn,ik),apwalm(:,:,:,ias,jspn), &
                 evecfv(:,ist,jspn),wfmt1(:,ist,jspn))
                done(ist,jspn)=.true.
              end if
! add to spinor wavefunction
              wfmt2(1:np,ispn)=wfmt2(1:np,ispn)+z1*wfmt1(1:np,ist,jspn)
            end if
          end do
        end do
      else
! spin-unpolarised wavefunction
        call wavefmt(1,ias,ngk(1,ik),apwalm(:,:,:,ias,1),evecfv(:,j,1),wfmt2)
      end if
! compute the gradient of the wavefunction
      do ispn=1,nspinor
        call gradzfmt(nr,nri,rsp(:,is),wfmt2(:,ispn),npmtmax,gwfmt)
        do i=1,3
! convert gradient from spherical harmonics to spherical coordinates
          call zbsht(nr,nri,gwfmt(:,i),zfmt)
! add to total
          gwf2mt(1:np,ias)=gwf2mt(1:np,ias) &
           +wo*(dble(zfmt(1:np))**2+aimag(zfmt(1:np))**2)
        end do
      end do
    end do
  end do
end do
deallocate(apwalm,wfmt1,wfmt2,gwfmt,zfmt)
!---------------------------!
!     interstitial part     !
!---------------------------!
allocate(zfft1(ngtot,nspinor),zfft2(ngtot))
do j=1,nstsv
  if (abs(occsv(j,ik)).lt.epsocc) cycle
  wo=wkpt(ik)*occsv(j,ik)/omega
  zfft1(:,:)=0.d0
  if (tevecsv) then
! generate spinor wavefunction from second-variational eigenvectors
    i=0
    do ispn=1,nspinor
      jspn=jspnfv(ispn)
      do ist=1,nstfv
        i=i+1
        z1=evecsv(i,j)
        if (abs(dble(z1))+abs(aimag(z1)).gt.epsocc) then
          do igk=1,ngk(jspn,ik)
            ifg=igfft(igkig(igk,jspn,ik))
            zfft1(ifg,ispn)=zfft1(ifg,ispn)+z1*evecfv(igk,ist,jspn)
          end do
        end if
      end do
    end do
  else
! spin-unpolarised wavefunction
    do igk=1,ngk(1,ik)
      ifg=igfft(igkig(igk,1,ik))
      zfft1(ifg,1)=evecfv(igk,j,1)
    end do
  end if
! compute gradient of wavefunction
  do ispn=1,nspinor
    jspn=jspnfv(ispn)
    do i=1,3
      zfft2(:)=0.d0
      do igk=1,ngk(jspn,ik)
        ifg=igfft(igkig(igk,jspn,ik))
        zfft2(ifg)=zi*vgkc(i,igk,jspn,ik)*zfft1(ifg,ispn)
      end do
! Fourier transform gradient to real-space
      call zfftifc(3,ngridg,1,zfft2)
      do ir=1,ngtot
        gwf2ir(ir)=gwf2ir(ir)+wo*(dble(zfft2(ir))**2+aimag(zfft2(ir))**2)
      end do
    end do
  end do
end do
deallocate(zfft1,zfft2)
return
end subroutine

