
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genwfpw(vpl,ngp,igpig,vgpl,vgpc,gpc,tpgpc,sfacgp,nhp,vhpc,hpc, &
 tphpc,sfachp,wfpw)
use modmain
use modpw
implicit none
! arguments
real(8), intent(in) :: vpl(3)
integer, intent(in) :: ngp(nspnfv)
integer, intent(in) :: igpig(ngkmax,nspnfv)
real(8), intent(in) :: vgpl(3,ngkmax,nspnfv)
real(8), intent(in) :: vgpc(3,ngkmax,nspnfv)
real(8), intent(in) :: gpc(ngkmax,nspnfv)
real(8), intent(in) :: tpgpc(2,ngkmax,nspnfv)
complex(8), intent(in) :: sfacgp(ngkmax,natmtot,nspnfv)
integer, intent(in) :: nhp(nspnfv)
real(8), intent(in) :: vhpc(3,nhkmax,nspnfv)
real(8), intent(in) :: hpc(nhkmax,nspnfv)
real(8), intent(in) :: tphpc(2,nhkmax,nspnfv)
complex(8), intent(in) :: sfachp(nhkmax,natmtot,nspnfv)
complex(8), intent(out) :: wfpw(nhkmax,nspinor,nstsv)
! local variables
integer ispn0,ispn1,ispn,jspn
integer ist,is,ia,ias,i
integer nrc,nrci,iro,irc
integer lmax,l,m,lm,igp,ihp
real(8) t0,t1,t2
complex(8) zsum1,zsum2
complex(8) z1,z2,z3,z4
! automatic arrays
integer idx(nstsv)
real(8) fr1(nrcmtmax),fr2(nrcmtmax),gr(nrcmtmax)
complex(8) ylm(lmmaxvr)
! allocatable arrays
real(8), allocatable :: jl(:,:)
complex(8), allocatable :: apwalm(:,:,:,:,:)
complex(8), allocatable :: evecfv(:,:,:),evecsv(:,:)
complex(8), allocatable :: wfmt(:,:,:,:,:),wfir(:,:,:)
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
allocate(evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv))
allocate(wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfir(ngkmax,nspinor,nstsv))
! get the eigenvectors from file
call getevecfv(filext,vpl,vgpl,evecfv)
call getevecsv(filext,vpl,evecsv)
! find the matching coefficients
do ispn=1,nspnfv
  call match(ngp(ispn),gpc(:,ispn),tpgpc(:,:,ispn),sfacgp(:,:,ispn), &
   apwalm(:,:,:,:,ispn))
end do
! index to all states
do ist=1,nstsv
  idx(ist)=ist
end do
! calculate the second-variational wavefunctions for all states
call genwfsv(.true.,.true.,nstsv,idx,ngp,igpig,apwalm,evecfv,evecsv,wfmt, &
 ngkmax,wfir)
deallocate(apwalm,evecfv,evecsv)
! zero the plane wave coefficients
wfpw(:,:,:)=0.d0
!---------------------------!
!     interstitial part     !
!---------------------------!
do jspn=1,nspnfv
  if (spinsprl) then
    ispn0=jspn; ispn1=jspn
  else
    ispn0=1; ispn1=nspinor
  end if
  i=1
  do ihp=1,nhp(jspn)
    do igp=i,ngp(jspn)
      t1=abs(vhpc(1,ihp,jspn)-vgpc(1,igp,jspn)) &
        +abs(vhpc(2,ihp,jspn)-vgpc(2,igp,jspn)) &
        +abs(vhpc(3,ihp,jspn)-vgpc(3,igp,jspn))
      if (t1.lt.epslat) then
        do ist=1,nstsv
          do ispn=ispn0,ispn1
            wfpw(ihp,ispn,ist)=wfir(igp,ispn,ist)
          end do
        end do
        if (igp.eq.i) i=i+1
        exit
      end if
    end do
  end do
end do
!-------------------------!
!     muffin-tin part     !
!-------------------------!
allocate(jl(0:lmaxvr,nrcmtmax))
t0=fourpi/sqrt(omega)
! remove continuation of interstitial function into muffin-tin
do jspn=1,nspnfv
  if (spinsprl) then
    ispn0=jspn; ispn1=jspn
  else
    ispn0=1; ispn1=nspinor
  end if
! loop over G+p-vectors
  do igp=1,ngp(jspn)
! generate the conjugate spherical harmonics Y_lm*(G+p)
    call genylm(lmaxvr,tpgpc(:,igp,jspn),ylm)
    ylm(:)=conjg(ylm(:))
! loop over species
    do is=1,nspecies
      nrc=nrcmt(is)
      nrci=nrcmtinr(is)
! generate spherical Bessel functions
      do irc=1,nrc
        t1=gpc(igp,jspn)*rcmt(irc,is)
        call sbessel(lmaxvr,t1,jl(:,irc))
      end do
! loop over atoms
      do ia=1,natoms(is)
        ias=idxas(ia,is)
        z1=t0*sfacgp(igp,ias,jspn)
        do ist=1,nstsv
          do ispn=ispn0,ispn1
            z2=z1*wfir(igp,ispn,ist)
            lm=0
            do l=0,lmaxvr
              if (l.le.lmaxinr) then
                iro=1
              else
                iro=nrci+1
              end if
              z3=z2*zil(l)
              do m=-l,l
                lm=lm+1
                z4=z3*ylm(lm)
                wfmt(lm,iro:nrc,ias,ispn,ist)=wfmt(lm,iro:nrc,ias,ispn,ist) &
                 -z4*jl(l,iro:nrc)
              end do
            end do
          end do
        end do
      end do
    end do
  end do
end do
! Fourier transform the muffin-tin wavefunctions
do jspn=1,nspnfv
  if (spinsprl) then
    ispn0=jspn; ispn1=jspn
  else
    ispn0=1; ispn1=nspinor
  end if
! loop over H+p-vectors
  do ihp=1,nhp(jspn)
! generate the spherical harmonics Y_lm(H+p)
    call genylm(lmaxvr,tphpc(:,ihp,jspn),ylm)
    do is=1,nspecies
      nrc=nrcmt(is)
      nrci=nrcmtinr(is)
! generate spherical Bessel functions
      do irc=1,nrc
        t1=hpc(ihp,jspn)*rcmt(irc,is)
        call sbessel(lmaxvr,t1,jl(:,irc))
      end do
      do ia=1,natoms(is)
        ias=idxas(ia,is)
! conjugate structure factor
        z1=t0*conjg(sfachp(ihp,ias,jspn))
! loop over states
        do ist=1,nstsv
          do ispn=ispn0,ispn1
            do irc=1,nrc
              if (irc.le.nrci) then
                lmax=lmaxinr
              else
                lmax=lmaxvr
              end if
              zsum1=jl(0,irc)*wfmt(1,irc,ias,ispn,ist)*ylm(1)
              lm=1
              do l=1,lmax
                lm=lm+1
                zsum2=wfmt(lm,irc,ias,ispn,ist)*ylm(lm)
                do m=1-l,l
                  lm=lm+1
                  zsum2=zsum2+wfmt(lm,irc,ias,ispn,ist)*ylm(lm)
                end do
                zsum1=zsum1+jl(l,irc)*zilc(l)*zsum2
              end do
              zsum1=zsum1*rcmt(irc,is)**2
              fr1(irc)=dble(zsum1)
              fr2(irc)=aimag(zsum1)
            end do
            call fderiv(-1,nrc,rcmt(:,is),fr1,gr)
            t1=gr(nrc)
            call fderiv(-1,nrc,rcmt(:,is),fr2,gr)
            t2=gr(nrc)
! add to the H+p wavefunction
            wfpw(ihp,ispn,ist)=wfpw(ihp,ispn,ist)+z1*cmplx(t1,t2,8)
          end do
        end do
      end do
    end do
  end do
end do
deallocate(jl,wfmt,wfir)
return
end subroutine

