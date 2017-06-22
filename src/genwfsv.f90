
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genwfsv
! !INTERFACE:
subroutine genwfsv(tsh,tgp,nst,idx,ngp,igpig,apwalm,evecfv,evecsv,wfmt,ld,wfir)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   tsh    : .true. if wfmt should be in spherical harmonic basis (in,logical)
!   tgp    : .true. if wfir should be in G+p-space, otherwise in real-space
!            (in,logical)
!   nst    : number of states to be calculated (in,integer)
!   idx    : index to states which are to be calculated (in,integer(nst))
!   ngp    : number of G+p-vectors (in,integer(nspnfv))
!   igpig  : index from G+p-vectors to G-vectors (in,integer(ngkmax,nspnfv))
!   apwalm : APW matching coefficients
!            (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
!   evecfv : first-variational eigenvectors (in,complex(nmatmax,nstfv,nspnfv))
!   evecsv : second-variational eigenvectors (in,complex(nstsv,nstsv))
!   wfmt   : muffin-tin part of the wavefunctions for every state in spherical
!            coordinates (out,complex(lmmaxvr,nrcmtmax,natmtot,nspinor,nst))
!   ld     : leading dimension, at least ngkmax if tgp is .true. and ngtot if
!            tgp is .false. (in,integer)
!   wfir   : interstitial part of the wavefunctions for every state
!            (out,complex(ld,nspinor,nst))
! !DESCRIPTION:
!   Calculates the second-variational spinor wavefunctions in both the
!   muffin-tin and interstitial regions for every state of a particular
!   $k$-point. The wavefunctions in both regions are stored on a real-space
!   grid. A coarse radial mesh is assumed in the muffin-tins with with angular
!   momentum cut-off of {\tt lmaxvr}. If {\tt tocc} is {\tt .true.}, then only
!   those states with occupancy greater than {\tt epsocc} are calculated.
!
! !REVISION HISTORY:
!   Created November 2004 (Sharma)
!   Updated for spin-spirals, June 2010 (JKD)
!EOP
!BOC
implicit none
! arguments
logical, intent(in) :: tsh,tgp
integer, intent(in) :: nst,idx(nst)
integer, intent(in) :: ngp(nspnfv),igpig(ngkmax,nspnfv)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv)
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv)
complex(8), intent(out) :: wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nst)
integer, intent(in) :: ld
complex(8), intent(out) :: wfir(ld,nspinor,nst)
! local variables
integer ist,ispn,jspn
integer is,ia,ias,nrc,nrci
integer igp,ifg,i,j,k
real(8) t1
complex(8) zq(2),z1
! automatic arrays
logical done(nstfv,nspnfv)
! allocatable arrays
complex(8), allocatable :: wfmt1(:,:,:,:),wfmt2(:,:)
!--------------------------------!
!     muffin-tin wavefunction    !
!--------------------------------!
if (tevecsv) allocate(wfmt1(lmmaxvr,nrcmtmax,nstfv,nspnfv))
if (.not.tsh) allocate(wfmt2(lmmaxvr,nrcmtmax))
do is=1,nspecies
  nrc=nrcmt(is)
  nrci=nrcmtinr(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! de-phasing factor for spin-spirals
    if (spinsprl.and.ssdph) then
      t1=-0.5d0*dot_product(vqcss(:),atposc(:,ia,is))
      zq(1)=cmplx(cos(t1),sin(t1),8)
      zq(2)=conjg(zq(1))
    end if
    done(:,:)=.false.
! loop only over required states
    do j=1,nst
! index to state in evecsv
      k=idx(j)
      if (tevecsv) then
! generate spinor wavefunction from second-variational eigenvectors
        i=0
        do ispn=1,nspinor
          jspn=jspnfv(ispn)
          if (tsh) then
            wfmt(:,:,ias,ispn,j)=0.d0
          else
            wfmt2(:,:)=0.d0
          end if
          do ist=1,nstfv
            i=i+1
            z1=evecsv(i,k)
            if (spinsprl.and.ssdph) z1=z1*zq(ispn)
            if (abs(dble(z1))+abs(aimag(z1)).gt.epsocc) then
              if (.not.done(ist,jspn)) then
                call wavefmt(lradstp,ias,ngp(jspn),apwalm(:,:,:,:,jspn), &
                 evecfv(:,ist,jspn),wfmt1(:,:,ist,jspn))
                done(ist,jspn)=.true.
              end if
! add to spinor wavefunction
              if (tsh) then
                call zfmtadd(nrc,nrci,z1,wfmt1(:,:,ist,jspn), &
                 wfmt(:,:,ias,ispn,j))
              else
                call zfmtadd(nrc,nrci,z1,wfmt1(:,:,ist,jspn),wfmt2)
              end if
            end if
          end do
! convert to spherical coordinates if required
          if (.not.tsh) then
            call zbsht(nrc,nrci,wfmt2,wfmt(:,:,ias,ispn,j))
          end if
        end do
      else
! spin-unpolarised wavefunction
        if (tsh) then
          call wavefmt(lradstp,ias,ngp,apwalm,evecfv(:,k,1),wfmt(:,:,ias,1,j))
        else
          call wavefmt(lradstp,ias,ngp,apwalm,evecfv(:,k,1),wfmt2)
! convert from to spherical coordinates
          call zbsht(nrc,nrci,wfmt2,wfmt(:,:,ias,1,j))
        end if
      end if
! end loop over second-variational states
    end do
! end loops over atoms and species
  end do
end do
if (tevecsv) deallocate(wfmt1)
if (.not.tsh) deallocate(wfmt2)
!-----------------------------------!
!     interstitial wavefunction     !
!-----------------------------------!
t1=1.d0/sqrt(omega)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(i,k,ispn,jspn,ist) &
!$OMP PRIVATE(z1,igp,ifg)
!$OMP DO
do j=1,nst
  k=idx(j)
  wfir(:,:,j)=0.d0
  if (tevecsv) then
! generate spinor wavefunction from second-variational eigenvectors
    i=0
    do ispn=1,nspinor
      jspn=jspnfv(ispn)
      do ist=1,nstfv
        i=i+1
        z1=evecsv(i,k)
        if (abs(dble(z1))+abs(aimag(z1)).gt.epsocc) then
          if (tgp) then
! wavefunction in G+p-space
            call zaxpy(ngp(jspn),z1,evecfv(:,ist,jspn),1,wfir(:,ispn,j),1)
          else
! wavefunction in real-space
            z1=t1*z1
            do igp=1,ngp(jspn)
              ifg=igfft(igpig(igp,jspn))
              wfir(ifg,ispn,j)=wfir(ifg,ispn,j)+z1*evecfv(igp,ist,jspn)
            end do
          end if
        end if
      end do
    end do
  else
! spin-unpolarised wavefunction
    if (tgp) then
      call zcopy(ngp(1),evecfv(:,k,1),1,wfir(:,1,j),1)
    else
      do igp=1,ngp(1)
        ifg=igfft(igpig(igp,1))
        wfir(ifg,1,j)=t1*evecfv(igp,k,1)
      end do
    end if
  end if
! Fourier transform wavefunction to real-space if required
  if (.not.tgp) then
    do ispn=1,nspinor
      call zfftifc(3,ngridg,1,wfir(:,ispn,j))
    end do
  end if
end do
!$OMP END DO
!$OMP END PARALLEL
return
end subroutine
!EOC

