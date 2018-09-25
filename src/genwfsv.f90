
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genwfsv
! !INTERFACE:
subroutine genwfsv(tsh,tgp,nst,idx,ngdg,igf,ngp,igpig,apwalm,evecfv,evecsv, &
 wfmt,ld,wfir)
! !USES:
use modmain
use modomp
! !INPUT/OUTPUT PARAMETERS:
!   tsh    : .true. if wfmt should be in spherical harmonic basis (in,logical)
!   tgp    : .true. if wfir should be in G+p-space, otherwise in real-space
!            (in,logical)
!   nst    : number of states to be calculated (in,integer)
!   idx    : index to states which are to be calculated (in,integer(nst))
!   ngdg   : G-vector grid sizes (in,integer(3))
!   igf    : map from G-vector index to FFT array (in,integer(*))
!   ngp    : number of G+p-vectors (in,integer(nspnfv))
!   igpig  : index from G+p-vectors to G-vectors (in,integer(ngkmax,nspnfv))
!   apwalm : APW matching coefficients
!            (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
!   evecfv : first-variational eigenvectors (in,complex(nmatmax,nstfv,nspnfv))
!   evecsv : second-variational eigenvectors (in,complex(nstsv,nstsv))
!   wfmt   : muffin-tin part of the wavefunctions for every state in spherical
!            coordinates (out,complex(npcmtmax,natmtot,nspinor,nst))
!   ld     : leading dimension of wfir (in,integer)
!   wfir   : interstitial part of the wavefunctions for every state
!            (out,complex(ld,nspinor,nst))
! !DESCRIPTION:
!   Calculates the second-variational spinor wavefunctions in both the
!   muffin-tin and interstitial regions for every state of a particular
!   $k$-point. A coarse radial mesh is assumed in the muffin-tins with angular
!   momentum cut-off of {\tt lmaxo}.
!
! !REVISION HISTORY:
!   Created November 2004 (Sharma)
!   Updated for spin-spirals, June 2010 (JKD)
!   Packed muffin-tins, April 2016 (JKD)
!EOP
!BOC
implicit none
! arguments
logical, intent(in) :: tsh,tgp
integer, intent(in) :: nst,idx(nst),ngdg(3),igf(*)
integer, intent(in) :: ngp(nspnfv),igpig(ngkmax,nspnfv)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv)
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv)
complex(8), intent(out) :: wfmt(npcmtmax,natmtot,nspinor,nst)
integer, intent(in) :: ld
complex(8), intent(out) :: wfir(ld,nspinor,nst)
! local variables
integer ist,ispn,jspn
integer is,ia,ias,i,j,k
integer nrc,nrci,npc
integer igp,ifg,nthd
real(8) t1
complex(8) zq(2),z1
! automatic arrays
logical done(nstfv,nspnfv)
! allocatable arrays
complex(8), allocatable :: wfmt1(:,:,:),wfmt2(:)
!--------------------------------!
!     muffin-tin wavefunction    !
!--------------------------------!
if (tevecsv) allocate(wfmt1(npcmtmax,nstfv,nspnfv))
if (.not.tsh) allocate(wfmt2(npcmtmax))
do is=1,nspecies
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  npc=npcmt(is)
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
          wfmt(1:npc,ias,ispn,j)=0.d0
          do ist=1,nstfv
            i=i+1
            z1=evecsv(i,k)
            if (abs(dble(z1))+abs(aimag(z1)).gt.epsocc) then
              if (spinsprl.and.ssdph) z1=z1*zq(ispn)
              if (.not.done(ist,jspn)) then
                if (tsh) then
                  call wavefmt(lradstp,ias,ngp(jspn),apwalm(:,:,:,ias,jspn), &
                   evecfv(:,ist,jspn),wfmt1(:,ist,jspn))
                else
                  call wavefmt(lradstp,ias,ngp(jspn),apwalm(:,:,:,ias,jspn), &
                   evecfv(:,ist,jspn),wfmt2)
! convert to spherical coordinates
                  call zbsht(nrc,nrci,wfmt2,wfmt1(:,ist,jspn))
                end if
                done(ist,jspn)=.true.
              end if
! add to spinor wavefunction
              call zaxpy(npc,z1,wfmt1(:,ist,jspn),1,wfmt(:,ias,ispn,j),1)
            end if
          end do
        end do
      else
! spin-unpolarised wavefunction
        if (tsh) then
          call wavefmt(lradstp,ias,ngp,apwalm(:,:,:,ias,1),evecfv(:,k,1), &
           wfmt(:,ias,1,j))
        else
          call wavefmt(lradstp,ias,ngp,apwalm(:,:,:,ias,1),evecfv(:,k,1),wfmt2)
! convert to spherical coordinates
          call zbsht(nrc,nrci,wfmt2,wfmt(:,ias,1,j))
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
call omp_hold(nst,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(i,k,ispn,jspn,ist) &
!$OMP PRIVATE(z1,igp,ifg) &
!$OMP NUM_THREADS(nthd)
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
              ifg=igf(igpig(igp,jspn))
              wfir(ifg,ispn,j)=wfir(ifg,ispn,j)+z1*evecfv(igp,ist,jspn)
            end do
          end if
        end if
      end do
! Fourier transform wavefunction to real-space if required
      if (.not.tgp) call zfftifc(3,ngdg,1,wfir(:,ispn,j))
    end do
  else
! spin-unpolarised wavefunction
    if (tgp) then
      call zcopy(ngp(1),evecfv(:,k,1),1,wfir(:,1,j),1)
    else
      do igp=1,ngp(1)
        ifg=igf(igpig(igp,1))
        wfir(ifg,1,j)=t1*evecfv(igp,k,1)
      end do
      call zfftifc(3,ngdg,1,wfir(:,1,j))
    end if
  end if
end do
!$OMP END DO
!$OMP END PARALLEL
call omp_free(nthd)
return
end subroutine
!EOC

