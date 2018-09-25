
! Copyright (C) 2017 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genzvmatk(zvmt,zvir,nst,ngp,igpig,wfmt,wfir,wfgp,vmat)
use modmain
use modomp
implicit none
! arguments
complex(8), intent(in) :: zvmt(npcmtmax,natmtot),zvir(ngtot)
integer, intent(in) :: nst,ngp,igpig(ngkmax)
complex(8), intent(in) :: wfmt(npcmtmax,natmtot,nspinor,nst)
complex(8), intent(in) :: wfir(ngtot,nspinor,nst)
complex(8), intent(in) :: wfgp(ngkmax,nspinor,nst)
complex(8), intent(out) :: vmat(nst,nst)
!**** nst -> nsvukmax?

! local variables
integer ist,jst,ispn
integer is,ias,nrc,nrci
integer npc,igp,ifg,nthd
! allocatable arrays
complex(8), allocatable :: wfmt1(:),wfir1(:),z(:)
! external functions
complex(8) zfcmtinp,zdotc
external zfcmtinp,zdotc
! zero the matrix elements
vmat(:,:)=0.d0
!-------------------------!
!     muffin-tin part     !
!-------------------------!
call omp_hold(nst,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(wfmt1,ias,is,nrc) &
!$OMP PRIVATE(nrci,npc,ispn,ist) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do jst=1,nst
  allocate(wfmt1(npcmtmax))
  do ias=1,natmtot
    is=idxis(ias)
    nrc=nrcmt(is)
    nrci=nrcmti(is)
    npc=npcmt(is)
    do ispn=1,nspinor
! apply complex potential to wavefunction
      wfmt1(1:npc)=zvmt(1:npc,ias)*wfmt(1:npc,ias,ispn,jst)
      do ist=1,nst
! compute inner product (functions are in spherical coordinates)
        vmat(ist,jst)=vmat(ist,jst)+zfcmtinp(nrc,nrci,rcmt(:,is),r2cmt(:,is), &
         wfmt(:,ias,ispn,ist),wfmt1)
      end do
    end do
  end do
  deallocate(wfmt1)
end do
!$OMP END DO
!$OMP END PARALLEL
call omp_free(nthd)
!---------------------------!
!     interstitial part     !
!---------------------------!
call omp_hold(nst,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(wfir1,z,ispn) &
!$OMP PRIVATE(igp,ifg,ist) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do jst=1,nst
  allocate(wfir1(ngtot),z(ngkmax))
  do ispn=1,nspinor
! apply potential to wavefunction
    wfir1(:)=zvir(:)*wfir(:,ispn,jst)
! Fourier transform to G+p-space
    call zfftifc(3,ngridg,-1,wfir1)
    do igp=1,ngp
      ifg=igfft(igpig(igp))
      z(igp)=wfir1(ifg)
    end do
    do ist=1,nst
! compute inner product
      vmat(ist,jst)=vmat(ist,jst)+zdotc(ngp,wfgp(:,ispn,ist),1,z,1)
    end do
  end do
  deallocate(wfir1,z)
end do
!$OMP END DO
!$OMP END PARALLEL
call omp_free(nthd)
return
end subroutine

