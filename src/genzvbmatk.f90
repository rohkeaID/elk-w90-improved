
! Copyright (C) 2018 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genzvbmatk(zvmt,zvir,zbmt,zbir,nst,ngp,igpig,wfmt,wfir,wfgp,vbmat)
use modmain
use modomp
implicit none
! arguments
complex(8), intent(in) :: zvmt(npcmtmax,natmtot),zvir(ngtot)
complex(8), intent(in) :: zbmt(npcmtmax,natmtot,ndmag),zbir(ngtot,ndmag)
integer, intent(in) :: nst,ngp,igpig(ngkmax)
complex(8), intent(in) :: wfmt(npcmtmax,natmtot,nspinor,nst)
complex(8), intent(in) :: wfir(ngtot,nspinor,nst)
complex(8), intent(in) :: wfgp(ngkmax,nspinor,nst)
complex(8), intent(out) :: vbmat(nst,nst)
! local variables
integer ist,jst,ispn
integer is,ias,nrc,nrci,npc
integer igp,ifg,nthd
! allocatable arrays
complex(8), allocatable :: wfmt1(:,:),wfir1(:,:),z(:)
! external functions
complex(8) zfcmtinp,zdotc
external zfcmtinp,zdotc
! zero the matrix elements
vbmat(:,:)=0.d0
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
  allocate(wfmt1(npcmtmax,nspinor))
  do ias=1,natmtot
    is=idxis(ias)
    nrc=nrcmt(is)
    nrci=nrcmti(is)
    npc=npcmt(is)
    do ispn=1,nspinor
      call zcopy(npc,wfmt(:,ias,ispn,jst),1,wfmt1(:,ispn),1)
    end do
! apply local potential and magnetic field to spinor wavefunction
    if (ncmag) then
! non-collinear case
      call zvbmk1(npc,zvmt(:,ias),zbmt(:,ias,1),zbmt(:,ias,2),zbmt(:,ias,3), &
       wfmt1(:,1),wfmt1(:,2))
    else
! collinear case
      call zvbmk2(npc,zvmt(:,ias),zbmt(:,ias,1),wfmt1(:,1),wfmt1(:,2))
    end if
    do ist=1,nst
      do ispn=1,nspinor
! compute inner product (functions are in spherical coordinates)
        vbmat(ist,jst)=vbmat(ist,jst)+zfcmtinp(nrc,nrci,rcmt(:,is), &
         r2cmt(:,is),wfmt(:,ias,ispn,ist),wfmt1(:,ispn))
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
  allocate(wfir1(ngtot,nspinor),z(ngkmax))
  do ispn=1,nspinor
    call zcopy(ngtot,wfir(:,ispn,jst),1,wfir1(:,ispn),1)
  end do
! apply local potential and magnetic field to spinor wavefunction
  if (ncmag) then
! non-collinear case
    call zvbmk1(ngtot,zvir,zbir(:,1),zbir(:,2),zbir(:,3),wfir1(:,1),wfir1(:,2))
  else
! collinear case
    call zvbmk2(ngtot,zvir,zbir,wfir1(:,1),wfir1(:,2))
  end if
  do ispn=1,nspinor
! Fourier transform to G+p-space
    call zfftifc(3,ngridg,-1,wfir1(:,ispn))
    do igp=1,ngp
      ifg=igfft(igpig(igp))
      z(igp)=wfir1(ifg,ispn)
    end do
    do ist=1,nst
      vbmat(ist,jst)=vbmat(ist,jst)+zdotc(ngp,wfgp(:,ispn,ist),1,z,1)
    end do
  end do
  deallocate(wfir1,z)
end do
!$OMP END DO
!$OMP END PARALLEL
call omp_free(nthd)
return

contains

subroutine zvbmk1(n,zv,zb1,zb2,zb3,wf1,wf2)
implicit none
! arguments
integer, intent(in) :: n
complex(8), intent(in) :: zv(n),zb1(n),zb2(n),zb3(n)
complex(8), intent(inout) :: wf1(n),wf2(n)
! local variables
integer i
complex(8) v,b1,b2,b3,z1,z2
do i=1,n
  v=zv(i)
  b1=zb1(i); b2=zb2(i); b3=zb3(i)
  b2=cmplx(-aimag(b2),dble(b2),8)
  z1=wf1(i); z2=wf2(i)
  wf1(i)=(v+b3)*z1+(b1-b2)*z2
  wf2(i)=(v-b3)*z2+(b1+b2)*z1
end do
return
end subroutine

subroutine zvbmk2(n,zv,zb,wf1,wf2)
implicit none
! arguments
integer, intent(in) :: n
complex(8), intent(in) :: zv(n),zb(n)
complex(8), intent(inout) :: wf1(n),wf2(n)
! local variables
integer i
complex(8) v,b3
do i=1,n
  v=zv(i)
  b3=zb(i)
  wf1(i)=(v+b3)*wf1(i)
  wf2(i)=(v-b3)*wf2(i)
end do
return
end subroutine

end subroutine

