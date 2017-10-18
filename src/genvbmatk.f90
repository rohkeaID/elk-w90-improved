
! Copyright (C) 2014 K. Krieger, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genvbmatk(vmt,vir,bmt,bir,ngp,igpig,wfmt,ld,wfir,vbmat)
use modmain
implicit none
! arguments
real(8), intent(in) :: vmt(npcmtmax,natmtot),vir(ngtot)
real(8), intent(in) :: bmt(npcmtmax,natmtot,ndmag),bir(ngtot,ndmag)
integer, intent(in) :: ngp(nspnfv),igpig(ngkmax,nspnfv)
complex(8), intent(in) :: wfmt(npcmtmax,natmtot,nspinor,nstsv)
integer, intent(in) :: ld
complex(8), intent(in) :: wfir(ld,nspinor,nstsv)
complex(8), intent(out) :: vbmat(nstsv,nstsv)
! local variables
integer ist,jst,ispn,jspn
integer is,ias,nrc,nrci,ir
integer igp,ifg
real(8) t0,t1
complex(8) z1,z2,z3
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
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(wfmt1,ias,is,nrc,nrci) &
!$OMP PRIVATE(ist,ispn,z1)
!$OMP DO
do jst=1,nstsv
  allocate(wfmt1(npcmtmax,nspinor))
  do ias=1,natmtot
    is=idxis(ias)
    nrc=nrcmt(is)
    nrci=nrcmti(is)
! apply local potential and magnetic field to spinor wavefunction
    if (ncmag) then
! non-collinear case
      call vbmt1(npcmt(is),vmt(:,ias),bmt(:,ias,1),bmt(:,ias,2),bmt(:,ias,3), &
       wfmt(:,ias,1,jst),wfmt(:,ias,2,jst),wfmt1(:,1),wfmt1(:,2))
    else
! collinear case
      call vbmt2(npcmt(is),vmt(:,ias),bmt(:,ias,1),wfmt(:,ias,1,jst), &
       wfmt(:,ias,2,jst),wfmt1(:,1),wfmt1(:,2))
    end if
    do ist=1,jst
      do ispn=1,nspinor
! compute inner product (functions are in spherical coordinates)
        z1=zfcmtinp(nrc,nrci,rcmt(:,is),r2cmt(:,is),wfmt(:,ias,ispn,ist), &
         wfmt1(:,ispn))
        vbmat(ist,jst)=vbmat(ist,jst)+z1
      end do
    end do
  end do
  deallocate(wfmt1)
end do
!$OMP END DO
!$OMP END PARALLEL
!---------------------------!
!     interstitial part     !
!---------------------------!
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(wfir1,z,ispn,jspn,igp,ifg) &
!$OMP PRIVATE(ir,t0,t1,z1,z2,z3,ist)
!$OMP DO
do jst=1,nstsv
  allocate(wfir1(ngtot,nspinor),z(ngkmax))
! Fourier transform wavefunction to real-space
  do ispn=1,nspinor
    jspn=jspnfv(ispn)
    wfir1(:,ispn)=0.d0
    do igp=1,ngp(jspn)
      ifg=igfft(igpig(igp,jspn))
      wfir1(ifg,ispn)=wfir(igp,ispn,jst)
    end do
    call zfftifc(3,ngridg,1,wfir1(:,ispn))
  end do
! apply local potential and magnetic field to spinor wavefunction
  if (ncmag) then
! non-collinear case
    do ir=1,ngtot
      t0=vir(ir)
      z3=cmplx(bir(ir,1),bir(ir,2),8)
      t1=bir(ir,3)
      z1=wfir1(ir,1)
      z2=wfir1(ir,2)
      wfir1(ir,1)=(t0+t1)*z1+conjg(z3)*z2
      wfir1(ir,2)=(t0-t1)*z2+z3*z1
    end do
  else
! collinear case
    do ir=1,ngtot
      t0=vir(ir)
      t1=bir(ir,1)
      wfir1(ir,1)=(t0+t1)*wfir1(ir,1)
      wfir1(ir,2)=(t0-t1)*wfir1(ir,2)
    end do
  end if
  do ispn=1,nspinor
    jspn=jspnfv(ispn)
! Fourier transform to G+p-space
    call zfftifc(3,ngridg,-1,wfir1(:,ispn))
    do igp=1,ngp(jspn)
      ifg=igfft(igpig(igp,jspn))
      z(igp)=wfir1(ifg,ispn)
    end do
    do ist=1,jst
      z1=zdotc(ngp(jspn),wfir(:,ispn,ist),1,z,1)
      vbmat(ist,jst)=vbmat(ist,jst)+z1
    end do
  end do
  deallocate(wfir1,z)
end do
!$OMP END DO
!$OMP END PARALLEL
! lower triangular part
do ist=1,nstsv
  do jst=1,ist-1
    vbmat(ist,jst)=conjg(vbmat(jst,ist))
  end do
end do
return

contains

subroutine vbmt1(n,v,b1,b2,b3,wf11,wf12,wf21,wf22)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: v(n),b1(n),b2(n),b3(n)
complex(8), intent(in) :: wf11(n),wf12(n)
complex(8), intent(out) :: wf21(n),wf22(n)
! local variables
integer i
real(8) t0,t1
complex(8) z1,z2,z3
do i=1,n
  t0=v(i)
  z3=cmplx(b1(i),b2(i),8)
  t1=b3(i)
  z1=wf11(i)
  z2=wf12(i)
  wf21(i)=(t0+t1)*z1+conjg(z3)*z2
  wf22(i)=(t0-t1)*z2+z3*z1
end do
return
end subroutine

subroutine vbmt2(n,v,b,wf11,wf12,wf21,wf22)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: v(n),b(n)
complex(8), intent(in) :: wf11(n),wf12(n)
complex(8), intent(out) :: wf21(n),wf22(n)
! local variables
integer i
real(8) t0,t1
do i=1,n
  t0=v(i)
  t1=b(i)
  wf21(i)=(t0+t1)*wf11(i)
  wf22(i)=(t0-t1)*wf12(i)
end do
return
end subroutine

end subroutine
