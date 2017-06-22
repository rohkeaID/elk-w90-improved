
! Copyright (C) 2014 K. Krieger, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genvbmatk(vmt,vir,bmt,bir,ngp,igpig,wfmt,ld,wfir,vbmat)
use modmain
implicit none
! arguments
real(8), intent(in) :: vmt(lmmaxvr,nrcmtmax,natmtot),vir(ngtot)
real(8), intent(in) :: bmt(lmmaxvr,nrcmtmax,natmtot,ndmag),bir(ngtot,ndmag)
integer, intent(in) :: ngp(nspnfv),igpig(ngkmax,nspnfv)
complex(8), intent(in) :: wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv)
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
complex(8), allocatable :: wfmt1(:,:,:),wfir1(:,:),z(:)
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
  allocate(wfmt1(lmmaxvr,nrcmtmax,nspinor))
  do ias=1,natmtot
    is=idxis(ias)
    nrc=nrcmt(is)
    nrci=nrcmtinr(is)
! apply local potential and magnetic field to spinor wavefunction
    if (ncmag) then
! non-collinear case
      call vbmt1(nrc,nrci,vmt(:,:,ias),bmt(:,:,ias,1),bmt(:,:,ias,2), &
       bmt(:,:,ias,3),wfmt(:,:,ias,1,jst),wfmt(:,:,ias,2,jst),wfmt1(:,:,1), &
       wfmt1(:,:,2))
    else
! collinear case
      call vbmt2(nrc,nrci,vmt(:,:,ias),bmt(:,:,ias,1),wfmt(:,:,ias,1,jst), &
       wfmt(:,:,ias,2,jst),wfmt1(:,:,1),wfmt1(:,:,2))
    end if
    do ist=1,jst
! compute inner product
      do ispn=1,nspinor
        z1=zfcmtinp(nrc,nrci,rcmt(:,is),r2cmt(:,is),wfmt(:,:,ias,ispn,ist), &
         wfmt1(:,:,ispn))
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

subroutine vbmt1(nr,nri,vmt,bmt1,bmt2,bmt3,wfmt11,wfmt12,wfmt21,wfmt22)
implicit none
! arguments
integer, intent(in) :: nr,nri
real(8), intent(in) :: vmt(lmmaxvr,nr)
real(8), intent(in) :: bmt1(lmmaxvr,nr),bmt2(lmmaxvr,nr),bmt3(lmmaxvr,nr)
complex(8), intent(in) :: wfmt11(lmmaxvr,nr),wfmt12(lmmaxvr,nr)
complex(8), intent(out) :: wfmt21(lmmaxvr,nr),wfmt22(lmmaxvr,nr)
! local variables
integer ir,itp
real(8) t0,t1
complex(8) z1,z2,z3
do ir=1,nri
  do itp=1,lmmaxinr
    t0=vmt(itp,ir)
    z3=cmplx(bmt1(itp,ir),bmt2(itp,ir),8)
    t1=bmt3(itp,ir)
    z1=wfmt11(itp,ir)
    z2=wfmt12(itp,ir)
    wfmt21(itp,ir)=(t0+t1)*z1+conjg(z3)*z2
    wfmt22(itp,ir)=(t0-t1)*z2+z3*z1
  end do
end do
do ir=nri+1,nr
  do itp=1,lmmaxvr
    t0=vmt(itp,ir)
    z3=cmplx(bmt1(itp,ir),bmt2(itp,ir),8)
    t1=bmt3(itp,ir)
    z1=wfmt11(itp,ir)
    z2=wfmt12(itp,ir)
    wfmt21(itp,ir)=(t0+t1)*z1+conjg(z3)*z2
    wfmt22(itp,ir)=(t0-t1)*z2+z3*z1
  end do
end do
return
end subroutine

subroutine vbmt2(nr,nri,vmt,bmt,wfmt11,wfmt12,wfmt21,wfmt22)
implicit none
! arguments
integer, intent(in) :: nr,nri
real(8), intent(in) :: vmt(lmmaxvr,nr),bmt(lmmaxvr,nr)
complex(8), intent(in) :: wfmt11(lmmaxvr,nr),wfmt12(lmmaxvr,nr)
complex(8), intent(out) :: wfmt21(lmmaxvr,nr),wfmt22(lmmaxvr,nr)
! local variables
integer ir,itp
real(8) t0,t1
do ir=1,nri
  do itp=1,lmmaxinr
    t0=vmt(itp,ir)
    t1=bmt(itp,ir)
    wfmt21(itp,ir)=(t0+t1)*wfmt11(itp,ir)
    wfmt22(itp,ir)=(t0-t1)*wfmt12(itp,ir)
  end do
end do
do ir=nri+1,nr
  do itp=1,lmmaxvr
    t0=vmt(itp,ir)
    t1=bmt(itp,ir)
    wfmt21(itp,ir)=(t0+t1)*wfmt11(itp,ir)
    wfmt22(itp,ir)=(t0-t1)*wfmt12(itp,ir)
  end do
end do
return
end subroutine

end subroutine
