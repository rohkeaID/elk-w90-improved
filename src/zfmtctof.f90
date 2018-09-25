
! Copyright (C) 2013 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine zfmtctof(zfmt)
use modmain
use modomp
implicit none
! arguments
complex(8), intent(inout) :: zfmt(npmtmax,natmtot)
! local variables
integer is,ias,lm
integer nr,nri,nro
integer iro,ir,npi
integer nrc,nrci,nrco
integer irco,irc,npci
integer i,nthd
! allocatable arrays
real(8), allocatable :: fi1(:),fi2(:),fo1(:),fo2(:)
complex(8), allocatable :: zfmt1(:)
if (lradstp.eq.1) return
call omp_hold(natmtot,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(fi1,fi2,fo1,fo2,zfmt1) &
!$OMP PRIVATE(is,nr,nri,nro,iro,npi) &
!$OMP PRIVATE(nrc,nrci,nrco,irco,npci) &
!$OMP PRIVATE(lm,i,irc,ir) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ias=1,natmtot
  allocate(fi1(nrcmtmax),fi2(nrcmtmax))
  allocate(fo1(nrmtmax),fo2(nrmtmax))
  allocate(zfmt1(npmtmax))
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  nro=nr-nri
  iro=nri+1
  npi=npmti(is)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  nrco=nrc-nrci
  irco=nrci+1
  npci=npcmti(is)
! interpolate up to lmaxi over entire muffin-tin
  do lm=1,lmmaxi
    i=lm
    do irc=1,nrci
      fi1(irc)=dble(zfmt(i,ias))
      fi2(irc)=aimag(zfmt(i,ias))
      i=i+lmmaxi
    end do
    do irc=irco,nrc
      fi1(irc)=dble(zfmt(i,ias))
      fi2(irc)=aimag(zfmt(i,ias))
      i=i+lmmaxo
    end do
    call rfinterp(nrc,rcmt(:,is),fi1,nr,rsp(:,is),fo1)
    call rfinterp(nrc,rcmt(:,is),fi2,nr,rsp(:,is),fo2)
    i=lm
    do ir=1,nri
      zfmt1(i)=cmplx(fo1(ir),fo2(ir),8)
      i=i+lmmaxi
    end do
    do ir=iro,nr
      zfmt1(i)=cmplx(fo1(ir),fo2(ir),8)
      i=i+lmmaxo
    end do
  end do
! interpolate up to lmaxo on outer part of muffin-tin
  do lm=lmmaxi+1,lmmaxo
    i=npci+lm
    do irc=irco,nrc
      fi1(irc)=dble(zfmt(i,ias))
      fi2(irc)=aimag(zfmt(i,ias))
      i=i+lmmaxo
    end do
    call rfinterp(nrco,rcmt(irco,is),fi1(irco),nro,rsp(iro,is),fo1(iro))
    call rfinterp(nrco,rcmt(irco,is),fi2(irco),nro,rsp(iro,is),fo2(iro))
    i=npi+lm
    do ir=iro,nr
      zfmt1(i)=cmplx(fo1(ir),fo2(ir),8)
      i=i+lmmaxo
    end do
  end do
  call zcopy(npmt(is),zfmt1,1,zfmt(:,ias),1)
  deallocate(fi1,fi2,fo1,fo2,zfmt1)
end do
!$OMP END DO
!$OMP END PARALLEL
call omp_free(nthd)
return
end subroutine

