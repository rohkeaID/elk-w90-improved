
! Copyright (C) 2007 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gendmat(tspndg,tlmdg,lmin,lmax,ias,ngp,apwalm,evecfv,evecsv,ld,dmat)
use modmain
implicit none
! arguments
logical, intent(in) :: tspndg,tlmdg
integer, intent(in) :: lmin,lmax
integer, intent(in) :: ias
integer, intent(in) :: ngp(nspnfv)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv)
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv)
complex(8), intent(in) :: evecsv(nstsv,nstsv)
integer, intent(in) :: ld
complex(8), intent(out) :: dmat(ld,nspinor,ld,nspinor,nstsv)
! local variables
integer ist,ispn,jspn,is,ia
integer nrc,nrci,nro,iro,irc
integer l,m1,m2,lm1,lm2,i,j
real(8) a,b,t1
complex(8) zq(2),z1
! automatic arrays
logical done(nstfv,nspnfv)
real(8) fr1(nrcmtmax),fr2(nrcmtmax)
! allocatable arrays
complex(8), allocatable :: wfmt1(:,:,:,:),wfmt2(:,:,:)
! external functions
real(8) fintgt
external fintgt
if (lmin.lt.0) then
  write(*,*)
  write(*,'("Error(gendmat): lmin < 0 : ",I8)') lmin
  write(*,*)
  stop
end if
if (lmax.gt.lmaxvr) then
  write(*,*)
  write(*,'("Error(gendmat): lmax > lmaxvr : ",2I8)') lmax,lmaxvr
  write(*,*)
  stop
end if
! allocate local arrays
allocate(wfmt1(lmmaxvr,nrcmtmax,nstfv,nspnfv))
allocate(wfmt2(lmmaxvr,nrcmtmax,nspinor))
! species and atom numbers
is=idxis(ias)
ia=idxia(ias)
nrc=nrcmt(is)
nrci=nrcmtinr(is)
! de-phasing factor for spin-spirals
if (spinsprl.and.ssdph) then
  t1=-0.5d0*dot_product(vqcss(:),atposc(:,ia,is))
  zq(1)=cmplx(cos(t1),sin(t1),8)
  zq(2)=conjg(zq(1))
end if
! zero the density matrix
dmat(:,:,:,:,:)=0.d0
done(:,:)=.false.
! begin loop over second-variational states
do j=1,nstsv
  if (tevecsv) then
! generate spinor wavefunction from second-variational eigenvectors
    wfmt2(:,:,:)=0.d0
    i=0
    do ispn=1,nspinor
      jspn=jspnfv(ispn)
      do ist=1,nstfv
        i=i+1
        z1=evecsv(i,j)
        if (spinsprl.and.ssdph) z1=z1*zq(ispn)
        if (abs(dble(z1))+abs(aimag(z1)).gt.epsocc) then
          if (.not.done(ist,jspn)) then
            call wavefmt(lradstp,ias,ngp(jspn),apwalm(:,:,:,:,jspn), &
             evecfv(:,ist,jspn),wfmt1(:,:,ist,jspn))
            done(ist,jspn)=.true.
          end if
! add to spinor wavefunction
          call zfmtadd(nrc,nrci,z1,wfmt1(:,:,ist,jspn),wfmt2(:,:,ispn))
        end if
      end do
    end do
  else
! spin-unpolarised wavefunction
    call wavefmt(lradstp,ias,ngp,apwalm,evecfv(:,j,1),wfmt2)
  end if
  do ispn=1,nspinor
    do jspn=1,nspinor
      if (tspndg.and.(ispn.ne.jspn)) cycle
      do l=lmin,lmax
        if (l.le.lmaxinr) then
          nro=nrc
          iro=1
        else
          nro=nrc-nrci
          iro=nrci+1
        end if
        do m1=-l,l
          lm1=idxlm(l,m1)
          do m2=-l,l
            lm2=idxlm(l,m2)
            if (tlmdg.and.(lm1.ne.lm2)) cycle
            do irc=iro,nrc
              z1=wfmt2(lm1,irc,ispn)*conjg(wfmt2(lm2,irc,jspn))*r2cmt(irc,is)
              fr1(irc)=dble(z1)
              fr2(irc)=aimag(z1)
            end do
            a=fintgt(-2,nro,rcmt(iro,is),fr1(iro))
            b=fintgt(-2,nro,rcmt(iro,is),fr2(iro))
            dmat(lm1,ispn,lm2,jspn,j)=cmplx(a,b,8)
          end do
        end do
      end do
    end do
  end do
! end loop over second-variational states
end do
deallocate(wfmt1,wfmt2)
return
end subroutine

