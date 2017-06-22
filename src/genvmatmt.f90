
! Copyright (C) 2007 F. Bultmark, F. Cricchio and L. Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genvmatmt
 !INTERFACE:
subroutine genvmatmt
! !USES:
use modmain
use moddftu
! !DESCRIPTION:
!   Calculate the DFT+U potential matrix to be used in the second-variational
!   step. See {\it Phys. Rev. B} {\bf 52}, 5467 (1995) and {\it Phys. Rev. B}
!   {\bf 80}, 035121 (2009).
!
! !REVISION HISTORY:
!   Created November 2007 (FC,FB,LN,JKD)
!EOP
!BOC
! local variables
integer is,ia,ias,i
integer l,m1,m2,lm1,lm2
real(8) t1,t2
! allocatable arrays
real(8), allocatable :: enfll(:,:)
complex(8), allocatable :: vmfll(:,:,:,:,:)
select case(dftu)
case(0)
  vmatmt(:,:,:,:,:)=0.d0
case(1,2)
! fully localised limit (FLL) or around mean field (AFM)
  call genvmatmt_12
case(3)
! interpolation between the two [PRB 67, 153106 (2003)]
  allocate(enfll(natmmax,ndftu))
  allocate(vmfll(lmmaxdm,nspinor,lmmaxdm,nspinor,natmtot))
  dftu=1
  call genvmatmt_12
  vmfll(:,:,:,:,:)=vmatmt(:,:,:,:,:)
  enfll(:,:)=engyadu(:,:)
  dftu=2
  call genvmatmt_12
! read in the interpolation constant (alpha) if required
  if (readadu) call readalphadu
! reset dftu value
  dftu=3
  do i=1,ndftu
    is=idftu(1,i)
    l=idftu(2,i)
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      t1=alphadu(ia,i)
      t2=1.d0-t1
      lm1=l**2
      do m1=-l,l
        lm1=lm1+1
        lm2=0
        do m2=-l,l
          lm2=lm2+1
          vmatmt(lm1,:,lm2,:,ias)=t1*vmfll(lm1,:,lm2,:,ias) &
           +t2*vmatmt(lm1,:,lm2,:,ias)
        end do
      end do
      engyadu(ia,i)=t1*enfll(ia,i)+t2*engyadu(ia,i)
    end do
  end do
  deallocate(enfll,vmfll)
case default
  write(*,*)
  write(*,'("Error(genvmatmt): invalid dftu : ",I8)') dftu
  write(*,*)
  stop
end select
! add the fixed tensor moment field if required
call ftmfield
return
end subroutine

subroutine genvmatmt_12
use modmain
use moddftu
use modtest
implicit none
! local variables
integer ispn,jspn
integer is,ia,ias,i
integer l,m1,m2,m3,m4,nm
integer lm1,lm2,lm3,lm4
real(8) u,j,n,n0
real(8) mg(3),mg0(3),mg2
real(8) edc,sum1,sum2
complex(8) z1,z2
! automatic arrays
real(8) vee(-lmaxdm:lmaxdm,-lmaxdm:lmaxdm,-lmaxdm:lmaxdm,-lmaxdm:lmaxdm)
complex(8) dm(lmmaxdm,nspinor,lmmaxdm,nspinor)
complex(8) dms(nspinor,nspinor)
! zero the DFT+U potential for each atom
vmatmt(:,:,:,:,:)=0.d0
! zero the DFT+U energy for each atom
engyadu(:,:)=0.d0
! zero the interpolation constants
alphadu(:,:)=0.d0
! loop over DFT+U entries
do i=1,ndftu
  is=idftu(1,i)
  l=idftu(2,i)
  nm=2*l+1
! calculate u, v and the Coulomb matrix elements
  call genveedu(i,u,j,vee)
  if ((abs(u).lt.1.d-10).and.(abs(j).lt.1.d-10)) cycle
! begin loop over atoms
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! copy the density matrix for this atom
    dm(:,:,:,:)=dmatmt(:,:,:,:,ias)
! trace of density matrix for each spin
    dms(:,:)=0.d0
    do ispn=1,nspinor
      do jspn=1,nspinor
        do m1=-l,l
          lm1=idxlm(l,m1)
          dms(ispn,jspn)=dms(ispn,jspn)+dm(lm1,ispn,lm1,jspn)
        end do
      end do
    end do
! trace over spin
    n=dble(dms(1,1))
    if (spinpol) n=n+dble(dms(2,2))
    n0=n/dble(nspinor*nm)
! magnetisation
    if (spinpol) then
      mg(:)=0.d0
      mg(3)=dble(dms(1,1)-dms(2,2))
! non-collinear terms
      if (ncmag) then
        mg(1)=dble(dms(1,2)+dms(2,1))
        mg(2)=dble(zi*(dms(1,2)-dms(2,1)))
      end if
      mg0(:)=mg(:)/dble(nspinor*nm)
    end if
! around mean field (AFM) approach
    if (dftu.eq.2) then
! modify density matrices
      do m1=-l,l
        lm1=idxlm(l,m1)
        if (spinpol) then
          dm(lm1,1,lm1,1)=dm(lm1,1,lm1,1)-(n0+mg0(3))
          dm(lm1,2,lm1,2)=dm(lm1,2,lm1,2)-(n0-mg0(3))
! non-collinear terms
          if (ncmag) then
            dm(lm1,1,lm1,2)=dm(lm1,1,lm1,2)-(mg0(1)-zi*mg0(2))
            dm(lm1,2,lm1,1)=dm(lm1,2,lm1,1)-(mg0(1)+zi*mg0(2))
          end if
        else
! spin-unpolarised case
          dm(lm1,1,lm1,1)=dm(lm1,1,lm1,1)-n0
        end if
      end do
! determine alpha [PRB 67,153106 (2003)]
      sum1=0.d0
      do ispn=1,nspinor
        do m1=-l,l
          lm1=idxlm(l,m1)
          do jspn=1,nspinor
            do m2=-l,l
              lm2=idxlm(l,m2)
              sum1=sum1+dble(dm(lm1,ispn,lm2,jspn)*dm(lm2,jspn,lm1,ispn))
            end do
          end do
        end do
      end do
      if (spinpol) then
        mg2=mg(3)**2
        if (ncmag) mg2=mg2+mg(1)**2+mg(2)**2
      else
        mg2=0.d0
      end if
      sum2=n*(1.d0-0.5d0*n/dble(nm))-0.5d0*mg2/dble(nm)
      if (abs(sum2).gt.1.d-14) then
        alphadu(ia,i)=sum1/sum2
      else
        alphadu(ia,i)=0.d0
      end if
    end if
! calculation of DFT+U potential and energy
! begin loops over m1 and m2
    do m1=-l,l
      lm1=idxlm(l,m1)
      do m2=-l,l
        lm2=idxlm(l,m2)
! begin loops over m3 and m4
        do m3=-l,l
          lm3=idxlm(l,m3)
          do m4=-l,l
            lm4=idxlm(l,m4)
            do ispn=1,nspinor
              do jspn=1,nspinor
                z1=dm(lm2,ispn,lm1,ispn)*dm(lm4,jspn,lm3,jspn)
                z2=dm(lm4,jspn,lm1,ispn)*dm(lm2,ispn,lm3,jspn)
                engyadu(ia,i)=engyadu(ia,i)+dble(z1-z2)*vee(m1,m3,m2,m4)
                vmatmt(lm1,ispn,lm2,ispn,ias)=vmatmt(lm1,ispn,lm2,ispn,ias) &
                 +dm(lm4,jspn,lm3,jspn)*vee(m1,m3,m2,m4)
                vmatmt(lm1,ispn,lm4,jspn,ias)=vmatmt(lm1,ispn,lm4,jspn,ias) &
                 -dm(lm2,ispn,lm3,jspn)*vee(m1,m3,m2,m4)
              end do
            end do
! end loops over m3 and m4
          end do
        end do
! end loops over m1 and m2
      end do
    end do
! multiply energy by factor 1/2
    engyadu(ia,i)=0.5d0*engyadu(ia,i)
! fully localised limit (FLL) approach: double counting corrections
    if (dftu.eq.1) then
      if (spinpol) then
! spin-polarised
        if (ncmag) then
! non-collinear case
! correction to the energy
          edc=0.5d0*u*n*(n-1.d0)
          edc=edc-0.5d0*j*dble(dms(1,1)*(dms(1,1)-1.d0))
          edc=edc-0.5d0*j*dble(dms(2,2)*(dms(2,2)-1.d0))
          edc=edc-0.5d0*j*dble(dms(1,2)*dms(2,1))
          edc=edc-0.5d0*j*dble(dms(2,1)*dms(1,2))
! correction to the potential
          do m1=-l,l
            lm1=idxlm(l,m1)
            vmatmt(lm1,1,lm1,1,ias)=vmatmt(lm1,1,lm1,1,ias) &
             -u*(n-0.5d0)+j*(dms(1,1)-0.5d0)
            vmatmt(lm1,2,lm1,2,ias)=vmatmt(lm1,2,lm1,2,ias) &
             -u*(n-0.5d0)+j*(dms(2,2)-0.5d0)
            vmatmt(lm1,1,lm1,2,ias)=vmatmt(lm1,1,lm1,2,ias)+j*dms(1,2)
            vmatmt(lm1,2,lm1,1,ias)=vmatmt(lm1,2,lm1,1,ias)+j*dms(2,1)
          end do
        else
! collinear case
! correction to the energy
          edc=0.5d0*u*n*(n-1.d0)
          edc=edc-0.5d0*j*dble(dms(1,1)*(dms(1,1)-1.d0))
          edc=edc-0.5d0*j*dble(dms(2,2)*(dms(2,2)-1.d0))
! correction to the potential
          do m1=-l,l
            lm1=idxlm(l,m1)
            vmatmt(lm1,1,lm1,1,ias)=vmatmt(lm1,1,lm1,1,ias) &
             -u*(n-0.5d0)+j*(dms(1,1)-0.5d0)
            vmatmt(lm1,2,lm1,2,ias)=vmatmt(lm1,2,lm1,2,ias) &
             -u*(n-0.5d0)+j*(dms(2,2)-0.5d0)
          end do
        end if
      else
! spin-unpolarised
! correction to the energy
        edc=0.5d0*u*n*(n-1.d0)
        edc=edc-0.5d0*j*n*(n-1.d0)
! correction to the potential
        do m1=-l,l
          lm1=idxlm(l,m1)
          vmatmt(lm1,1,lm1,1,ias)=vmatmt(lm1,1,lm1,1,ias)-u*(n-0.5d0) &
           +j*(n-0.5d0)
        end do
      end if
      engyadu(ia,i)=engyadu(ia,i)-edc
    end if
! trace of dmatmt times vmatmt
    sum1=0.d0
    do ispn=1,nspinor
      do m1=-l,l
        lm1=idxlm(l,m1)
        do jspn=1,nspinor
          do m2=-l,l
            lm2=idxlm(l,m2)
            sum1=sum1+dble(dm(lm1,ispn,lm2,jspn)*vmatmt(lm2,jspn,lm1,ispn,ias))
          end do
        end do
      end do
    end do
! subtract contribution to the energy of DFT+U potential
    engyadu(ia,i)=engyadu(ia,i)-sum1
! end loop over atoms
  end do
! end loop over species
end do
! symmetrise the potential
call symdmat(lmaxdm,lmmaxdm,vmatmt)
! write potential matrix to test file
call writetest(800,'DFT+U energy for each atom',nv=natmmax*ndftu,tol=1.d-4, &
 rva=engyadu)
! write U and J parameters to test file
call writetest(810,'U and J parameters',nv=2*maxdftu,tol=1.d-4,rva=ujdu)
return
end subroutine
!EOC

