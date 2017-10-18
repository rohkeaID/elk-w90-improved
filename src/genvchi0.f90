
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genvchi0(ik,icmp,scsr,vqpl,igq0,gqc,jlgqr,ylmgq,sfacgq,vchi0)
use modmain
implicit none
! local variables
integer, intent(in) :: ik,icmp
real(8), intent(in) :: scsr
real(8), intent(in) :: vqpl(3)
integer, intent(in) :: igq0
real(8), intent(in) :: gqc(ngrf),jlgqr(njcmax,nspecies,ngrf)
complex(8), intent(in) :: ylmgq(lmmaxo,ngrf)
complex(8), intent(in) :: sfacgq(ngrf,natmtot)
complex(8), intent(inout) :: vchi0(nwrf,ngrf,ngrf)
! local variables
integer isym,jk,jkq,iw
integer nst,nstq,ist,jst
integer kst,lst,ig,jg
real(8) vkql(3),eij,t1,t2
complex(8) z1,z2
! automatic arrays
integer idx(nstsv),idxq(nstsv)
! allocatable arrays
complex(8), allocatable :: pmat(:,:,:)
complex(8), allocatable :: wfmt(:,:,:,:),wfir(:,:,:)
complex(8), allocatable :: wfmtq(:,:,:,:),wfirq(:,:,:)
complex(8), allocatable :: zrhomt(:,:),zrhoir(:)
complex(8), allocatable :: zrhoig(:),zw(:)
! k+q-vector in lattice coordinates
vkql(:)=vkl(:,ik)+vqpl(:)
! equivalent reduced k-points for k and k+q
call findkpt(vkl(:,ik),isym,jk)
call findkpt(vkql,isym,jkq)
! count and index states at k and k+q in energy window
nst=0
do ist=1,nstsv
  if (abs(evalsv(ist,jk)-efermi).lt.emaxrf) then
    nst=nst+1
    idx(nst)=ist
  end if
end do
nstq=0
do jst=1,nstsv
  if (abs(evalsv(jst,jkq)-efermi).lt.emaxrf) then
    nstq=nstq+1
    idxq(nstq)=jst
  end if
end do
! generate the wavefunctions for all states at k and k+q in energy window
allocate(wfmt(npcmtmax,natmtot,nspinor,nst),wfir(ngtot,nspinor,nst))
call genwfsvp(.false.,.false.,nst,idx,vkl(:,ik),wfmt,ngtot,wfir)
allocate(wfmtq(npcmtmax,natmtot,nspinor,nstq),wfirq(ngtot,nspinor,nstq))
call genwfsvp(.false.,.false.,nstq,idxq,vkql,wfmtq,ngtot,wfirq)
! read the momentum matrix elements from file
allocate(pmat(nstsv,nstsv,3))
call getpmat(.false.,vkl(:,ik),pmat)
! divide by unit cell volume
t1=1.d0/omega
pmat(:,:,:)=t1*pmat(:,:,:)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zrhomt,zrhoir,zrhoig,zw) &
!$OMP PRIVATE(jst,kst,lst,t1,t2,eij) &
!$OMP PRIVATE(iw,ig,jg,z1,z2)
!$OMP DO
do ist=1,nst
  allocate(zrhomt(npcmtmax,natmtot),zrhoir(ngtot))
  allocate(zrhoig(ngrf),zw(nwrf))
  kst=idx(ist)
  do jst=1,nstq
    lst=idxq(jst)
    t1=wkptnr*omega*(occsv(kst,jk)-occsv(lst,jkq))
    if (abs(t1).lt.1.d-8) cycle
    eij=evalsv(kst,jk)-evalsv(lst,jkq)
! scissor operator
    if (abs(scsr).gt.1.d-8) then
      t2=eij
      if (eij.gt.0.d0) then
        eij=eij+scsr
      else
        eij=eij-scsr
      end if
      t2=eij/t2
! scale the momentum matrix elements
      pmat(kst,lst,:)=t2*pmat(kst,lst,:)
    end if
! frequency-dependent part in response function formula for all frequencies
    do iw=1,nwrf
      zw(iw)=t1/(eij+wrf(iw))
    end do
! compute the complex density in G+q-space
    call genzrho(.true.,.true.,wfmt(:,:,:,ist),wfir(:,:,ist),wfmtq(:,:,:,jst), &
     wfirq(:,:,jst),zrhomt,zrhoir)
    call zftzf(ngrf,jlgqr,ylmgq,ngrf,sfacgq,zrhomt,zrhoir,zrhoig)
!$OMP CRITICAL
!------------------------!
!     body of matrix     !
!------------------------!
    do jg=1,ngrf
      z1=conjg(zrhoig(jg))
      do ig=1,ngrf
        t1=gqc(ig)*gqc(jg)
        if (t1.gt.1.d-8) then
          z2=(fourpi/t1)*zrhoig(ig)*z1
          call zaxpy(nwrf,z2,zw,1,vchi0(:,ig,jg),1)
        end if
      end do
    end do
! special case of q = 0
    if ((gqc(igq0).lt.epslat).and.(abs(eij).gt.1.d-8)) then
!----------------------------------------!
!     head of matrix: G = G' = q = 0     !
!----------------------------------------!
      if (icmp.eq.0) then
! trace of dielectric tensor
        t1=sum(dble(pmat(kst,lst,:))**2+aimag(pmat(kst,lst,:))**2)/3.d0
      else
! particular macroscopic component
        t1=dble(pmat(kst,lst,icmp))**2+aimag(pmat(kst,lst,icmp))**2
      end if
      t1=fourpi*t1/eij**2
      vchi0(:,igq0,igq0)=vchi0(:,igq0,igq0)+t1*zw(:)
!-------------------------!
!     wings of matrix     !
!-------------------------!
      t1=-fourpi/eij
      if (icmp.eq.0) then
        z1=(t1/3.d0)*(pmat(kst,lst,1)+pmat(kst,lst,2)+pmat(kst,lst,3))
      else
        z1=t1*pmat(kst,lst,icmp)
      end if
! G = q = 0
      do ig=2,ngrf
        z2=zrhoig(ig)*conjg(z1)/gqc(ig)
        call zaxpy(nwrf,z2,zw,1,vchi0(:,ig,igq0),1)
      end do
! G' = q = 0
      do jg=2,ngrf
        z2=z1*conjg(zrhoig(jg))/gqc(jg)
        call zaxpy(nwrf,z2,zw,1,vchi0(:,igq0,jg),1)
      end do
    end if
!$OMP END CRITICAL
! end loop over jst
  end do
  deallocate(zrhomt,zrhoir,zrhoig,zw)
! end loop over ist
end do
!$OMP END DO
!$OMP END PARALLEL
deallocate(pmat,wfmt,wfmtq,wfir,wfirq)
return
end subroutine

