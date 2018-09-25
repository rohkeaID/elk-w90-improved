
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gendwfsv(tsh,tgp,nst,idx,ngp,ngpq,igpqig,apwalmq,dapwalm,evecfv, &
 devecfv,evecsv,devecsv,dwfmt,ld,dwfir)
use modmain
implicit none
! arguments
logical, intent(in) :: tsh,tgp
integer, intent(in) :: nst,idx(nst)
integer, intent(in) :: ngp(nspnfv),ngpq(nspnfv)
integer, intent(in) :: igpqig(ngkmax,nspnfv)
complex(8), intent(in) :: apwalmq(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv)
complex(8), intent(in) :: dapwalm(ngkmax,apwordmax,lmmaxapw,nspnfv)
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv)
complex(8), intent(in) :: devecfv(nmatmax,nstfv,nspnfv)
complex(8), intent(in) :: evecsv(nstsv,nstsv),devecsv(nstsv,nstsv)
complex(8), intent(out) :: dwfmt(npcmtmax,natmtot,nspinor,nst)
integer, intent(in) :: ld
complex(8), intent(out) :: dwfir(ld,nspinor,nst)
! local variables
integer ist,ispn,jspn
integer is,ia,ias,nrc,nrci
integer npc,igp,ifg,i,j,k
real(8) t1
complex(8) z1
! automatic arrays
logical done(nstfv),ddone(nstfv)
! allocatable arrays
complex(8), allocatable :: wfmt1(:,:),wfmt2(:),dwfmt1(:,:)
!---------------------------------------------!
!     muffin-tin wavefunction derivatives     !
!---------------------------------------------!
if (tevecsv) then
  allocate(wfmt1(npcmtmax,nstfv),dwfmt1(npcmtmax,nstfv))
end if
if (.not.tsh) allocate(wfmt2(npcmtmax))
do is=1,nspecies
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  npc=npcmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    done(:)=.false.
    do j=1,nst
      k=idx(j)
      if (tevecsv) then
        i=0
        do ispn=1,nspinor
          jspn=jspnfv(ispn)
          dwfmt(1:npc,ias,ispn,j)=0.d0
          do ist=1,nstfv
            i=i+1
            z1=devecsv(i,k)
!***** check if tq0 is needed here
            if (abs(dble(z1))+abs(aimag(z1)).gt.epsocc) then
              if (.not.done(ist)) then
                if (tsh) then
                  call wavefmt(lradstp,ias,ngp(jspn),apwalmq(:,:,:,ias,jspn), &
                   evecfv(:,ist,jspn),wfmt1(:,ist))
                else
                  call wavefmt(lradstp,ias,ngp(jspn),apwalmq(:,:,:,ias,jspn), &
                   evecfv(:,ist,jspn),wfmt2)
                  call zbsht(nrc,nrci,wfmt2,wfmt1(:,ist))
                end if
                done(ist)=.true.
              end if
              call zaxpy(npc,z1,wfmt1(:,ist),1,dwfmt(:,ias,ispn,j),1)
            end if
            z1=evecsv(i,k)
            if (abs(dble(z1))+abs(aimag(z1)).gt.epsocc) then
              if (.not.ddone(ist)) then
                if (tsh) then
                  call dwavefmt(lradstp,ias,ngp(jspn),ngpq(jspn), &
                   apwalmq(:,:,:,:,jspn),dapwalm(:,:,:,jspn), &
                   evecfv(:,ist,jspn),devecfv(:,ist,jspn),dwfmt1(:,ist))
                else
                  call dwavefmt(lradstp,ias,ngp(jspn),ngpq(jspn), &
                   apwalmq(:,:,:,:,jspn),dapwalm(:,:,:,jspn), &
                   evecfv(:,ist,jspn),devecfv(:,ist,jspn),wfmt2)
                  call zbsht(nrc,nrci,wfmt2,dwfmt1(:,ist))
                end if
                ddone(ist)=.true.
              end if
              call zaxpy(npc,z1,dwfmt1(:,ist),1,dwfmt(:,ias,ispn,j),1)
            end if
          end do
        end do
      else
        if (tsh) then
          call dwavefmt(lradstp,ias,ngp,ngpq,apwalmq,dapwalm,evecfv(:,k,1), &
           devecfv(:,k,1),dwfmt(:,ias,1,j))
        else
          call dwavefmt(lradstp,ias,ngp,ngpq,apwalmq,dapwalm,evecfv(:,k,1), &
           devecfv(:,k,1),wfmt2)
          call zbsht(nrc,nrci,wfmt2,dwfmt(:,ias,1,j))
        end if
      end if
    end do
  end do
end do
if (tevecsv) deallocate(wfmt1,dwfmt1)
if (.not.tsh) deallocate(wfmt2)
!-----------------------------------------------!
!     interstitial wavefunction derivatives     !
!-----------------------------------------------!
t1=1.d0/sqrt(omega)
do j=1,nst
  k=idx(j)
  dwfir(:,:,j)=0.d0
  if (tevecsv) then
    i=0
    do ispn=1,nspinor
      jspn=jspnfv(ispn)
      do ist=1,nstfv
        i=i+1
        z1=devecsv(i,k)
        if (abs(dble(z1))+abs(aimag(z1)).gt.epsocc) then
          if (tgp) then
            do igp=1,ngp(jspn)
              dwfir(igp,ispn,j)=dwfir(igp,ispn,j)+z1*evecfv(igp,ist,jspn)
            end do
          else
            z1=t1*z1
            do igp=1,ngp(jspn)
              ifg=igfft(igpqig(igp,jspn))
              dwfir(ifg,ispn,j)=dwfir(ifg,ispn,j)+z1*evecfv(igp,ist,jspn)
            end do
          end if
        end if
        z1=evecsv(i,k)
        if (abs(dble(z1))+abs(aimag(z1)).gt.epsocc) then
          if (tgp) then
            do igp=1,ngpq(jspn)
              dwfir(igp,ispn,j)=dwfir(igp,ispn,j)+z1*devecfv(igp,ist,jspn)
            end do
          else
            z1=t1*z1
            do igp=1,ngpq(jspn)
              ifg=igfft(igpqig(igp,jspn))
              dwfir(ifg,ispn,j)=dwfir(ifg,ispn,j)+z1*devecfv(igp,ist,jspn)
            end do
          end if
        end if
      end do
    end do
  else
    if (tgp) then
      do igp=1,ngpq(1)
        dwfir(igp,1,j)=devecfv(igp,k,1)
      end do
    else
      do igp=1,ngpq(1)
        ifg=igfft(igpqig(igp,1))
        dwfir(ifg,1,j)=t1*devecfv(igp,k,1)
      end do
    end if
  end if
  if (.not.tgp) then
    do ispn=1,nspinor
      call zfftifc(3,ngridg,1,dwfir(:,ispn,j))
    end do
  end if
end do
return
end subroutine

