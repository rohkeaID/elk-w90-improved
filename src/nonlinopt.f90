
! Copyright (C) 2002-2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: nonlinopt
! !INTERFACE:
subroutine nonlinopt
! !USES:
use modmain
! !DESCRIPTION:
!   Calculates susceptibility tensor for non-linear optical second-harmonic
!   generation (SHG). The terms (ztm) are numbered according to Eqs. (49)-(51)
!   of the article {\it Physica Scripta} {\bf T109}, 128 (2004). Other good
!   references are {\it Phys. Rev. B} {\bf 48}, 11705 (1993) and
!   {\it Phys. Rev. B} {\bf 53}, 10751 (1996).
!
! !REVISION HISTORY:
!   Rewrote earlier version, June 2010 (Sharma)
!EOP
!BOC
implicit none
! local variables
integer ik,jk,ist,jst,kst
integer iw,a,b,c,l
! smallest eigenvalue difference allowed in denominator
real(8), parameter :: etol=1.d-4
real(8) eji,eki,ekj,t1
complex(8) pii(3),dji(3),vji(3),vik(3),vkj(3)
complex(8) eta,ztm(3,3),z1
character(256) fname
! allocatable arrays
real(8), allocatable :: w(:)
complex(8), allocatable :: pmat(:,:,:)
complex(8), allocatable :: chiw(:,:),chi2w(:,:)
! initialise universal variables
call init0
call init1
! read Fermi energy from file
call readfermi
! get the eigenvalues and occupancies from file
do ik=1,nkpt
  call getevalsv(filext,vkl(:,ik),evalsv(:,ik))
  call getoccsv(filext,vkl(:,ik),occsv(:,ik))
end do
! generate energy grid (starting from zero)
allocate(w(nwplot))
t1=wplot(2)/dble(nwplot)
do iw=1,nwplot
  w(iw)=t1*dble(iw-1)
end do
! allocate response function arrays
allocate(chiw(nwplot,3))
allocate(chi2w(nwplot,2))
! i divided by the complex relaxation time
eta=cmplx(0.d0,swidth,8)
! begin loop over components
do l=1,noptcomp
  a=optcomp(1,l)
  b=optcomp(2,l)
  c=optcomp(3,l)
  chiw(:,:)=0.d0
  chi2w(:,:)=0.d0
! parallel loop over non-reduced k-points
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(pmat,jk,ist,jst,kst) &
!$OMP PRIVATE(eji,eki,ekj,t1,z1) &
!$OMP PRIVATE(pii,dji,vji,vik,vkj,ztm)
!$OMP DO
  do ik=1,nkptnr
    allocate(pmat(nstsv,nstsv,3))
!$OMP CRITICAL
    write(*,'("Info(nonlinopt): ",I6," of ",I6," k-points")') ik,nkptnr
!$OMP END CRITICAL
! equivalent reduced k-point
    jk=ivkik(ivk(1,ik),ivk(2,ik),ivk(3,ik))
! read momentum matrix elements from file
    call getpmat(.false.,vkl(:,ik),pmat)
! scissor correct the matrix elements
    do ist=1,nstsv
      if (evalsv(ist,jk).lt.efermi) then
        do jst=1,nstsv
          if (evalsv(jst,jk).gt.efermi) then
            eji=evalsv(jst,jk)-evalsv(ist,jk)
            t1=(eji+scissor)/eji
            pmat(ist,jst,:)=t1*pmat(ist,jst,:)
          end if
        end do
      end if
    end do
    z1=(wkptnr*occmax/omega)*zi
! loop over valence states
    do ist=1,nstsv
      if (evalsv(ist,jk).lt.efermi) then
        pii(:)=pmat(ist,ist,:)
! loop over conduction states
        do jst=1,nstsv
          if (evalsv(jst,jk).gt.efermi) then
            eji=evalsv(jst,jk)-evalsv(ist,jk)+scissor
            vji(:)=pmat(jst,ist,:)
            dji(:)=pmat(jst,jst,:)-pii(:)
! loop over intermediate states
            do kst=1,nstsv
              if ((kst.ne.ist).and.(kst.ne.jst)) then
                eki=evalsv(kst,jk)-evalsv(ist,jk)
                ekj=evalsv(kst,jk)-evalsv(jst,jk)
                if (evalsv(kst,jk).gt.efermi) then
                  eki=eki+scissor
                else
                  ekj=ekj-scissor
                end if
                vkj(:)=pmat(kst,jst,:)
                vik(:)=pmat(ist,kst,:)
! interband terms
                t1=-eji*eki*(-ekj)*(eki+ekj)
                if (abs(t1).gt.etol) then
                  t1=1.d0/t1
                else
                  t1=0.d0
                end if
                ztm(1,1)=z1*conjg(vji(a))*(conjg(vkj(b))*conjg(vik(c)) &
                 +conjg(vik(b))*conjg(vkj(c)))*t1
                t1=eji*(-eki)*ekj*(-eki-eji)
                if (abs(t1).gt.etol) then
                  t1=1.d0/t1
                else
                  t1=0.d0
                end if
                ztm(1,2)=0.5d0*z1*vkj(c)*(vji(a)*vik(b)+vik(a)*vji(b))*t1
                t1=eji*(-eki)*ekj*(ekj-eji)
                if (abs(t1).gt.etol) then
                  t1=1.d0/t1
                else
                  t1=0.d0
                end if
                ztm(1,3)=0.5d0*z1*vik(b)*(vkj(c)*vji(a)+vji(c)*vkj(a))*t1
! intraband terms
                t1=(-eki)*ekj*eji**3
                if (abs(t1).gt.etol) then
                  t1=1.d0/t1
                else
                  t1=0.d0
                end if
                ztm(2,1)=0.5d0*z1*(eki*vik(b)*(vkj(c)*vji(a)+vji(c)*vkj(a)) &
                 +ekj*vkj(c)*(vji(a)*vik(b)+vik(a)*vji(b)))*t1
                t1=((-eji)*eki*(-ekj)*eji**2)/(-ekj-eki)
                if (abs(t1).gt.etol) then
                  t1=1.d0/t1
                else
                  t1=0.d0
                end if
                ztm(2,3)=z1*conjg(vji(a))*(conjg(vkj(b))*conjg(vik(c)) &
                 +conjg(vik(b))*conjg(vkj(c)))*t1
! modulation terms
                t1=ekj*(-eki)*eji**3
                if (abs(t1).gt.etol) then
                  t1=1.d0/t1
                else
                  t1=0.d0
                end if
                ztm(3,1)=0.25d0*z1*(-eki)*vkj(a)*(vji(b)*vik(c)+vik(b)*vji(c))*t1
                ztm(3,2)=0.25d0*z1*ekj*vik(a)*(vkj(b)*vji(c)+vji(b)*vkj(c))*t1
!$OMP CRITICAL
                do iw=1,nwplot
! 2w interband
                  chi2w(iw,1)=chi2w(iw,1)+ztm(1,1)/(eji-2.d0*w(iw)+eta)
! 2w intraband
                  chi2w(iw,2)=chi2w(iw,2)+ztm(2,3)/(eji-2.d0*w(iw)+eta)
! w interband
                  chiw(iw,1)=chiw(iw,1)-(ztm(1,2)-ztm(1,3))/(eji-w(iw)+eta)
! w intraband
                  chiw(iw,2)=chiw(iw,2)+ztm(2,1)/(eji-w(iw)+eta)
! w modulation
                  chiw(iw,3)=chiw(iw,3)+(ztm(3,1)-ztm(3,2))/(eji-w(iw)+eta)
                end do
!$OMP END CRITICAL
              end if
! end loop over kst
            end do
            ztm(2,2)=4.d0*z1*conjg(vji(a))*(dji(b)*vji(c)+vji(b)*dji(c))/eji**4
            ztm(3,3)=0.25d0*z1*vji(a)*(vji(b)*dji(c)+dji(b)*vji(c))/eji**4
!$OMP CRITICAL
            do iw=1,nwplot
! 2w intraband
              chi2w(iw,2)=chi2w(iw,2)+ztm(2,2)/(eji-2.d0*w(iw)+eta)
! w modulation
              chiw(iw,3)=chiw(iw,3)+ztm(3,3)/(eji-w(iw)+eta)
            end do
!$OMP END CRITICAL
! end loop over jst
          end if
        end do
! end loop over ist
      end if
    end do
    deallocate(pmat)
! end loop over k-points
  end do
!$OMP END DO
!$OMP END PARALLEL
! write to files
  write(fname,'("CHI_INTER2w_",3I1,".OUT")') a,b,c
  open(51,file=trim(fname),action='WRITE',form='FORMATTED')
  write(fname,'("CHI_INTRA2w_",3I1,".OUT")') a,b,c
  open(52,file=trim(fname),action='WRITE',form='FORMATTED')
  write(fname,'("CHI_INTERw_",3I1,".OUT")') a,b,c
  open(53,file=trim(fname),action='WRITE',form='FORMATTED')
  write(fname,'("CHI_INTRAw_",3I1,".OUT")') a,b,c
  open(54,file=trim(fname),action='WRITE',form='FORMATTED')
  write(fname,'("CHI_",3I1,".OUT")') a,b,c
  open(55,file=trim(fname),action='WRITE',form='FORMATTED')
  do iw=1,nwplot
    write(51,'(2G18.10)') w(iw),dble(chi2w(iw,1))
    write(52,'(2G18.10)') w(iw),dble(chi2w(iw,2))
    write(53,'(2G18.10)') w(iw),dble(chiw(iw,1))
    write(54,'(2G18.10)') w(iw),dble(chiw(iw,2))
    t1=dble(chi2w(iw,1)+chi2w(iw,2)+chiw(iw,1)+chiw(iw,2)+chiw(iw,3))
    write(55,'(2G18.10)') w(iw),t1
  end do
  write(51,'("     ")')
  write(52,'("     ")')
  write(53,'("     ")')
  write(54,'("     ")')
  write(55,'("     ")')
  do iw=1,nwplot
    write(51,'(2G18.10)') w(iw),aimag(chi2w(iw,1))
    write(52,'(2G18.10)') w(iw),aimag(chi2w(iw,2))
    write(53,'(2G18.10)') w(iw),aimag(chiw(iw,1))
    write(54,'(2G18.10)') w(iw),aimag(chiw(iw,2))
    t1=aimag(chi2w(iw,1)+chi2w(iw,2)+chiw(iw,1)+chiw(iw,2)+chiw(iw,3))
    write(55,'(2G18.10)') w(iw),t1
  end do
  close(51); close(52); close(53); close(54); close(55)
! end loop over components
end do
write(*,*)
write(*,'("Info(nonlinopt):")')
write(*,'(" susceptibility tensor written to CHI_abc.OUT")')
write(*,'(" interband contributions written to CHI_INTERx_abc.OUT")')
write(*,'(" intraband contributions written to CHI_INTRAx_abc.OUT")')
write(*,'(" for components")')
do l=1,noptcomp
  write(*,'("  a = ",I1,", b = ",I1,", c = ",I1)') optcomp(1:3,l)
end do
deallocate(w,chiw,chi2w)
return
end subroutine
!EOC

