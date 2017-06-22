
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine elnes
use modmain
use modtest
implicit none
! local variables
integer ik,jk,ikq,ist,jst
integer isym,n,nsk(3),iw
real(8) vkql(3),v(3)
real(8) qc,wd,dw,w,t1
! allocatable arrays
real(8), allocatable :: e(:,:,:),f(:,:,:),ddcs(:)
complex(8), allocatable :: expmt(:,:,:),emat(:,:)
! initialise universal variables
call init0
call init1
! check q-vector is commensurate with k-point grid
v(:)=dble(ngridk(:))*vecql(:)
v(:)=abs(v(:)-nint(v(:)))
if ((v(1).gt.epslat).or.(v(2).gt.epslat).or.(v(3).gt.epslat)) then
  write(*,*)
  write(*,'("Error(elnes): q-vector incommensurate with k-point grid")')
  write(*,'(" ngridk : ",3I6)') ngridk
  write(*,'(" vecql : ",3G18.10)') vecql
  write(*,*)
  stop
end if
! allocate local arrays
allocate(e(nstsv,nstsv,nkptnr))
allocate(f(nstsv,nstsv,nkptnr))
allocate(ddcs(nwplot))
allocate(expmt(lmmaxvr,nrcmtmax,natmtot))
! read in the density and potentials from file
call readstate
! read Fermi energy from file
call readfermi
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! get the second-variational eigenvalues and occupancies from file
do ik=1,nkpt
  call getevalsv(filext,vkl(:,ik),evalsv(:,ik))
  call getoccsv(filext,vkl(:,ik),occsv(:,ik))
end do
! generate the phase factor function exp(iq.r) in the muffin-tins
call genexpmt(vecqc,expmt)
e(:,:,:)=0.d0
f(:,:,:)=0.d0
! begin parallel loop over non-reduced k-points
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(emat,jk,vkql,isym) &
!$OMP PRIVATE(ikq,ist,jst,t1)
!$OMP DO
do ik=1,nkptnr
  allocate(emat(nstsv,nstsv))
!$OMP CRITICAL
  write(*,'("Info(elnes): ",I6," of ",I6," k-points")') ik,nkptnr
!$OMP END CRITICAL
! equivalent reduced k-point
  jk=ivkik(ivk(1,ik),ivk(2,ik),ivk(3,ik))
! k+q-vector in lattice coordinates
  vkql(:)=vkl(:,ik)+vecql(:)
! index to k+q-vector
  call findkpt(vkql,isym,ikq)
! compute < i,k+q | exp(iq.r) | j,k > matrix elements
  call genexpmat(vkl(:,ik),expmt,emat)
! add to the double differential scattering cross-section
  do jst=1,nstsv
    if (evalsv(jst,jk).lt.emaxelnes) then
      do ist=1,nstsv
        e(ist,jst,ik)=evalsv(ist,ikq)-evalsv(jst,jk)
        t1=dble(emat(ist,jst))**2+aimag(emat(ist,jst))**2
        f(ist,jst,ik)=t1*occsv(jst,jk)*(occmax-occsv(ist,ikq))
      end do
    end if
  end do
  deallocate(emat)
end do
!$OMP END DO
!$OMP END PARALLEL
! number of subdivisions used for interpolation
nsk(:)=max(ngrkf/ngridk(:),1)
n=nstsv*nstsv
! integrate over the Brillouin zone
call brzint(nswplot,ngridk,nsk,ivkiknr,nwplot,wplot,n,n,e,f,ddcs)
qc=sqrt(vecqc(1)**2+vecqc(2)**2+vecqc(3)**2)
t1=2.d0/(omega*occmax)
if (qc.gt.epslat) t1=t1/qc**4
ddcs(:)=t1*ddcs(:)
open(50,file='ELNES.OUT',action='WRITE',form='FORMATTED')
wd=wplot(2)-wplot(1)
dw=wd/dble(nwplot)
do iw=1,nwplot
  w=dw*dble(iw-1)+wplot(1)
  write(50,'(2G18.10)') w,ddcs(iw)
end do
close(50)
write(*,*)
write(*,'("Info(elnes):")')
write(*,'(" ELNES double differential cross-section written to ELNES.OUT")')
! write ELNES distribution to test file
call writetest(140,'ELNES cross-section',nv=nwplot,tol=1.d-2,rva=ddcs)
deallocate(e,f,ddcs,expmt)
return
end subroutine

