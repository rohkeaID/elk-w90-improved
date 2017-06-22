
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: init1
! !INTERFACE:
subroutine init1
! !USES:
use modmain
use moddftu
use modtest
use modvars
use modstore
! !DESCRIPTION:
!   Generates the $k$-point set and then allocates and initialises global
!   variables which depend on the $k$-point set.
!
! !REVISION HISTORY:
!   Created January 2004 (JKD)
!EOP
!BOC
implicit none
! local variables
logical lsym(48)
integer is,ia,ias,nppt
integer io,ilo,i1,i2,i3
integer ik,igk,isym,jspn
integer l1,l2,l3,m1,m2,m3
integer lm1,lm2,lm3
real(8) vl(3),vc(3),t1
real(8) ts0,ts1
! external functions
complex(8) gauntyry
external gauntyry

call timesec(ts0)

!---------------------!
!     k-point set     !
!---------------------!
! check if the system is an isolated molecule
if (molecule) then
  ngridk(:)=1
  vkloff(:)=0.d0
  autokpt=.false.
end if
! store the point group symmetries for reducing the k-point set
if (reducek.eq.0) then
  nsymkpt=1
  symkpt(:,:,1)=symlat(:,:,1)
else
  lsym(:)=.false.
  do isym=1,nsymcrys
    if (reducek.eq.2) then
! check symmetry is symmorphic
      t1=abs(vtlsymc(1,isym))+abs(vtlsymc(2,isym))+abs(vtlsymc(3,isym))
      if (t1.gt.epslat) goto 10
! check also that the spin rotation is the same as the spatial rotation
      if (spinpol) then
        if (lspnsymc(isym).ne.lsplsymc(isym)) goto 10
      end if
    end if
    lsym(lsplsymc(isym))=.true.
10 continue
  end do
  nsymkpt=0
  do isym=1,nsymlat
    if (lsym(isym)) then
      nsymkpt=nsymkpt+1
      symkpt(:,:,nsymkpt)=symlat(:,:,isym)
    end if
  end do
end if
if ((task.eq.20).or.(task.eq.21)) then
! generate k-points along a line for band structure plots
  call plotpt1d(bvec,nvp1d,npp1d,vvlp1d,vplp1d,dvp1d,dpp1d)
  nkpt=npp1d
  if (allocated(vkl)) deallocate(vkl)
  allocate(vkl(3,nkpt))
  if (allocated(vkc)) deallocate(vkc)
  allocate(vkc(3,nkpt))
  do ik=1,nkpt
    vkl(:,ik)=vplp1d(:,ik)
    call r3mv(bvec,vkl(:,ik),vkc(:,ik))
  end do
  nkptnr=nkpt
else if (task.eq.25) then
! effective mass calculation
  nkpt=(2*ndspem+1)**3
  if (allocated(ivk)) deallocate(ivk)
  allocate(ivk(3,nkpt))
  if (allocated(vkl)) deallocate(vkl)
  allocate(vkl(3,nkpt))
  if (allocated(vkc)) deallocate(vkc)
  allocate(vkc(3,nkpt))
! map vector to [0,1)
  call r3frac(epslat,vklem)
  ik=0
  do i3=-ndspem,ndspem
    do i2=-ndspem,ndspem
      do i1=-ndspem,ndspem
        ik=ik+1
        ivk(1,ik)=i1; ivk(2,ik)=i2; ivk(3,ik)=i3
        vc(1)=dble(i1); vc(2)=dble(i2); vc(3)=dble(i3)
        vc(:)=vc(:)*deltaem
        call r3mv(binv,vc,vl)
        vkl(:,ik)=vklem(:)+vl(:)
        call r3mv(bvec,vkl(:,ik),vkc(:,ik))
      end do
    end do
  end do
  nkptnr=nkpt
else
! determine the k-point grid automatically from radkpt if required
  if (autokpt) then
    t1=radkpt/twopi
    ngridk(:)=int(t1*sqrt(bvec(1,:)**2+bvec(2,:)**2+bvec(3,:)**2))+1
  end if
! set up the default k-point box
  kptboxl(:,1)=vkloff(:)/dble(ngridk(:))
  if (task.eq.102) kptboxl(:,1)=0.d0
  kptboxl(:,2)=kptboxl(:,1)
  kptboxl(:,3)=kptboxl(:,1)
  kptboxl(:,4)=kptboxl(:,1)
  kptboxl(1,2)=kptboxl(1,2)+1.d0
  kptboxl(2,3)=kptboxl(2,3)+1.d0
  kptboxl(3,4)=kptboxl(3,4)+1.d0
! k-point set and box for Fermi surface plots
  if ((task.eq.100).or.(task.eq.101).or.(task.eq.102)) then
    ngridk(:)=np3d(:)
    if (task.ne.102) kptboxl(:,:)=vclp3d(:,:)
  end if
! allocate the k-point set arrays
  if (allocated(ivkik)) deallocate(ivkik)
  allocate(ivkik(0:ngridk(1)-1,0:ngridk(2)-1,0:ngridk(3)-1))
  if (allocated(ivkiknr)) deallocate(ivkiknr)
  allocate(ivkiknr(0:ngridk(1)-1,0:ngridk(2)-1,0:ngridk(3)-1))
  nkptnr=ngridk(1)*ngridk(2)*ngridk(3)
  if (allocated(ivk)) deallocate(ivk)
  allocate(ivk(3,nkptnr))
  if (allocated(vkl)) deallocate(vkl)
  allocate(vkl(3,nkptnr))
  if (allocated(vkc)) deallocate(vkc)
  allocate(vkc(3,nkptnr))
  if (allocated(wkpt)) deallocate(wkpt)
  allocate(wkpt(nkptnr))
! generate the k-point set
  call genppts(.false.,nsymkpt,symkpt,ngridk,nkptnr,epslat,bvec,kptboxl,nkpt, &
   ivkik,ivkiknr,ivk,vkl,vkc,wkpt,wkptnr)
! write to VARIABLES.OUT
  call writevars('nsymkpt',iv=nsymkpt)
  call writevars('symkpt',nv=9*nsymkpt,iva=symkpt)
  call writevars('ngridk',nv=3,iva=ngridk)
  call writevars('vkloff',nv=3,rva=vkloff)
  call writevars('nkpt',iv=nkpt)
  call writevars('ivkik',nv=nkptnr,iva=ivkik)
  call writevars('ivk',nv=3*nkptnr,iva=ivk)
  call writevars('vkl',nv=3*nkptnr,rva=vkl)
  call writevars('wkpt',nv=nkpt,rva=wkpt)
end if
! write the k-points to test file
call writetest(910,'k-points (Cartesian)',nv=3*nkpt,tol=1.d-8,rva=vkc)

!---------------------!
!     G+k-vectors     !
!---------------------!
if ((xctype(1).lt.0).or.(task.eq.5).or.(task.eq.10).or.(task.eq.205).or. &
 (task.eq.300)) then
  nppt=nkptnr
else
  nppt=nkpt
end if
! find the maximum number of G+k-vectors
call findngkmax(nkpt,vkc,nspnfv,vqcss,ngvec,vgc,gkmax,ngkmax)
! allocate the G+k-vector arrays
if (allocated(ngk)) deallocate(ngk)
allocate(ngk(nspnfv,nppt))
if (allocated(igkig)) deallocate(igkig)
allocate(igkig(ngkmax,nspnfv,nppt))
if (allocated(vgkl)) deallocate(vgkl)
allocate(vgkl(3,ngkmax,nspnfv,nppt))
if (allocated(vgkc)) deallocate(vgkc)
allocate(vgkc(3,ngkmax,nspnfv,nppt))
if (allocated(gkc)) deallocate(gkc)
allocate(gkc(ngkmax,nspnfv,nppt))
if (allocated(tpgkc)) deallocate(tpgkc)
allocate(tpgkc(2,ngkmax,nspnfv,nppt))
if (allocated(sfacgk)) deallocate(sfacgk)
allocate(sfacgk(ngkmax,natmtot,nspnfv,nppt))
do ik=1,nppt
  do jspn=1,nspnfv
    vl(:)=vkl(:,ik)
    vc(:)=vkc(:,ik)
! spin-spiral case
    if (spinsprl) then
      if (jspn.eq.1) then
        vl(:)=vl(:)+0.5d0*vqlss(:)
        vc(:)=vc(:)+0.5d0*vqcss(:)
      else
        vl(:)=vl(:)-0.5d0*vqlss(:)
        vc(:)=vc(:)-0.5d0*vqcss(:)
      end if
    end if
! generate the G+k-vectors
    call gengkvec(ngvec,ivg,vgc,vl,vc,gkmax,ngkmax,ngk(jspn,ik), &
     igkig(:,jspn,ik),vgkl(:,:,jspn,ik),vgkc(:,:,jspn,ik))
! generate the spherical coordinates of the G+k-vectors
    do igk=1,ngk(jspn,ik)
      call sphcrd(vgkc(:,igk,jspn,ik),gkc(igk,jspn,ik),tpgkc(:,igk,jspn,ik))
    end do
! generate structure factors for G+k-vectors
    call gensfacgp(ngk(jspn,ik),vgkc(:,:,jspn,ik),ngkmax,sfacgk(:,:,jspn,ik))
  end do
end do
! write to VARIABLES.OUT
call writevars('nspnfv',iv=nspnfv)
call writevars('ngk',nv=nspnfv*nkpt,iva=ngk)
do ik=1,nkpt
  do jspn=1,nspnfv
    call writevars('igkig',l=jspn,m=ik,nv=ngk(jspn,ik),iva=igkig(:,jspn,ik))
  end do
end do

!---------------------------------!
!     APWs and local-orbitals     !
!---------------------------------!
! allocate linearisation energy arrays
if (allocated(apwe)) deallocate(apwe)
allocate(apwe(maxapword,0:lmaxapw,natmtot))
if (allocated(lorbe)) deallocate(lorbe)
allocate(lorbe(maxlorbord,maxlorb,natmtot))
nlomax=0
lolmax=0
apwordmax=0
do is=1,nspecies
! find the maximum APW order
  do l1=0,lmaxapw
    apwordmax=max(apwordmax,apword(l1,is))
  end do
! set the APW linearisation energies to the default
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do l1=0,lmaxapw
      do io=1,apword(l1,is)
        apwe(io,l1,ias)=apwe0(io,l1,is)
      end do
    end do
  end do
! find the maximum number of local-orbitals
  nlomax=max(nlomax,nlorb(is))
! set the local-orbital linearisation energies to the default
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ilo=1,nlorb(is)
      lolmax=max(lolmax,lorbl(ilo,is))
      do io=1,lorbord(ilo,is)
        lorbe(io,ilo,ias)=lorbe0(io,ilo,is)
      end do
    end do
  end do
end do
lolmmax=(lolmax+1)**2
! generate the local-orbital index
call genidxlo
! allocate radial function arrays
if (allocated(apwfr)) deallocate(apwfr)
allocate(apwfr(nrmtmax,2,apwordmax,0:lmaxapw,natmtot))
if (allocated(apwdfr)) deallocate(apwdfr)
allocate(apwdfr(apwordmax,0:lmaxapw,natmtot))
if (allocated(lofr)) deallocate(lofr)
allocate(lofr(nrmtmax,2,nlomax,natmtot))

!-------------------------!
!     DFT+U variables     !
!-------------------------!
if (dftu.ne.0) then
! allocate energy arrays to calculate Slater integrals with Yukawa potential
  if (allocated(fdue)) deallocate(fdue)
  allocate(fdue(0:lmaxdm,natmtot))
! allocate radial functions to calculate Slater integrals with Yukawa potential
  if (allocated(fdufr)) deallocate(fdufr)
  allocate(fdufr(nrmtmax,0:lmaxdm,natmtot))
end if

!---------------------------------------!
!     eigenvalue equation variables     !
!---------------------------------------!
! total number of empty states (M. Meinert)
nempty=nint(nempty0*max(natmtot,1))
if (nempty.lt.1) nempty=1
! number of first-variational states
nstfv=int(chgval/2.d0)+nempty+1
! overlap and Hamiltonian matrix sizes
if (allocated(nmat)) deallocate(nmat)
allocate(nmat(nspnfv,nkpt))
nmatmax=0
do ik=1,nkpt
  do jspn=1,nspnfv
    nmat(jspn,ik)=ngk(jspn,ik)+nlotot
    if (nstfv.gt.nmat(jspn,ik)) then
      write(*,*)
      write(*,'("Error(init1): number of first-variational states larger than &
       &matrix size")')
      write(*,'("Increase rgkmax or decrease nempty")')
      write(*,*)
      stop
    end if
    nmatmax=max(nmatmax,nmat(jspn,ik))
  end do
end do
! number of second-variational states
nstsv=nstfv*nspinor
! allocate second-variational arrays
if (allocated(evalsv)) deallocate(evalsv)
allocate(evalsv(nstsv,nkpt))
if (allocated(occsv)) deallocate(occsv)
allocate(occsv(nstsv,nkpt))
occsv(:,:)=0.d0
! allocate overlap and Hamiltonian integral arrays
if (allocated(oalo)) deallocate(oalo)
allocate(oalo(apwordmax,nlomax,natmtot))
if (allocated(ololo)) deallocate(ololo)
allocate(ololo(nlomax,nlomax,natmtot))
if (allocated(haa)) deallocate(haa)
allocate(haa(lmmaxvr,apwordmax,0:lmaxmat,apwordmax,0:lmaxmat,natmtot))
if (allocated(hloa)) deallocate(hloa)
allocate(hloa(lmmaxvr,apwordmax,0:lmaxmat,nlomax,natmtot))
if (allocated(hlolo)) deallocate(hlolo)
allocate(hlolo(lmmaxvr,nlomax,nlomax,natmtot))
! allocate and generate complex Gaunt coefficient array
if (allocated(gntyry)) deallocate(gntyry)
allocate(gntyry(lmmaxmat,lmmaxvr,lmmaxmat))
do l1=0,lmaxmat
  do m1=-l1,l1
    lm1=idxlm(l1,m1)
    do l2=0,lmaxvr
      do m2=-l2,l2
        lm2=idxlm(l2,m2)
        do l3=0,lmaxmat
          do m3=-l3,l3
            lm3=idxlm(l3,m3)
            gntyry(lm1,lm2,lm3)=gauntyry(l1,l2,l3,m1,m2,m3)
          end do
        end do
      end do
    end do
  end do
end do
! write to VARIABLES.OUT
call writevars('nempty',iv=nempty)
call writevars('nstfv',iv=nstfv)
call writevars('nlotot',iv=nlotot)
call writevars('nstsv',iv=nstsv)

call timesec(ts1)
timeinit=timeinit+ts1-ts0

return
end subroutine
!EOC

