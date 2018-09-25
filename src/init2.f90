
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine init2
use modmain
use modrdm
use modvars
implicit none
! local variables
logical lsym(48)
integer isym,iv(3)
integer is,ia,ic
integer ist,m,n
real(8) ts0,ts1
real(8) boxl(3,0:3),t1

call timesec(ts0)

!---------------------!
!     q-point set     !
!---------------------!
! check if the system is an isolated molecule
if (molecule) ngridq(:)=1
! store the point group symmetries for reducing the q-point set
if (reduceq.eq.0) then
  nsymqpt=1
  symqpt(:,:,1)=symlat(:,:,1)
else
  lsym(:)=.false.
  do isym=1,nsymcrys
    lsym(lsplsymc(isym))=.true.
  end do
  nsymqpt=0
  do isym=1,nsymlat
    if (lsym(isym)) then
      nsymqpt=nsymqpt+1
      symqpt(:,:,nsymqpt)=symlat(:,:,isym)
    end if
  end do
end if
if ((task.eq.105).or.(task.eq.180).or.(task.eq.185).or. &
 (task.eq.320).or.(task.eq.330).or.(task.eq.331)) then
! equal k- and q-point grids for nesting function, BSE and linear-reposnse TDDFT
  ngridq(:)=ngridk(:)
else if ((xctype(1).lt.0).or.(task.eq.5).or.(task.eq.300).or.(task.eq.600).or. &
 (task.eq.620)) then
! allow the q-point grid to be smaller than the k-point grid for OEP,
! Hartree-Fock, GW and RDMFT
  if ((ngridq(1).le.0).or.(ngridq(2).le.0).or.(ngridq(3).le.0)) then
    ngridq(:)=ngridk(:)
  end if
else
  ngridq(:)=abs(ngridq(:))
end if
! check that the q-point and k-point grids are commensurate for some tasks
if ((xctype(1).lt.0).or.(task.eq.5).or.(task.eq.205).or.(task.eq.240).or. &
 (task.eq.300).or.(task.eq.600).or.(task.eq.620)) then
  iv(:)=mod(ngridk(:),ngridq(:))
  if ((iv(1).ne.0).or.(iv(2).ne.0).or.(iv(3).ne.0)) then
    write(*,*)
    write(*,'("Error(init2): k-point grid incommensurate with q-point grid")')
    write(*,'(" ngridk : ",3I6)') ngridk
    write(*,'(" ngridq : ",3I6)') ngridq
    write(*,*)
    stop
  end if
end if
! allocate the q-point arrays
if (allocated(iqmap)) deallocate(iqmap)
allocate(iqmap(0:ngridq(1)-1,0:ngridq(2)-1,0:ngridq(3)-1))
if (allocated(iqmapnr)) deallocate(iqmapnr)
allocate(iqmapnr(0:ngridq(1)-1,0:ngridq(2)-1,0:ngridq(3)-1))
nqptnr=ngridq(1)*ngridq(2)*ngridq(3)
if (allocated(ivq)) deallocate(ivq)
allocate(ivq(3,nqptnr))
if (allocated(vql)) deallocate(vql)
allocate(vql(3,nqptnr))
if (allocated(vqc)) deallocate(vqc)
allocate(vqc(3,nqptnr))
if (allocated(wqpt)) deallocate(wqpt)
allocate(wqpt(nqptnr))
! set up the q-point box (offset should always be zero)
boxl(:,:)=0.d0
boxl(1,1)=1.d0; boxl(2,2)=1.d0; boxl(3,3)=1.d0
! generate the q-point set
! (note that the vectors vql and vqc are in the first Brillouin zone)
call genppts(.true.,nsymqpt,symqpt,ngridq,nqptnr,epslat,bvec,boxl,nqpt,iqmap, &
 iqmapnr,ivq,vql,vqc,wqpt,wqptnr)
! write the q-points to QPOINTS.OUT
call writeqpts
! find the index for q = 0
do iq0=1,nqpt
  t1=sum(abs(vql(:,iq0)))
  if (t1.lt.epslat) goto 10
end do
write(*,*)
write(*,'("Error(init2): q = 0 not in q-point set")')
write(*,*)
stop
10 continue
! find the maximum number of spherical Bessel functions on the coarse radial
! mesh over all species
njcmax=1
do is=1,nspecies
  n=(lmaxi+1)*nrcmti(is)+(lmaxo+1)*(nrcmt(is)-nrcmti(is))
  if (n.gt.njcmax) njcmax=n
end do
! write to VARIABLES.OUT
call writevars('nsymqpt',iv=nsymqpt)
call writevars('symqpt',nv=9*nsymqpt,iva=symqpt)
call writevars('ngridq',nv=3,iva=ngridq)
call writevars('nqpt',iv=nqpt)
call writevars('iqmap',nv=nqptnr,iva=iqmap)
call writevars('ivq',nv=3*nqptnr,iva=ivq)
call writevars('vql',nv=3*nqptnr,rva=vql)
call writevars('wqpt',nv=nqpt,rva=wqpt)

!--------------------------------------------------------!
!     OEP, Hartree-Fock, RDMFT, BSE and GW variables     !
!--------------------------------------------------------!
if ((xctype(1).lt.0).or.(task.eq.5).or.(task.eq.180).or.(task.eq.185).or. &
 (task.eq.188).or.(task.eq.205).or.(task.eq.300).or.(task.eq.320).or. &
 (task.eq.330).or.(task.eq.331).or.(task.eq.600).or.(task.eq.620)) then
! determine the regularised Coulomb Green's function for small q
  call gengclq
! output the Coulomb Green's function to GCLQ.OUT
  call writegclq
end if
if (xctype(1).lt.0) then
! initialise OEP residual magnitude
  resoep=1.d0
! find maximum core states over all species
  ncrmax=0
  do is=1,nspecies
    do ia=1,natoms(is)
      ic=0
      do ist=1,nstsp(is)
        if (spcore(ist,is)) then
          do m=-ksp(ist,is),ksp(ist,is)-1
            ic=ic+1
          end do
        end if
      end do
      ncrmax=max(ncrmax,ic)
    end do
  end do
! allocate and zero the exchange potential and magnetic field
  if (allocated(vxmt)) deallocate(vxmt)
  allocate(vxmt(npcmtmax,natmtot))
  vxmt(:,:)=0.d0
  if (allocated(vxir)) deallocate(vxir)
  allocate(vxir(ngtot))
  vxir(:)=0.d0
  if (spinpol) then
    if (allocated(bxmt)) deallocate(bxmt)
    allocate(bxmt(npcmtmax,natmtot,ndmag))
    bxmt(:,:,:)=0.d0
    if (allocated(bxir)) deallocate(bxir)
    allocate(bxir(ngtot,ndmag))
    bxir(:,:)=0.d0
  end if
end if
if (task.eq.300) then
  if (allocated(vclmat)) deallocate(vclmat)
  allocate(vclmat(nstsv,nstsv,nkpt))
  if (allocated(dkdc)) deallocate(dkdc)
  allocate(dkdc(nstsv,nstsv,nkpt))
end if

!----------------------------------------------------------!
!     G-vector variables for coarse grid (G < 2*gkmax)     !
!----------------------------------------------------------!
! generate the G-vectors
call gengvc
! generate the characteristic function
call gencfrc

call timesec(ts1)
timeinit=timeinit+ts1-ts0

return
end subroutine

