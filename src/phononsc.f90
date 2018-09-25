
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine phononsc
use modmain
use modphonon
use modstore
use modmpi
implicit none
! local variables
integer is,ia,ja,ias,jas
integer ip,nph,i,p
real(8) a,b,t1
real(8) ft(3,maxatoms*maxspecies)
complex(8) z1,z2
! allocatable arrays
real(8), allocatable :: vsmt_(:,:),vsir_(:)
complex(8), allocatable :: dyn(:,:)
! store original parameters
natoms_(:)=natoms(:)
avec_(:,:)=avec(:,:)
atposl_(:,:,:)=atposl(:,:,:)
bfcmt0_(:,:,:)=bfcmt0(:,:,:)
mommtfix_(:,:,:)=mommtfix(:,:,:)
tshift_=tshift
tforce_=tforce
autokpt_=autokpt
primcell_=primcell
ngridk_(:)=ngridk(:)
! no shifting of atomic basis allowed
tshift=.false.
! require forces
tforce=.true.
! determine k-point grid size from radkpt
autokpt=.true.
! no reduction to primitive cell
primcell=.false.
! initialise universal variables
call init0
! initialise q-point dependent variables
call init2
! store original parameters
natmtot_=natmtot
idxis_(1:natmtot)=idxis(1:natmtot)
bvec_(:,:)=bvec(:,:)
binv_(:,:)=binv(:,:)
atposc_(:,:,:)=atposc(:,:,:)
ngridg_(:)=ngridg(:)
ngtot_=ngtot
if (allocated(ivg_)) deallocate(ivg_)
allocate(ivg_(3,ngtot_))
ivg_(:,:)=ivg(:,:)
if (allocated(igfft_)) deallocate(igfft_)
allocate(igfft_(ngtot_))
igfft_(:)=igfft(:)
! allocate the Kohn-Sham potential derivative arrays
if (allocated(dvsmt)) deallocate(dvsmt)
allocate(dvsmt(npmtmax,natmtot))
if (allocated(dvsir)) deallocate(dvsir)
allocate(dvsir(ngtot))
! allocate supercell offset vector array
if (allocated(vscph)) deallocate(vscph)
allocate(vscph(3,nqptnr))
! allocate dynamical matrix column
allocate(dyn(3,natmtot))
! begin new phonon task
10 continue
natoms(:)=natoms_(:)
! find a dynamical matrix to calculate
call dyntask(80,filext)
! if nothing more to do then restore original input parameters and return
if (iqph.eq.0) then
  filext='.OUT'
  natoms(:)=natoms_(:)
  avec(:,:)=avec_(:,:)
  atposl(:,:,:)=atposl_(:,:,:)
  bfcmt0(:,:,:)=bfcmt0_(:,:,:)
  mommtfix(:,:,:)=mommtfix_(:,:,:)
  tshift=tshift_
  tforce=tforce_
  autokpt=autokpt_
  primcell=primcell_
  ngridk(:)=ngridk_(:)
  deallocate(ivg_,igfft_)
  return
end if
if (mp_mpi) write(*,'("Info(phononsc): working on ",A)') 'DYN'//trim(filext)
! phonon dry run: just generate empty DYN files
if (task.eq.202) goto 10
! zero the dynamical matrix row
dyn(:,:)=0.d0
! zero the Kohn-Sham potential derivative
dvsmt(:,:)=0.d0
dvsir(:)=0.d0
! check to see if mass is considered infinite
if (spmass(isph).le.0.d0) goto 20
! loop over phases: 0 = cos and 1 = sin displacements
if ((ivq(1,iqph).eq.0).and.(ivq(2,iqph).eq.0).and.(ivq(3,iqph).eq.0)) then
  nph=0
else
  nph=1
end if
! initialise or read the charge density and potentials from file
if (task.eq.200) then
  trdstate=.false.
else
  trdstate=.true.
end if
! loop over cos and sin displacements
do p=0,nph
! generate the supercell with negative displacement
  call genscph(p,-deltaph)
! run the ground-state calculation
  call gndstate
! subsequent calculations will read in this supercell potential
  trdstate=.true.
! store the total force for the first displacement
  do ias=1,natmtot
    ft(:,ias)=forcetot(:,ias)
  end do
! store the Kohn-Sham potential for the first displacement
  allocate(vsmt_(npmtmax,natmtot),vsir_(ngtot))
  vsmt_(:,:)=vsmt(:,:)
  vsir_(:)=vsir(:)
! generate the supercell again with positive displacement
  call genscph(p,deltaph)
! run the ground-state calculation again
  call gndstate
! compute the complex Kohn-Sham potential derivative with implicit q-phase
  call phscdvs(p,vsmt_,vsir_)
  deallocate(vsmt_,vsir_)
! Fourier transform the force differences to obtain the dynamical matrix
  z1=1.d0/(dble(nscph)*2.d0*deltaph)
! multiply by i for sin-like displacement
  if (p.eq.1) z1=z1*zi
  ias=0
  jas=0
  do is=1,nspecies
    ja=0
    do ia=1,natoms_(is)
      ias=ias+1
      do i=1,nscph
        ja=ja+1
        jas=jas+1
        t1=-dot_product(vqc(:,iqph),vscph(:,i))
        z2=z1*cmplx(cos(t1),sin(t1),8)
        do ip=1,3
          t1=-(forcetot(ip,jas)-ft(ip,jas))
          dyn(ip,ias)=dyn(ip,ias)+z2*t1
        end do
      end do
    end do
  end do
end do
20 continue
! write dynamical matrix row to file
if (mp_mpi) then
  ias=0
  do is=1,nspecies
    do ia=1,natoms_(is)
      ias=ias+1
      do ip=1,3
        a=dble(dyn(ip,ias))
        b=aimag(dyn(ip,ias))
        if (abs(a).lt.1.d-12) a=0.d0
        if (abs(b).lt.1.d-12) b=0.d0
        write(80,'(2G18.10," : is = ",I4,", ia = ",I4,", ip = ",I4)') a,b,is, &
         ia,ip
      end do
    end do
  end do
  close(80)
! write the complex Kohn-Sham potential derivative to file
  natoms(:)=natoms_(:)
  natmtot=natmtot_
  idxis(1:natmtot)=idxis_(1:natmtot)
  ngridg(:)=ngridg_(:)
  call writedvs(filext)
! delete the non-essential files
  call phscdelete
end if
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
goto 10
end subroutine

