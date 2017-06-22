
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine phonon
use modmain
use modphonon
use modpw
use modmpi
implicit none
! local variables
integer ik,jk,jspn,idm
integer is,ia,ias,ip
integer iv(3),nwork,n
! use Broyden mixing only
integer, parameter :: mtype=3
real(8) ddv,a,b
character(256) fext
! allocatable arrays
real(8), allocatable :: evalfv(:,:),devalsv(:)
real(8), allocatable :: v(:),work(:)
complex(8), allocatable :: dyn(:,:)
complex(8), allocatable :: apwalm(:,:,:,:,:),apwalmq(:,:,:,:,:)
complex(8), allocatable :: dapwalm(:,:,:,:),dapwalmq(:,:,:,:)
complex(8), allocatable :: evecfv(:,:,:),devalfv(:,:),devecfv(:,:,:)
complex(8), allocatable :: evecsv(:,:),devecsv(:,:)
! initialise universal variables
call init0
call init1
call init2
call init4
! check k-point grid is commensurate with q-point grid
iv(:)=mod(ngridk(:),ngridq(:))
if ((iv(1).ne.0).or.(iv(2).ne.0).or.(iv(3).ne.0)) then
  write(*,*)
  write(*,'("Error(phonon): k-point grid incommensurate with q-point grid")')
  write(*,'(" ngridk : ",3I6)') ngridk
  write(*,'(" ngridq : ",3I6)') ngridq
  write(*,*)
  stop
end if
if (spinpol) then
  write(*,*)
  write(*,'("Error(phonon): spin-polarised phonons not yet available")')
  write(*,*)
  stop
end if
! check spin-spiral de-phasing is not used
if (spinsprl.and.ssdph) then
  write(*,*)
  write(*,'("Error(phonon): ssdph should be .false. for DFPT phonons")')
  write(*,*)
  stop
end if
! check for zero atoms
if (natmtot.eq.0) return
! read in the density and potentials
call readstate
! read Fermi energy from file
call readfermi
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! compute the overlap radial integrals
call olprad
! compute the Hamiltonian radial integrals
call hmlrad
! generate the spin-orbit coupling radial functions
call gensocfr
! generate the first- and second-variational eigenvectors and eigenvalues
call genevfsv
! find the occupation numbers
call occupy
! size of mixing vector (complex array)
n=2*(lmmaxvr*nrmtmax*natmtot+ngtot)
if (spinpol) n=n+2*ndmag*(lmmaxvr*nrcmtmax*natmtot+ngtot)
! allocate mixing arrays
if (allocated(v)) deallocate(v)
allocate(v(n))
! determine the size of the mixer work array
nwork=-1
call mixerifc(mtype,n,v,ddv,nwork,v)
allocate(work(nwork))
allocate(dyn(3,natmtot))
! begin new phonon task
10 continue
call dyntask(80,fext)
! if nothing more to do then return
if (iqph.eq.0) return
if (mp_mpi) then
  write(*,'("Info(phonon): working on ",A)') 'DYN'//trim(fext)
! open RMSDDVS.OUT
  open(65,file='RMSDDVS'//trim(fext),action='WRITE',form='FORMATTED')
end if
! zero the dynamical matrix row
dyn(:,:)=0.d0
! check to see if mass is considered infinite
if (spmass(isph).le.0.d0) goto 20
! generate the G+k+q-vectors and related quantities
call gengkqvec(iqph)
! generate the G+q-vectors and related quantities
call gengqvec(iqph)
! generate the characteristic function derivative
call gendcfun
! generate the gradient of the Kohn-Sham potential
call gengvsmt
! initialise the potential derivative
drhomt(:,:,:)=0.d0
drhoir(:)=0.d0
if (spinpol) then
  dmagmt(:,:,:,:)=0.d0
  dmagir(:,:)=0.d0
end if
call dpotks
call gendvsig
! initialise the mixer
iscl=0
call phmixpack(.true.,n,v)
call mixerifc(mtype,n,v,ddv,nwork,work)
! initialise the Fermi energy and occupancy derivatives
defermi=0.d0
doccsv(:,:)=0.d0
! begin the self-consistent loop
do iscl=1,maxscl
! compute the Hamiltonian radial integral derivatives
  call dhmlrad
! zero the density and magnetisation derivatives
  drhomt(:,:,:)=0.d0
  drhoir(:)=0.d0
  if (spinpol) then
    dmagmt(:,:,:,:)=0.d0
    dmagir(:,:)=0.d0
  end if
! parallel loop over k-points
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(evalfv,devalsv,apwalm,apwalmq,dapwalm,dapwalmq) &
!$OMP PRIVATE(evecfv,devalfv,devecfv,evecsv,devecsv) &
!$OMP PRIVATE(jk,jspn)
!$OMP DO
  do ik=1,nkptnr
! distribute among MPI processes
    if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
    allocate(evalfv(nstfv,nspnfv),devalsv(nstsv))
    allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
    allocate(apwalmq(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
    allocate(dapwalm(ngkmax,apwordmax,lmmaxapw,nspnfv))
    allocate(dapwalmq(ngkmax,apwordmax,lmmaxapw,nspnfv))
    allocate(evecfv(nmatmax,nstfv,nspnfv))
    allocate(devalfv(nstfv,nspnfv),devecfv(nmatmax,nstfv,nspnfv))
    allocate(evecsv(nstsv,nstsv),devecsv(nstsv,nstsv))
! equivalent reduced k-point
    jk=ivkik(ivk(1,ik),ivk(2,ik),ivk(3,ik))
! compute the matching coefficients and derivatives
    do jspn=1,nspnfv
      call match(ngk(jspn,ik),gkc(:,jspn,ik),tpgkc(:,:,jspn,ik), &
       sfacgk(:,:,jspn,ik),apwalm(:,:,:,:,jspn))
      call dmatch(iasph,ipph,ngk(jspn,ik),vgkc(:,:,jspn,ik), &
       apwalm(:,:,:,:,jspn),dapwalm(:,:,:,jspn))
      call match(ngkq(jspn,ik),gkqc(:,jspn,ik),tpgkqc(:,:,jspn,ik), &
       sfacgkq(:,:,jspn,ik),apwalmq(:,:,:,:,jspn))
      call dmatch(iasph,ipph,ngkq(jspn,ik),vgkqc(:,:,jspn,ik), &
       apwalmq(:,:,:,:,jspn),dapwalmq(:,:,:,jspn))
    end do
! get the first- and second-variational eigenvalues and eigenvectors from file
    call getevalfv(filext,vkl(:,ik),evalfv)
    call getevecfv(filext,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
    call getevecsv(filext,vkl(:,ik),evecsv)
! solve the first-variational eigenvalue equation derivative
    do jspn=1,nspnfv
      call deveqnfv(ngk(jspn,ik),ngkq(jspn,ik),igkig(:,jspn,ik), &
       igkqig(:,jspn,ik),vgkc(:,:,jspn,ik),vgkqc(:,:,jspn,ik),evalfv(:,jspn), &
       apwalm(:,:,:,:,jspn),apwalmq(:,:,:,:,jspn),dapwalm(:,:,:,jspn), &
       dapwalmq(:,:,:,jspn),evecfv(:,:,jspn),devalfv(:,jspn),devecfv(:,:,jspn))
    end do
    if (spinsprl) then
! solve the spin-spiral second-variational eigenvalue equation derivative
!      call deveqnss(ngk(:,ik),ngkq(:,ik),igkig(:,:,ik),igkqig(:,:,ik),apwalm, &
!       dapwalm,devalfv,evecfv,evecfvq,devecfv,evalsv(:,jk),evecsv,devecsv)
    else
! solve the second-variational eigenvalue equation derivative
!      call deveqnsv(ngk(1,ik),ngkq(1,ik),igkig(:,1,ik),igkqig(:,1,ik), &
!       vgkqc(:,:,1,ik),apwalm,dapwalm,devalfv,evecfv,evecfvq,devecfv, &
!       evalsv(:,jk),evecsv,devecsv)
    end if

!*******
devecsv=0.d0
devalsv(:)=dble(devalfv(:,1))
!*******

! write the eigenvalue/vector derivatives to file
    call putdevalfv(ik,devalfv)
    call putdevalsv(ik,devalsv)
    call putdevecfv(ik,devecfv)
    call putdevecsv(ik,devecsv)
! add to the density and magnetisation derivatives
    call drhomagk(ngk(:,ik),ngkq(:,ik),igkig(:,:,ik),igkqig(:,:,ik), &
     occsv(:,jk),doccsv(:,ik),apwalm,apwalmq,dapwalm,evecfv,devecfv,evecsv, &
     devecsv)
    deallocate(evalfv,devalsv,apwalm,apwalmq,dapwalm,dapwalmq)
    deallocate(evecfv,devalfv,devecfv,evecsv,devecsv)
  end do
!$OMP END DO
!$OMP END PARALLEL
! convert from a coarse to a fine radial mesh
  call zfmtctof(drhomt)
  do idm=1,ndmag
    call zfmtctof(dmagmt(:,:,:,idm))
  end do
! add densities from each process and redistribute
  if (np_mpi.gt.1) then
    n=lmmaxvr*nrmtmax*natmtot
    call mpi_allreduce(mpi_in_place,drhomt,n,mpi_double_complex,mpi_sum, &
     mpi_comm_kpt,ierror)
    call mpi_allreduce(mpi_in_place,drhoir,ngtot,mpi_double_complex,mpi_sum, &
     mpi_comm_kpt,ierror)
    if (spinpol) then
      call mpi_allreduce(mpi_in_place,dmagmt,n,mpi_double_complex,mpi_sum, &
       mpi_comm_kpt,ierror)
      n=ngtot*ndmag
      call mpi_allreduce(mpi_in_place,dmagir,n,mpi_double_complex,mpi_sum, &
       mpi_comm_kpt,ierror)
    end if
  end if
! synchronise MPI processes
  call mpi_barrier(mpi_comm_kpt,ierror)
! add gradient contribution to density derivative
  call gradrhomt
! compute the Kohn-Sham potential derivative
  call dpotks
! pack interstitial and muffin-tin potential and field into one array
  call phmixpack(.true.,n,v)
! mix in the old potential and field with the new
  call mixerifc(mtype,n,v,ddv,nwork,work)
! make sure every MPI process has a numerically identical potential
  if (np_mpi.gt.1) then
    call mpi_bcast(v,n,mpi_double_precision,0,mpi_comm_kpt,ierror)
    call mpi_bcast(ddv,1,mpi_double_precision,0,mpi_comm_kpt,ierror)
  end if
! unpack potential and field
  call phmixpack(.false.,n,v)
  if (mp_mpi) then
    write(65,'(G18.10)') ddv
    call flushifc(65)
  end if
! check for convergence
  if (iscl.ge.2) then
    if (ddv.lt.epspot) goto 20
  end if
! Fourier transform Kohn-Sham potential derivative to G+q-space
  call gendvsig
! compute the occupation number derivatives
  call doccupy
! end the self-consistent loop
end do
write(*,*)
write(*,'("Warning(phonon): failed to reach self-consistency after ",I4,&
 &" loops")') maxscl
20 continue
! close the RMSDDVS.OUT file
if (mp_mpi) close(65)
! generate the dynamical matrix row from force derivatives
call dforce(dyn)
! synchronise MPI processes
call mpi_barrier(mpi_comm_kpt,ierror)
! write dynamical matrix row to file
if (mp_mpi) then
  do ias=1,natmtot
    is=idxis(ias)
    ia=idxia(ias)
    do ip=1,3
      a=dble(dyn(ip,ias))
      b=aimag(dyn(ip,ias))
      if (abs(a).lt.1.d-12) a=0.d0
      if (abs(b).lt.1.d-12) b=0.d0
      write(80,'(2G18.10," : is = ",I4,", ia = ",I4,", ip = ",I4)') a,b,is,ia,ip
    end do
  end do
  close(80)
! write the complex Kohn-Sham potential derivative to file
  call writedvs(fext)
! delete the non-essential files
  call phdelete
end if
! synchronise MPI processes
call mpi_barrier(mpi_comm_kpt,ierror)
goto 10
end subroutine

