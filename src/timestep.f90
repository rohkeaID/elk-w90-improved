
! Copyright (C) 2014 K. Krieger, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine timestep
use modmain
use modtddft
use modmpi
use modomp
implicit none
! local variables
integer ik,is,ias
integer i,j,nthd
real(8) dt,t1
complex(8) z1
! allocatable arrays
real(8), allocatable :: vmt(:,:),vir(:),rfmt(:),w(:)
complex(8), allocatable :: evecsv(:,:),evectv(:,:),evecsvt(:,:)
complex(8), allocatable :: a(:,:),b(:,:),c(:,:)
if (itimes.ge.ntimes) then
  write(*,*)
  write(*,'("Error(timestep): itimes >= ntimes : ",2I8)') itimes,ntimes
  write(*,*)
  stop
end if
! convert muffin-tin Kohn-Sham potential to spherical coordinates
allocate(vmt(npcmtmax,natmtot))
call omp_hold(natmtot,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(rfmt,is) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ias=1,natmtot
  allocate(rfmt(npcmtmax))
  is=idxis(ias)
  call rfmtftoc(nrmt(is),nrmti(is),vsmt(:,ias),rfmt)
  call rbsht(nrcmt(is),nrcmti(is),rfmt,vmt(:,ias))
  deallocate(rfmt)
end do
!$OMP END DO
!$OMP END PARALLEL
call omp_free(nthd)
! multiply interstitial potential by characteristic function
allocate(vir(ngtot))
vir(:)=vsir(:)*cfunir(:)
! loop over k-points
call omp_hold(nkpt/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(w,evecsv,evectv,evecsvt) &
!$OMP PRIVATE(a,b,c,i,j,dt,t1,z1) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  allocate(w(nstsv))
  allocate(evecsv(nstsv,nstsv),evectv(nstsv,nstsv),evecsvt(nstsv,nstsv))
  allocate(a(nstsv,nstsv),b(nstsv,nstsv),c(nstsv,nstsv))
! generate the Hamiltonian matrix in the ground-state second-variational basis
  call genhmlt(ik,vmt,vir,evectv)
! diagonalise the Hamiltonian to get third-variational eigenvectors
  if (spinpol.and.(.not.ncmag)) then
! collinear case requires block diagonalisation
    call eveqnz(nstfv,nstsv,evectv,w)
    i=nstfv+1
    call eveqnz(nstfv,nstsv,evectv(i,i),w(i))
    do i=1,nstfv
      do j=1,nstfv
        evectv(i,j+nstfv)=0.d0
        evectv(i+nstfv,j)=0.d0
      end do
    end do
  else
! non-collinear or spin-unpolarised: full diagonalisation
    call eveqnz(nstsv,nstsv,evectv,w)
  end if
! read in ground-state eigenvectors
  call getevecsv('.OUT',ik,vkl(:,ik),evecsv)
! convert third-variational eigenvectors to first-variational basis
  call zgemm('N','N',nstsv,nstsv,nstsv,zone,evecsv,nstsv,evectv,nstsv,zzero,a, &
   nstsv)
  deallocate(evecsv,evectv)
! time propagate instantaneous eigenvectors across one time step
  dt=times(itimes+1)-times(itimes)
  do i=1,nstsv
    t1=-w(i)*dt
    z1=cmplx(cos(t1),sin(t1),8)
    b(:,i)=z1*a(:,i)
  end do
! read in time-dependent Kohn-Sham eigenvectors
  call getevecsv(filext,ik,vkl(:,ik),evecsvt)
  call zgemm('C','N',nstsv,nstsv,nstsv,zone,a,nstsv,evecsvt,nstsv,zzero,c,nstsv)
  call zgemm('N','N',nstsv,nstsv,nstsv,zone,b,nstsv,c,nstsv,zzero,evecsvt,nstsv)
! write the new eigenvectors to file
  call putevecsv(filext,ik,evecsvt)
  deallocate(w,a,b,c,evecsvt)
end do
!$OMP END DO
!$OMP END PARALLEL
call omp_free(nthd)
deallocate(vmt,vir)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
return
end subroutine

