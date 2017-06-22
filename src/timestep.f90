
! Copyright (C) 2014 K. Krieger, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine timestep
use modmain
use modtddft
use modmpi
implicit none
! local variables
integer ik,is,ias,i,j
integer lwork,info
real(8) dt,t1
complex(8) z1
! allocatable arrays
real(8), allocatable :: vmt(:,:,:),vir(:)
real(8), allocatable :: w(:),rwork(:)
complex(8), allocatable :: evecsv(:,:),evectv(:,:),evecsvt(:,:)
complex(8), allocatable :: a(:,:),b(:,:),c(:,:),work(:)
if (itimes.ge.ntimes) then
  write(*,*)
  write(*,'("Error(timestep): itimes >= ntimes : ",2I8)') itimes,ntimes
  write(*,*)
  stop
end if
! convert muffin-tin Kohn-Sham potential to spherical coordinates
allocate(vmt(lmmaxvr,nrcmtmax,natmtot))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(is)
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  call rbsht(nrcmt(is),nrcmtinr(is),lradstp,vsmt(:,:,ias),1,vmt(:,:,ias))
end do
!$OMP END DO
!$OMP END PARALLEL
! multiply interstitial potential by characteristic function
allocate(vir(ngtot))
vir(:)=vsir(:)*cfunir(:)
lwork=2*nstsv
! loop over k-points
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(w,rwork,work,evecsv,evectv,evecsvt) &
!$OMP PRIVATE(a,b,c,info,i,j,dt,t1,z1)
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  allocate(w(nstsv),rwork(3*nstsv),work(lwork))
  allocate(evecsv(nstsv,nstsv),evectv(nstsv,nstsv),evecsvt(nstsv,nstsv))
  allocate(a(nstsv,nstsv),b(nstsv,nstsv),c(nstsv,nstsv))
! generate the Hamiltonian matrix in the ground-state second-variational basis
  call genhmlt(ik,vmt,vir,evectv)
! diagonalise the Hamiltonian to get third-variational eigenvectors
  if (spinpol.and.(.not.ncmag)) then
! collinear case requires block diagonalisation
    call zheev('V','U',nstfv,evectv,nstsv,w,work,lwork,rwork,info)
    if (info.eq.0) then
      i=nstfv+1
      call zheev('V','U',nstfv,evectv(i,i),nstsv,w(i),work,lwork,rwork,info)
    end if
    do i=1,nstfv
      do j=1,nstfv
        evectv(i,j+nstfv)=0.d0
        evectv(i+nstfv,j)=0.d0
      end do
    end do
  else
! non-collinear or spin-unpolarised: full diagonalisation
    call zheev('V','U',nstsv,evectv,nstsv,w,work,lwork,rwork,info)
  end if
  if (info.ne.0) then
    write(*,*)
    write(*,'("Error(timestep): diagonalisation of the third-variational &
     &Hamiltonian failed")')
    write(*,'(" ZHEEV returned INFO = ",I8)') info
    stop
  end if
! read in ground-state eigenvectors
  call getevecsv('.OUT',vkl(:,ik),evecsv)
! convert third-variational eigenvectors to first-variational basis
  call zgemm('N','N',nstsv,nstsv,nstsv,zone,evecsv,nstsv,evectv,nstsv,zzero,a, &
   nstsv)
  deallocate(evecsv,evectv,rwork,work)
! time propagate instantaneous eigenvectors across one time step
  dt=times(itimes+1)-times(itimes)
  do i=1,nstsv
    t1=-w(i)*dt
    z1=cmplx(cos(t1),sin(t1),8)
    b(:,i)=z1*a(:,i)
  end do
! read in time-dependent Kohn-Sham eigenvectors
  call getevecsv(filext,vkl(:,ik),evecsvt)
  call zgemm('C','N',nstsv,nstsv,nstsv,zone,a,nstsv,evecsvt,nstsv,zzero,c,nstsv)
  call zgemm('N','N',nstsv,nstsv,nstsv,zone,b,nstsv,c,nstsv,zzero,evecsvt,nstsv)
! write the new eigenvectors to file
  call putevecsv(filext,ik,evecsvt)
  deallocate(w,a,b,c,evecsvt)
end do
!$OMP END DO
!$OMP END PARALLEL
deallocate(vmt,vir)
! synchronise MPI processes
call mpi_barrier(mpi_comm_kpt,ierror)
return
end subroutine

