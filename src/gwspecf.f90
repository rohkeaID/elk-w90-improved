
! Copyright (C) 2016 A. Davydov, A. Sanna, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine gwspecf
use modmain
use modgw
use modmpi
use modomp
use modtest
implicit none
! local variables
integer ik,iw,nthd
real(8) dw,w
! allocatable arrays
real(8), allocatable :: wr(:),sft(:),sf(:)
complex(8), allocatable :: se(:,:,:)
! initialise universal variables
call init0
call init1
call init2
call init3
! read Fermi energy from file
call readfermi
! get the eigenvalues and occupancies from file
do ik=1,nkpt
  call getevalsv(filext,ik,vkl(:,ik),evalsv(:,ik))
  call getoccsv(filext,ik,vkl(:,ik),occsv(:,ik))
end do
! real axis frequencies
allocate(wr(nwplot))
dw=(wplot(2)-wplot(1))/dble(nwplot)
do iw=1,nwplot
  wr(iw)=dw*dble(iw-1)+wplot(1)
end do
! allocate and zero the total spectral function
allocate(sft(nwplot))
sft(:)=0.d0
if (mp_mpi) write(*,*)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
! loop over reduced k-point set
call omp_hold(nkpt/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(sf,se) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  allocate(sf(nwplot))
  allocate(se(nstsv,nstsv,0:nwfm))
!$OMP CRITICAL(gwspecf_1)
  write(*,'("Info(gwspecf): ",I6," of ",I6," k-points")') ik,nkpt
!$OMP END CRITICAL(gwspecf_1)
! get the self-energy at the fermionic frequencies from file
  call getgwsefm(ik,se)
! solve the Dyson equation on the real axis
  call dysonr(ik,wr,se,sf)
! write the spectral function to file
!$OMP CRITICAL(gwspecf_2)
  call writegwsf(ik,sf)
!$OMP END CRITICAL(gwspecf_2)
! add to the total spectral function
!$OMP CRITICAL(gwspecf_3)
  sft(:)=sft(:)+wkpt(ik)*sf(:)
!$OMP END CRITICAL(gwspecf_3)
  deallocate(sf,se)
end do
!$OMP END DO
!$OMP END PARALLEL
call omp_free(nthd)
! add total spectral function from each process
if (np_mpi.gt.1) then
  call mpi_allreduce(mpi_in_place,sft,nwplot,mpi_double_precision,mpi_sum, &
   mpicom,ierror)
end if
! write the total spectral function to file (MPI master process only)
if (mp_mpi) then
  open(50,file='GWTSF.OUT',form='FORMATTED')
  dw=(wplot(2)-wplot(1))/dble(nwplot)
  do iw=1,nwplot
    w=dw*dble(iw-1)+wplot(1)
    write(50,'(2G18.10)') w,sft(iw)
  end do
  close(50)
  write(*,*)
  write(*,'("Info(gw):")')
  write(*,'(" GW spectral functions and Kohn-Sham eigenvalues written to &
   &GWSF_Kkkkkkk.OUT")')
  write(*,'(" for all k-points")')
  write(*,*)
  write(*,'(" Total GW spectral function written to GWTSF.OUT")')
  write(*,*)
  write(*,'(" Fermi energy for the Kohn-Sham eigenvalues is at zero in plots")')
  write(*,'(" Fermi energy for the GW spectral function is undetermined")')
  write(*,*)
  write(*,'(" Spectral function units are states/Hartree/unit cell")')
end if
! write the total GW spectral function to test file
call writetest(610,'total GW spectral function',nv=nwplot,tol=5.d-2,rva=sft)
deallocate(wr,sft)
return
end subroutine

