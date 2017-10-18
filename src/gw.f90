
! Copyright (C) 2016 A. Davydov, A. Sanna, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine gw
use modmain
use modgw
use modmpi
use modtest
implicit none
! local variables
integer ik,iw
real(8) dw,w
! allocatable arrays
real(8), allocatable :: vmt(:,:),vir(:)
real(8), allocatable :: bmt(:,:,:),bir(:,:)
real(8), allocatable :: sft(:),sf(:)
complex(8), allocatable :: swfm(:,:,:)
! initialise universal variables
call init0
call init1
call init2
call init3
! read density and potentials from file
call readstate
! read Fermi energy from file
call readfermi
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! get the eigenvalues and occupancies from file
do ik=1,nkpt
  call getevalsv(filext,ik,vkl(:,ik),evalsv(:,ik))
  call getoccsv(filext,ik,vkl(:,ik),occsv(:,ik))
end do
! write the momentum matrix elements in the second-variational basis to file
call genpmat(.false.,.true.)
! generate the inverse dielectric function and write to file
call epsinv
! compute the matrix elements of -V_xc and -B_xc
allocate(vmt(npcmtmax,natmtot),vir(ngtot))
if (spinpol) then
  allocate(bmt(npcmtmax,natmtot,ndmag),bir(ngtot,ndmag))
end if
call gwlocal(vmt,vir,bmt,bir)
! allocate and zero the total spectral function
allocate(sft(nwplot))
sft(:)=0.d0
! loop over reduced k-point set
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(sf,swfm)
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  allocate(sf(nwplot))
  allocate(swfm(nstsv,nstsv,0:nwfm))
!$OMP CRITICAL
  write(*,'("Info(gw): ",I6," of ",I6," k-points")') ik,nkpt
!$OMP END CRITICAL
! determine the self-energy at the fermionic frequencies
  call selfengyk(ik,vmt,vir,bmt,bir,swfm)
! solve the Dyson equation on the real axis
  call dysonr(ik,swfm,sf)
! write the spectral function to file
!$OMP CRITICAL
  call writesfgw(ik,sf)
!$OMP END CRITICAL
! add to the total spectral function
!$OMP CRITICAL
  sft(:)=sft(:)+wkpt(ik)*sf(:)
!$OMP END CRITICAL
  deallocate(sf,swfm)
end do
!$OMP END DO
!$OMP END PARALLEL
! add total spectral function from each process
if (np_mpi.gt.1) then
  call mpi_allreduce(mpi_in_place,sft,nwplot,mpi_double_precision,mpi_sum, &
   mpi_comm_kpt,ierror)
end if
! write the total spectral function to file (MPI master process only)
if (mp_mpi) then
  open(50,file='TSFGW.OUT',action='WRITE',form='FORMATTED')
  dw=(wplot(2)-wplot(1))/dble(nwplot)
  do iw=1,nwplot
    w=dw*dble(iw-1)+wplot(1)
    write(50,'(2G18.10)') w,sft(iw)
  end do
  write(*,*)
  write(*,'("Info(gw):")')
  write(*,'(" GW spectral functions and Kohn-Sham eigenvalues written to &
   &SFGW_Kkkkkkk.OUT")')
  write(*,'(" for all k-points")')
  write(*,*)
  write(*,'(" Total GW spectral function written to TSFGW.OUT")')
  write(*,*)
  write(*,'(" Fermi energy for the Kohn-Sham eigenvalues is at zero in plots")')
  write(*,'(" Fermi energy for the GW spectral function is undetermined")')
  write(*,*)
  write(*,'(" Spectral function units are states/Hartree/unit cell")')
end if
! write the total GW spectral function to test file
call writetest(600,'total GW spectral function',nv=nwplot,tol=5.d-2,rva=sft)
deallocate(vmt,vir,sft)
if (spinpol) deallocate(bmt,bir)
return
end subroutine

