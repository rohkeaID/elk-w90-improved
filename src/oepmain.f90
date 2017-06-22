
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine oepmain
use modmain
use modmpi
implicit none
! local variables
integer ik,idm,is,ias,it
integer nr,nri,nrc,nrci,n
real(8) tau,resp,t1
! allocatable arrays
real(8), allocatable :: dvxmt(:,:,:),dvxir(:)
real(8), allocatable :: dbxmt(:,:,:,:),dbxir(:,:)
real(8), allocatable :: rfmt1(:,:,:),rfmt2(:,:),rfir(:)
real(8), allocatable :: rvfmt(:,:,:,:),rvfir(:,:)
complex(8), allocatable :: vclcv(:,:,:,:),vclvv(:,:,:)
! external functions
real(8) rfinp
external rfinp
if (iscl.lt.1) return
! calculate Coulomb matrix elements
allocate(vclcv(ncrmax,natmtot,nstsv,nkpt),vclvv(nstsv,nstsv,nkpt))
call oepvcl(vclcv,vclvv)
! allocate local arrays
allocate(dvxmt(lmmaxvr,nrcmtmax,natmtot),dvxir(ngtot))
allocate(rfmt1(lmmaxvr,nrmtmax,natmtot),rfir(ngtot))
if (spinpol) then
  allocate(dbxmt(lmmaxvr,nrcmtmax,natmtot,ndmag))
  allocate(dbxir(ngtot,ndmag))
  allocate(rvfmt(lmmaxvr,nrmtmax,natmtot,ndmag))
  allocate(rvfir(ngtot,ndmag))
end if
! set the exchange potential and magnetic field to zero
vxmt(:,:,:)=0.d0
vxir(:)=0.d0
if (spinpol) then
  bxmt(:,:,:,:)=0.d0
  bxir(:,:)=0.d0
end if
resp=0.d0
! initial step size
tau=tauoep(1)
!------------------------------!
!     start iteration loop     !
!------------------------------!
do it=1,maxitoep
  if ((mod(it,10).eq.0).and.mp_mpi) then
    write(*,'("Info(oepmain): done ",I4," iterations of ",I4)') it,maxitoep
  end if
! zero the residuals
  dvxmt(:,:,:)=0.d0
  dvxir(:)=0.d0
  if (spinpol) then
    dbxmt(:,:,:,:)=0.d0
    dbxir(:,:)=0.d0
  end if
! calculate the k-dependent residuals
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
  do ik=1,nkpt
! distribute among MPI processes
    if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
    call oepresk(ik,vclcv,vclvv,dvxmt,dvxir,dbxmt,dbxir)
  end do
!$OMP END DO
!$OMP END PARALLEL
! add residuals from each process and redistribute
  if (np_mpi.gt.1) then
    n=lmmaxvr*nrcmtmax*natmtot
    call mpi_allreduce(mpi_in_place,dvxmt,n,mpi_double_precision,mpi_sum, &
     mpi_comm_kpt,ierror)
    call mpi_allreduce(mpi_in_place,dvxir,ngtot,mpi_double_precision, &
     mpi_sum,mpi_comm_kpt,ierror)
    if (spinpol) then
      n=n*ndmag
      call mpi_allreduce(mpi_in_place,dbxmt,n,mpi_double_precision,mpi_sum, &
       mpi_comm_kpt,ierror)
      n=ngtot*ndmag
      call mpi_allreduce(mpi_in_place,dbxir,n,mpi_double_precision,mpi_sum, &
       mpi_comm_kpt,ierror)
    end if
  end if
! convert muffin-tin residuals to spherical harmonics
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(is,nrc,nrci,idm)
!$OMP DO
  do ias=1,natmtot
    is=idxis(ias)
    nrc=nrcmt(is)
    nrci=nrcmtinr(is)
    call rfsht(nrc,nrci,1,dvxmt(:,:,ias),lradstp,rfmt1(:,:,ias))
    do idm=1,ndmag
      call rfsht(nrc,nrci,1,dbxmt(:,:,ias,idm),lradstp,rvfmt(:,:,ias,idm))
    end do
  end do
!$OMP END DO
!$OMP END PARALLEL
! symmetrise the residuals
  call symrf(lradstp,rfmt1,dvxir)
  if (spinpol) call symrvf(lradstp,rvfmt,dbxir)
! magnitude of residuals
  resoep=sqrt(abs(rfinp(lradstp,rfmt1,dvxir,rfmt1,dvxir)))
  do idm=1,ndmag
    t1=rfinp(lradstp,rvfmt(:,:,:,idm),dbxir(:,idm),rvfmt(:,:,:,idm), &
     dbxir(:,idm))
    resoep=resoep+sqrt(abs(t1))
  end do
  resoep=resoep/omega
! adjust step size
  if (it.gt.1) then
    if (resoep.gt.resp) then
      tau=tau*tauoep(2)
    else
      tau=tau*tauoep(3)
    end if
  end if
  resp=resoep
! update exchange potential and magnetic field
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(rfmt2,is,nrc,nrci,idm)
!$OMP DO
  do ias=1,natmtot
    allocate(rfmt2(lmmaxvr,nrcmtmax))
    is=idxis(ias)
    nrc=nrcmt(is)
    nrci=nrcmtinr(is)
! convert residual to spherical coordinates
    call rbsht(nrc,nrci,lradstp,rfmt1(:,:,ias),1,rfmt2)
! subtract from exchange potential
    call rfmtaddc(nrc,nrci,-tau,rfmt2,vxmt(:,:,ias))
! repeat for exchange magnetic field
    do idm=1,ndmag
      call rbsht(nrc,nrci,lradstp,rvfmt(:,:,ias,idm),1,rfmt2)
      call rfmtaddc(nrc,nrci,-tau,rfmt2,bxmt(:,:,ias,idm))
    end do
    deallocate(rfmt2)
  end do
!$OMP END DO
!$OMP END PARALLEL
  vxir(:)=vxir(:)-tau*dvxir(:)
  do idm=1,ndmag
    bxir(:,idm)=bxir(:,idm)-tau*dbxir(:,idm)
  end do
! end iteration loop
end do
! convert the exchange potential and field to spherical harmonics
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(is,nrc,nrci,idm)
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  nrc=nrcmt(is)
  nrci=nrcmtinr(is)
  call rfsht(nrc,nrci,1,vxmt(:,:,ias),lradstp,rfmt1(:,:,ias))
  do idm=1,ndmag
    call rfsht(nrc,nrci,1,bxmt(:,:,ias,idm),lradstp,rvfmt(:,:,ias,idm))
  end do
end do
!$OMP END DO
!$OMP END PARALLEL
! convert potential and field from a coarse to a fine radial mesh
call rfmtctof(rfmt1)
do idm=1,ndmag
  call rfmtctof(rvfmt(:,:,:,idm))
end do
! add to existing (density derived) correlation potential and field
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmtinr(is)
  call rfmtadd(nr,nri,1,rfmt1(:,:,ias),vxcmt(:,:,ias))
  do idm=1,ndmag
    call rfmtadd(nr,nri,1,rvfmt(:,:,ias,idm),bxcmt(:,:,ias,idm))
  end do
end do
vxcir(:)=vxcir(:)+vxir(:)
do idm=1,ndmag
  bxcir(:,idm)=bxcir(:,idm)+bxir(:,idm)
end do
! symmetrise the exchange potential and field
call symrf(1,vxcmt,vxcir)
if (spinpol) call symrvf(1,bxcmt,bxcir)
deallocate(rfmt1,rfir,vclcv,vclvv)
deallocate(dvxmt,dvxir)
if (spinpol) then
  deallocate(rvfmt,rvfir)
  deallocate(dbxmt,dbxir)
end if
return
end subroutine

