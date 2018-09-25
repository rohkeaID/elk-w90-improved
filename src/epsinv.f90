
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine epsinv
use modmain
use modmpi
use modomp
implicit none
! local variables
integer iq,ik,ig,iw,n
integer info,recl,nthd
! allocatable arrays
integer, allocatable :: ipiv(:)
integer(8), allocatable :: lock(:)
real(8), allocatable :: vgqc(:,:),gqc(:),gclgq(:),jlgqr(:,:,:)
complex(8), allocatable :: ylmgq(:,:),sfacgq(:,:)
complex(8), allocatable :: epsi(:,:,:),work(:)
! allocate local arrays
allocate(vgqc(3,ngrf),gqc(ngrf),gclgq(ngrf))
allocate(jlgqr(njcmax,nspecies,ngrf))
allocate(ylmgq(lmmaxo,ngrf),sfacgq(ngrf,natmtot))
allocate(epsi(ngrf,ngrf,nwrf))
! initialise the OpenMP locks
allocate(lock(nwrf))
do iw=1,nwrf
  call omp_init_lock(lock(iw))
end do
if (mp_mpi) then
! determine the record length for EPSINV.OUT
  inquire(iolength=recl) vql(:,1),ngrf,nwrf,epsi
! open EPSINV.OUT
  open(180,file='EPSINV.OUT',form='UNFORMATTED',access='DIRECT', &
   status='REPLACE',recl=recl)
end if
if (mp_mpi) write(*,*)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
! loop over q-points
do iq=1,nqpt
  if (mp_mpi) write(*,'("Info(epsinv): ",I6," of ",I6," q-points")') iq,nqpt
! generate the G+q-vectors and related quantities
  call gengqrf(vqc(:,iq),vgqc,gqc,jlgqr,ylmgq,sfacgq)
! generate the regularised Coulomb Green's function in G+q-space
  call gengclgq(.true.,iq,ngrf,gqc,gclgq)
! use the symmetric form
  gclgq(:)=sqrt(gclgq(:))
! zero the response function (stored in epsi)
  epsi(:,:,:)=0.d0
  call omp_hold(nkptnr/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
  do ik=1,nkptnr
! distribute among MPI processes
    if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
! compute v^1/2 chi0 v^1/2
    call genvchi0(.false.,ik,lock,0.d0,vql(:,iq),gclgq,jlgqr,ylmgq,sfacgq, &
     ngrf,epsi)
  end do
!$OMP END DO
!$OMP END PARALLEL
  call omp_free(nthd)
! add epsi from each process and redistribute
  if (np_mpi.gt.1) then
    n=nwrf*ngrf*ngrf
    call mpi_allreduce(mpi_in_place,epsi,n,mpi_double_complex,mpi_sum,mpicom, &
     ierror)
  end if
! negate and add delta(G,G')
  epsi(:,:,:)=-epsi(:,:,:)
  do ig=1,ngrf
    epsi(ig,ig,:)=epsi(ig,ig,:)+1.d0
  end do
!-------------------------------------!
!     invert epsilon over G-space     !
!-------------------------------------!
  call omp_hold(nwrf,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ipiv,work,info) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
  do iw=1,nwrf
    allocate(ipiv(ngrf),work(ngrf))
    call zgetrf(ngrf,ngrf,epsi(:,:,iw),ngrf,ipiv,info)
    if (info.eq.0) call zgetri(ngrf,epsi(:,:,iw),ngrf,ipiv,work,ngrf,info)
    if (info.ne.0) then
      write(*,*)
      write(*,'("Error(epsinv): unable to invert epsilon")')
      write(*,'(" for q-point ",I6)') iq
      write(*,'(" and frequency ",I6)') iw
      write(*,*)
      stop
    end if
    deallocate(ipiv,work)
  end do
!$OMP END DO
!$OMP END PARALLEL
  call omp_free(nthd)
! write inverse RPA epsilon to EPSINV.OUT
  if (mp_mpi) write(180,rec=iq) vql(:,iq),ngrf,nwrf,epsi
! end loop over q-points
end do
if (mp_mpi) close(180)
! destroy the OpenMP locks
do iw=1,nwrf
  call omp_destroy_lock(lock(iw))
end do
deallocate(lock)
deallocate(vgqc,gqc,gclgq,jlgqr)
deallocate(ylmgq,sfacgq,epsi)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
return
end subroutine

