
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine epsinv
use modmain
use modmpi
implicit none
! local variables
integer iq,ik,iw,n
integer igq0,ig,jg
integer info,recl
! allocatable arrays
integer, allocatable :: ipiv(:)
real(8), allocatable :: vgqc(:,:),gqc(:),jlgqr(:,:,:)
complex(8), allocatable :: ylmgq(:,:),sfacgq(:,:)
complex(8), allocatable :: vchi0(:,:,:),epsi(:,:,:),work(:)
! allocate local arrays
allocate(vgqc(3,ngrf),gqc(ngrf),jlgqr(njcmax,nspecies,ngrf))
allocate(ylmgq(lmmaxo,ngrf),sfacgq(ngrf,natmtot))
allocate(vchi0(nwrf,ngrf,ngrf),epsi(ngrf,ngrf,nwrf))
if (mp_mpi) then
! determine the record length for EPSINV.OUT
  inquire(iolength=recl) vql(:,1),ngrf,nwrf,epsi
! open EPSINV.OUT
  open(50,file='EPSINV.OUT',action='WRITE',form='UNFORMATTED',access='DIRECT', &
   status='REPLACE',recl=recl)
end if
! loop over q-points
do iq=1,nqpt
  if (mp_mpi) write(*,'("Info(epsinv): ",I6," of ",I6," q-points")') iq,nqpt
! generate the G+q-vectors and related quantities
  call gengqrf(vqc(:,iq),igq0,vgqc,gqc,jlgqr,ylmgq,sfacgq)
! zero the response function array
  vchi0(:,:,:)=0.d0
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
  do ik=1,nkptnr
! distribute among MPI processes
    if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
! compute v^1/2 chi0 v^1/2
    call genvchi0(ik,0,0.d0,vql(:,iq),igq0,gqc,jlgqr,ylmgq,sfacgq,vchi0)
  end do
!$OMP END DO
!$OMP END PARALLEL
! add vchi0 from each process and redistribute
  if (np_mpi.gt.1) then
    n=nwrf*ngrf*ngrf
    call mpi_allreduce(mpi_in_place,vchi0,n,mpi_double_complex,mpi_sum, &
     mpi_comm_kpt,ierror)
  end if
! negate and add delta(G,G')
  do ig=1,ngrf
    do jg=1,ngrf
      epsi(ig,jg,:)=-vchi0(:,ig,jg)
    end do
    epsi(ig,ig,:)=epsi(ig,ig,:)+1.d0
  end do
!-------------------------------------!
!     invert epsilon over G-space     !
!-------------------------------------!
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ipiv,work,info)
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
! write inverse RPA epsilon to EPSINV.OUT
  if (mp_mpi) write(50,rec=iq) vql(:,iq),ngrf,nwrf,epsi
! end loop over q-points
end do
if (mp_mpi) close(50)
deallocate(vgqc,gqc,jlgqr,ylmgq,sfacgq,vchi0,epsi)
return
end subroutine

