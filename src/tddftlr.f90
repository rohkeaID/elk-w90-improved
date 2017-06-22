
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine tddftlr
use modmain
use modtddft
use modmpi
implicit none
! local variables
integer, parameter :: maxit=500
integer ik,iw,n
integer igq0,ig,jg
integer it,info
real(8) v(3),t1
complex(8) vfxcp,z1
! allocatable arrays
integer, allocatable :: ipiv(:)
real(8), allocatable :: vgqc(:,:),gqc(:)
complex(8), allocatable :: ylmgq(:,:),sfacgq(:,:)
complex(8), allocatable :: vchi0(:,:,:),vfxc(:,:,:)
complex(8), allocatable :: eps0(:,:,:),eps(:,:,:)
complex(8), allocatable :: a(:,:),b(:,:),work(:)
! initialise global variables
call init0
call init1
call init2
call init3
! check q-vector is commensurate with k-point grid
v(:)=dble(ngridk(:))*vecql(:)
v(:)=abs(v(:)-nint(v(:)))
if ((v(1).gt.epslat).or.(v(2).gt.epslat).or.(v(3).gt.epslat)) then
  write(*,*)
  write(*,'("Error(tddftlr): q-vector incommensurate with k-point grid")')
  write(*,'(" ngridk : ",3I6)') ngridk
  write(*,'(" vecql : ",3G18.10)') vecql
  write(*,*)
  stop
end if
! read density and potentials from file
call readstate
! read Fermi energy from a file
call readfermi
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! get the eigenvalues and occupancies from file
do ik=1,nkpt
  call getevalsv(filext,vkl(:,ik),evalsv(:,ik))
  call getoccsv(filext,vkl(:,ik),occsv(:,ik))
end do
! generate the G+q-vectors and related quantities
allocate(vgqc(3,ngrf),gqc(ngrf))
allocate(ylmgq(lmmaxvr,ngrf),sfacgq(ngrf,natmtot))
call gengqrf(vecqc,igq0,vgqc,gqc,ylmgq,sfacgq)
! allocate local arrays
allocate(vchi0(nwrf,ngrf,ngrf),vfxc(ngrf,ngrf,nwrf))
allocate(eps0(ngrf,ngrf,nwrf),eps(ngrf,ngrf,nwrf))
! compute v^1/2 chi0 v^1/2 (the symmetric version of v chi0)
vchi0(:,:,:)=0.d0
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
do ik=1,nkptnr
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
!$OMP CRITICAL
  write(*,'("Info(tddftlr): ",I6," of ",I6," k-points")') ik,nkptnr
!$OMP END CRITICAL
! compute v^1/2 chi0 v^1/2
  call genvchi0(ik,optcomp(1,1),scissor,vecql,igq0,gqc,ylmgq,sfacgq,vchi0)
end do
!$OMP END DO
!$OMP END PARALLEL
! add vchi0 from each process and redistribute
if (np_mpi.gt.1) then
  n=nwrf*ngrf*ngrf
  call mpi_allreduce(mpi_in_place,vchi0,n,mpi_double_complex,mpi_sum, &
   mpi_comm_kpt,ierror)
end if
! calculate symmetric epsilon = 1 - v^1/2 chi0 v^1/2
do ig=1,ngrf
  do jg=1,ngrf
    eps0(ig,jg,:)=-vchi0(:,ig,jg)
    eps(ig,jg,:)=vchi0(:,ig,jg)
  end do
  eps0(ig,ig,:)=eps0(ig,ig,:)+1.d0
  eps(ig,ig,:)=eps(ig,ig,:)+1.d0
end do
vfxcp=0.d0
it=0
10 continue
! compute vchi0 v^(-1/2) f_xc v^(-1/2) vchi0
call genvfxc(gqc,vchi0,eps0,eps,vfxc)
allocate(ipiv(ngrf),a(ngrf,ngrf),b(ngrf,ngrf),work(ngrf))
! begin loop over frequencies
do iw=1,nwrf
! compute 1 - v^1/2 chi0 v^1/2 - v^(-1/2) f_xc v^(-1/2) vchi0
  a(:,:)=eps0(:,:,iw)-vfxc(:,:,iw)
! invert this matrix
  call zgetrf(ngrf,ngrf,a,ngrf,ipiv,info)
  if (info.eq.0) call zgetri(ngrf,a,ngrf,ipiv,work,ngrf,info)
  if (info.ne.0) then
    write(*,*)
    write(*,'("Error(tddftlr): unable to invert epsilon")')
    write(*,'(" for frequency ",I6)') iw
    write(*,*)
    stop
  end if
! left multiply by v^1/2 chi0 v^1/2
  b(:,:)=vchi0(iw,:,:)
  call zgemm('N','N',ngrf,ngrf,ngrf,zone,b,ngrf,a,ngrf,zzero,eps(:,:,iw),ngrf)
! compute epsilon^(-1) = 1 + v^1/2 chi v^1/2
  do ig=1,ngrf
    eps(ig,ig,iw)=1.d0+eps(ig,ig,iw)
  end do
end do
deallocate(ipiv,a,b,work)
! bootstrap f_xc
if (fxctype(1).eq.210) then
  it=it+1
  if (it.gt.maxit) then
    write(*,*)
    write(*,'("Error(tddftlr): bootstrap kernel failed to converge")')
    write(*,*)
    stop
  end if
  if (mod(it,10).eq.0) then
    write(*,'("Info(tddftlr): done ",I4," bootstrap iterations")') it
    write(*,'(" head of matrix v.f_xc : ",2G18.10)') vfxc(1,1,1)
  end if
! check for convergence
  t1=abs(vfxcp)-abs(vfxc(1,1,1))
  vfxcp=vfxc(1,1,1)
  if (abs(t1).gt.1.d-8) goto 10
end if
! write G = G' = 0 components to file
if (mp_mpi) then
  open(50,file="EPSILON_TDDFT.OUT",action='WRITE',form='FORMATTED')
  open(51,file="EELS_TDDFT.OUT",action='WRITE',form='FORMATTED')
  do iw=1,nwrf
    z1=1.d0/eps(1,1,iw)
    write(50,'(2G18.10)') dble(wrf(iw)),dble(z1)
    write(51,'(2G18.10)') dble(wrf(iw)),-dble(eps(1,1,iw))
  end do
  write(50,*)
  write(51,*)
  do iw=1,nwrf
    z1=1.d0/eps(1,1,iw)
    write(50,'(2G18.10)') dble(wrf(iw)),aimag(z1)
    write(51,'(2G18.10)') dble(wrf(iw)),-aimag(eps(1,1,iw))
  end do
  close(50)
  close(51)
  write(*,*)
  write(*,'("Info(tddftlr):")')
  write(*,'(" Dielectric tensor written to EPSILON_TDDFT.OUT")')
  write(*,'(" Electron loss function written to EELS_TDDFT.OUT")')
  write(*,'(" for component i, j = ",I1)') optcomp(1,1)
  write(*,'(" q-vector (lattice coordinates) : ")')
  write(*,'(3G18.10)') vecql
  write(*,'(" q-vector length : ",G18.10)') gqc(1)
end if
deallocate(vgqc,gqc,ylmgq,sfacgq)
deallocate(vchi0,vfxc,eps0,eps)
return
end subroutine
