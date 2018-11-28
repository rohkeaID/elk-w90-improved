
! Copyright (C) 2017-18 Arsenii Gerasimov, Yaroslav Kvashnin and Lars Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genw90unk
! !INTERFACE:
subroutine genw90unk
! !USES:
use modmain
use modw90
use modmpi
use modomp
! !DESCRIPTION:
!   
!
! !REVISION HISTORY:
!   Created August 2018 (Arsenii Gerasimov)
!EOP
!BOC
implicit none
! local variables
logical :: wann_spn_up = .false., wann_spn_dn = .false.
integer    nthd
integer    ikp
integer    j,is,ispn,m
integer    nrc,nrci
integer    fileID
integer    wann_spn_start,wann_spn_end
real(8)    gqc
! allocatable arrays
integer,    allocatable :: igpig(:,:)
integer,    allocatable :: ngp(:)
real(8),    allocatable :: jlgqr(:,:)
complex(8), allocatable :: wfmt(:,:,:,:),wfir(:,:,:)
complex(8), allocatable :: zwfmt(:,:,:,:)
complex(8), allocatable :: ylmgq(:),sfacgq(:)
complex(8), allocatable :: expmt(:,:)
! automatic arrays
real(8)        bqc(3),vgqc(3)
character(8)   fmt_ikp,fmt_is ! Format descriptor
character(8)   ikp_str,is_str
character(256) filename
!-------------------------------------------------------------------------------
!
call initw90

fmt_ikp = '(I5.5)' ! Format of k-points   in UNK filename (5 digits)
fmt_is  = '(I1.1)' ! Format of spin index in UNK filename (1 digit )

wann_spn_start = 1; wann_spn_end = nspinor
! Check whether both spinor projections need to be plotted
if ( any( wann_proj_spin .eq.  1) ) wann_spn_up = .true.
if ( any( wann_proj_spin .eq. -1) ) wann_spn_dn = .true.
! Special cases where one spin component is enough
if( spinpol .and. .not.(spinorb) ) then
  if ( wann_spn_up ) then
    if ( .not.(wann_spn_dn) ) wann_spn_end = 1
  else
    if ( wann_spn_dn ) then
      wann_spn_start = 2
      wann_spn_end   = 2
    end if
  end if
end if

if( mp_mpi ) then
  write(*,*)
  write(*,*) " Info(Wannier): Unk calculation"
end if

ngrf = 1 ! Corresponds to expmt
! Synchronise MPI processes
call mpi_barrier(mpicom,ierror)
! Loop over reduced k-point set
call omp_hold(nkpt/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ikp,j,is,ispn,nrc,nrci,m,bqc) &
!$OMP PRIVATE(ikp_str,is_str,filename,ngp,igpig) &
!$OMP PRIVATE(fileID,wfmt,wfir,zwfmt) &
!$OMP PRIVATE(vgqc,gqc,jlgqr,ylmgq,sfacgq,expmt) &
!$OMP NUM_THREADS(nthd)

allocate(wfmt(npcmtmax,natmtot,nspinor,wann_nband))
allocate(wfir(ngtot,nspinor,wann_nband))
allocate(zwfmt(npcmtmax,natmtot,nspinor,wann_nband))
allocate(jlgqr(njcmax,nspecies))
allocate(ylmgq(lmmaxo),sfacgq(natmtot))
allocate(expmt(npcmtmax,natmtot))
allocate(igpig(ngkmax,nspnfv))
allocate(ngp(nspnfv))

!$OMP DO
do ikp = 1,nkpt

  call genwfsvp(.false.,.false.,wann_nband,wann_bands,ngridg,igfft,vkl(:,ikp),ngp,igpig,wfmt,ngtot,wfir)

  ! k-vector in Cartesian coordinates
  call r3mv(bvec,vkl(:,ikp),bqc)
  ! Generate phase factor function exp(ikr) in the muffin-tins
  call gengqrf(bqc,vgqc,gqc,jlgqr,ylmgq,sfacgq)
  call genexpmt(1,jlgqr,ylmgq,1,sfacgq,expmt)

  ! Multiply psi_nk with exp(-ikr) to get the periodic part u_nk
  do m = 1,wann_nband
    do is = wann_spn_start,wann_spn_end
      wfmt(:,:,is,m) = wfmt(:,:,is,m)*dconjg(expmt(:,:))
    end do
  end do

  do j = 1,natmtot
    nrc = nrcmt(idxis(j))
    nrci = nrcmti(idxis(j))
    do is = wann_spn_start,wann_spn_end
      do m = 1,wann_nband
        call zfsht(nrc,nrci,wfmt(:,j,is,m),zwfmt(:,j,is,m))
      end do
    end do
  end do

  ! Dump data to files
  ispn = 0
  fileID = 600 + (ikp-1)*nspinor
  do is = wann_spn_start,wann_spn_end
    ispn = ispn + 1
    write(ikp_str,fmt_ikp) ikp
    write(is_str,fmt_is) ispn
    filename = 'UNK'//trim(ikp_str)//'.'//trim(is_str)
    fileID = fileID + is -1
    open(fileID,file=filename,action='WRITE',form='FORMATTED')

    write(fileID,*) np3d(1),np3d(2),np3d(3),ikp,wann_nband
    do m = 1,wann_nband
      call plotUNK(fileID,zwfmt(:,:,is,m),wfir(:,is,m))
    end do

    close(fileID)
  end do

!$OMP CRITICAL(genw90unk_)
  if ( mp_mpi ) then
    write(*,'("    Info(Wannier Unk): completed ",I6," of ",I6," k-points")')&
                                                                     &ikp,nkpt
  end if
!$OMP END CRITICAL(genw90unk_)

end do ! End loop over k points
!$OMP END DO

deallocate(ngp)
deallocate(igpig)
deallocate(expmt)
deallocate(ylmgq,sfacgq,jlgqr)
deallocate(zwfmt,wfmt,wfir)

!$OMP END PARALLEL
call omp_free(nthd)

if ( mp_mpi ) then
  write(*,*)
  write(*,*) " Info(Wannier): Unk for each k-point has been computed [ OK ]"
end if

reducek = reducek0

end subroutine
!EOC
