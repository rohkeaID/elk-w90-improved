! Copyright (C) 2018 Arsenii Gerasimov
! This file is distributed under the terms of the GNU General Public
! License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genw90unk
! !INTERFACE:
subroutine genw90unk
! !USES:
use modmain
use modw90
use modw90overlap
use modmpi
use modomp

! !DESCRIPTION:
!
!EOP
!BOC
implicit none
! local variables
integer ikp
integer i,j,is,m
integer nrc,nrci
integer nthd
integer fileID
real(8) vgqc(3),gqc
! allocatable arrays
complex(8), allocatable :: wfmt(:,:,:,:),wfir(:,:,:)
complex(8), allocatable :: zwfmt(:,:,:,:)
real(8),    allocatable :: jlgqr(:,:)
complex(8), allocatable :: ylmgq(:),sfacgq(:)
complex(8), allocatable :: expmt(:,:)
integer,    allocatable :: igpig(:,:)
integer,    allocatable :: ngp(:)
! automatic arrays
character(8)            :: fmt_ikp,fmt_is ! format descriptor
real(8)                 :: bqc(3)
character(8)            :: ikp_str,is_str
character(256)          :: filename

reducek0=reducek ! if reducek=1 was used in ground state calculations,
                 ! need to regenerate the eigenvectors set for the full BZ.
reducek=0

call setupw90lib

! read density and potentials from file
call readstate
! Fourier transform Kohn-Sham potential to G-space
call genvsig
! read Fermi energy from a file
call readfermi
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr

fmt_ikp = '(I5.5)' ! format of k-points in UNK filename an integer of width 5 with zeros at the left
fmt_is = '(I1.1)'

if(mp_mpi) then
  write(*,*)
  write(*,*) "Info(Wannier): Unk calculation"
end if

ngrf=1 ! corresponds to expmt
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
! loop over reduced k-point set
call omp_hold(nkpt/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ikp,i,j,is,nrc,nrci,m,bqc) &
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
do ikp=1,nkpt

  call genwfsvp(.false.,.false.,wann_nband,wann_bands,ngridg,igfft,vkl(:,ikp),ngp,igpig,wfmt,ngtot,wfir) 

  ! b-vector in Cartesian coordinates
  call r3mv(bvec,vkl(:,ikp),bqc)
  ! generate the phase factor function exp(ib.r) in the muffin-tins
  call gengqrf(bqc,vgqc,gqc,jlgqr,ylmgq,sfacgq)
  call genexpmt(1,jlgqr,ylmgq,1,sfacgq,expmt)

  do i=1,npcmtmax
    do j=1,natmtot
      do is = 1,nspinor
        do m = 1,wann_nband
          wfmt(i,j,is,m) = wfmt(i,j,is,m)*dconjg(expmt(i,j))
        end do
      end do
    end do
  end do

  do j=1,natmtot
    nrc=nrcmt(idxis(j))
    nrci=nrcmti(idxis(j))
    do is = 1,nspinor
      do m = 1,wann_nband
        call zfsht(nrc,nrci,wfmt(:,j,is,m),zwfmt(:,j,is,m))
      end do
    end do
  end do

  do is = 1,nspinor
    write(ikp_str,fmt_ikp) ikp ! converting integer to string using an 'internal file'
    write(is_str,fmt_is) is ! converting integer to string using an 'internal file'
    filename = 'UNK'//trim(ikp_str)//'.'//trim(is_str)
    fileID=600+ikp+is
    open(fileID,file=filename,action='WRITE',form='FORMATTED')
    
    write(fileID,*) np3d(1),np3d(2),np3d(3),ikp,wann_nband
    do m = 1,wann_nband
      call plotUNK(fileID,zwfmt(:,:,is,m),wfir(:,is,m))
    end do

    close(fileID)
  !$OMP CRITICAL(genw90unk_)
    if (mp_mpi) then
      write(*,'("  Info(Wannier Unk): completed ",I6," of ",I6," k-points")')&
                                                                     &ikp,nkpt
    end if
  !$OMP END CRITICAL(genw90unk_)

  end do

end do
!$OMP END DO

deallocate(ngp)
deallocate(igpig)
deallocate(expmt)
deallocate(ylmgq,sfacgq,jlgqr)
deallocate(zwfmt,wfmt,wfir)

!$OMP END PARALLEL
call omp_free(nthd)

if (mp_mpi) then
  write(*,*)
  write(*,*) "Info(Wannier): Unk for each k-point are computed"
end if

reducek=reducek0

end subroutine
!EOC

