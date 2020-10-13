
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
!   Generates and writes out {\tt UNKp.s} files with periodic part of
!   the Bloch functions on a regular real space grid for plotting Wannier
!   Functions in the Wannier90, where {\tt p} denotes $\bf k$-point and {\tt s}
!   corresponds to spins.
!
! !REVISION HISTORY:
!   Created August 2018 (Arsenii Gerasimov)
!EOP
!BOC
implicit none
! local variables
integer    wann_spn_start,wann_spn_end,nthd,ikp,j,is,ispn,m,nrc,nrci,fileID
real(8)    gqc
! allocatable arrays
integer,    allocatable :: ngp(:),igpig(:,:)
complex(8), allocatable :: wfmt(:,:,:,:),wfir(:,:,:),zwfmt(:,:,:,:)
real(8),    allocatable :: jlgqr(:,:)
complex(8), allocatable :: ylmgq(:),sfacgq(:),expmt(:,:)
! automatic arrays
character(8)   ::  fmt_ikp = '(I5.5)',fmt_is = '(I1.1)' ! Format descriptors
real(8)            bqc(3),vgqc(3)
character(8)       ikp_str,is_str
character(256)     filename
!-------------------------------------------------------------------------------

! Generate global variables for the Wannier90 interface
call initw90

! Generate wann_spn_start and wann_spn_end
call getw90ispin

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
!$OMP PRIVATE(ikp,j,is,ispn,nrc,nrci,m,bqc,ikp_str,is_str,filename,ngp,igpig) &
!$OMP PRIVATE(fileID,wfmt,wfir,zwfmt,vgqc,gqc,jlgqr,ylmgq,sfacgq,expmt) &
!$OMP NUM_THREADS(nthd)

allocate(ngp(nspnfv),igpig(ngkmax,nspnfv))
allocate(wfmt(npcmtmax,natmtot,nspinor,wann_nband))
allocate(wfir(ngtot,nspinor,wann_nband))
allocate(zwfmt(npcmtmax,natmtot,nspinor,wann_nband))
allocate(jlgqr(njcmax,nspecies),ylmgq(lmmaxo),sfacgq(natmtot))
allocate(expmt(npcmtmax,natmtot))

!$OMP DO
do ikp = 1,nkpt

  call genwfsvp(.false.,.false.,wann_nband,wann_bands,ngridg,igfft,vkl(:,ikp),&
                                                      ngp,igpig,wfmt,ngtot,wfir)

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

  !wfir = zzero

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
  fileID = 500 + ( ikp - 1 )*nspinor
  do is = wann_spn_start,wann_spn_end
    ispn = ispn + 1
    write(ikp_str,fmt_ikp) ikp; write(is_str, fmt_is ) ispn
    filename = 'UNK'//trim(ikp_str)//'.'//trim(is_str)
    fileID = fileID + is - 1
    open(fileID,file=filename,action='WRITE',form='FORMATTED')
    
    ! AG_UNK
    !write(fileID,*) np3d(1)+1,np3d(2)+1,np3d(3)+1,ikp,wann_nband
    write(fileID,*) np3d(:),ikp,wann_nband
    do m = 1,wann_nband
      call plotw90unk(fileID,zwfmt(:,:,is,m),wfir(:,is,m))
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

deallocate(ngp,igpig,ylmgq,sfacgq,jlgqr,expmt,zwfmt,wfmt,wfir)

!$OMP END PARALLEL
call omp_free(nthd)

if( allocated(wann_atomsymb) ) deallocate(wann_atomsymb,wann_atompos,       &
                                          nnlist,nncell,wann_proj_site,     &
                                          wann_proj_l,wann_proj_m,          &
                                          wann_proj_radial,wann_proj_zaxis, &
                                          wann_proj_xaxis,wann_proj_zona,   &
                                          wann_proj_exclude_bands_lib,      &
                                          wann_proj_spin,wann_proj_quantdir,&
                                          wann_proj_isrand)

if ( mp_mpi ) then
  write(*,*)
  write(*,*) " Info(Wannier): Unk for each k-point has been computed [ OK ]"
end if

reducek = reducek0

return


! LOCAL: subroutines
contains


! LOCAL: subroutine getw90ispin
subroutine getw90ispin
implicit none
! local variables
logical :: wann_spn_up = .false., wann_spn_dn = .false.
!-------------------------------------------------------------------------------
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

end subroutine getw90ispin

end subroutine
!EOC
