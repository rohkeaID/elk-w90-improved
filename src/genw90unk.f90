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
use modmpi
use modomp

! !DESCRIPTION:
!
!EOP
!BOC
implicit none
! local variables
integer ikp
integer i,j,is,ispn,m
integer nrc,nrci
integer nthd
integer fileID
integer redkfil,nstsv_,recl
real(8) vgqc(3),gqc
logical exists
real(8) t1
logical :: wann_spn_up = .false., wann_spn_dn = .false.
integer wann_spn_start,wann_spn_end
! allocatable arrays
complex(8), allocatable :: wfmt(:,:,:,:),wfir(:,:,:)
complex(8), allocatable :: zwfmt(:,:,:,:)
real(8),    allocatable :: jlgqr(:,:)
complex(8), allocatable :: ylmgq(:),sfacgq(:)
complex(8), allocatable :: expmt(:,:)
integer,    allocatable :: igpig(:,:)
integer,    allocatable :: ngp(:)
real(8),    allocatable :: evalsv_(:)
! automatic arrays
character(8)            :: fmt_ikp,fmt_is ! format descriptor
real(8)                 :: bqc(3)
character(8)            :: ikp_str,is_str
character(256)          :: filename
real(8)                 :: vkl_(3)

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

! check that EVECSV.OUT has all necessary k-points
allocate(evalsv_(nstsv))
redkfil=0
inquire(iolength=recl) vkl_,nstsv_,evalsv_
do ikp=1,nkpt
  exists=.false.
  t1=9.d99
  inquire(file='EVALSV'//trim(filext),exist=exists)
  if(exists) then
    open(70,file='EVALSV'//trim(filext),action='READ',form='UNFORMATTED', &
        access='DIRECT',recl=recl,err=101)
    read(70,rec=ikp,err=101) vkl_,nstsv_,evalsv_
    close(70)
    t1=abs(vkl(1,ikp)-vkl_(1))+abs(vkl(2,ikp)-vkl_(2))+abs(vkl(3,ikp)-vkl_(3))
  end if
101 continue
  if (.not.exists.or.t1.gt.epslat.or.nstsv.ne.nstsv_) then
    redkfil=1
    exit
  end if
end do
! If kpoint not found in saved eigen-values/vectors, then need to recompute
! EVEC*OUT.
if (redkfil.ne.0) then
! compute the overlap radial integrals
  call olprad
! compute the Hamiltonian radial integrals
  call hmlrad
! generate the spin-orbit coupling radial functions
  call gensocfr
! generate the first- and second-variational eigenvectors and eigenvalues
  call genevfsv
end if


fmt_ikp = '(I5.5)' ! format of k-points in UNK filename (5 digits) 
fmt_is = '(I1.1)'

if(mp_mpi) then
  write(*,*)
  write(*,*) "Info(Wannier): Unk calculation" 
end if

! check whether both spinor projections need to be plotted
if (any( wann_proj_spin .eq. 1))  wann_spn_up = .true.
if (any( wann_proj_spin .eq. -1)) wann_spn_dn = .true.
wann_spn_start = 1; wann_spn_end = nspinor
! special cases where one spin component is enough
if(spinpol .and. .not.(spinorb) .and. wann_spn_up .and. .not.(wann_spn_dn)) wann_spn_end = 1
if(spinpol .and. .not.(spinorb) .and. .not.(wann_spn_up) .and. wann_spn_dn) then
  wann_spn_start = 2
  wann_spn_end = 2
end if

ngrf=1 ! corresponds to expmt
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
! loop over reduced k-point set
call omp_hold(nkpt/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ikp,i,j,is,ispn,nrc,nrci,m,bqc) &
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

  ! k-vector in Cartesian coordinates
  call r3mv(bvec,vkl(:,ikp),bqc)
  ! generate the phase factor function exp(ik.r) in the muffin-tins
  call gengqrf(bqc,vgqc,gqc,jlgqr,ylmgq,sfacgq)
  call genexpmt(1,jlgqr,ylmgq,1,sfacgq,expmt)

! multiply psi_nk with exp(-ikr) to get the periodic part u_nk
  do i=1,npcmtmax
    do j=1,natmtot
      do is = wann_spn_start,wann_spn_end
        do m = 1,wann_nband
          wfmt(i,j,is,m) = wfmt(i,j,is,m)*dconjg(expmt(i,j))
        end do
      end do
    end do
  end do

  do j=1,natmtot
    nrc=nrcmt(idxis(j))
    nrci=nrcmti(idxis(j))
    do is = wann_spn_start,wann_spn_end
      do m = 1,wann_nband
        call zfsht(nrc,nrci,wfmt(:,j,is,m),zwfmt(:,j,is,m))
      end do
    end do
  end do
  
! dump data to files
  ispn = 0
  fileID = 600+(ikp-1)*nspinor
  do is = wann_spn_start,wann_spn_end
    ispn = ispn + 1
    write(ikp_str,fmt_ikp) ikp 
    write(is_str,fmt_is) ispn 
    filename = 'UNK'//trim(ikp_str)//'.'//trim(is_str)
    fileID=fileID+is-1
    open(fileID,file=filename,action='WRITE',form='FORMATTED')
    
    write(fileID,*) np3d(1),np3d(2),np3d(3),ikp,wann_nband
    do m = 1,wann_nband
      call plotUNK(fileID,zwfmt(:,:,is,m),wfir(:,is,m))
    end do

    close(fileID)
  end do

!$OMP CRITICAL(genw90unk_)
  if (mp_mpi) then
    write(*,'("  Info(Wannier Unk): completed ",I6," of ",I6," k-points")')&
                                                                   &ikp,nkpt
  end if
!$OMP END CRITICAL(genw90unk_)

end do!End loop over k points
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
  write(*,*) "Info(Wannier): Unk for each k-point has been computed [ OK ]"
end if

reducek=reducek0

end subroutine
!EOC
