! Copyright (C) 2015 Jon Lafuente and Manh Duc Le
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writew90mmn
! !INTERFACE:
subroutine genw90input
! !USES:
use modmain
use modw90
use modw90overlap
! !DESCRIPTION:
!   Writes the Mmn and Amn matrices required for Wannier90 to file.
!
! !REVISION HISTORY:
!   Created March 2015 (Jon Lafuente and Manh Duc Le)
!   Upgraded July 2017 (Arsenii Gerasimov)
!EOP
!BOC

implicit none
! local variables
character(20) :: atomFileName
integer ik,jk,n,m,ig,is,ias,i
integer innkp,ikp,inn,ispn
integer nrc,nrci,irco,ia,lm,l,ist,ir
integer reducek0,redkfil,nstsv_,recl
logical exists
real(8) t1
complex(8) :: z1
integer jst
! allocatable arrays
complex(8), allocatable :: expmt(:,:)
complex(8), allocatable :: mmn(:,:,:,:)
complex(8), allocatable :: amn(:,:,:,:)
complex(8), allocatable :: spn_x(:,:),spn_y(:,:),spn_z(:,:)
complex(8), allocatable :: wfmt(:,:,:,:),wfir(:,:,:)
complex(8), allocatable :: wfmtq(:,:,:,:),wfirq(:,:,:)
complex(8), allocatable :: twfmt(:,:,:,:),twfir(:,:,:),twfmt1(:)
real(8), allocatable :: evalsv_(:)
! automatic arrays
real(8) :: bqvec(3),bqc(3),vkql(3)
character(256) filename
character(10) dat,tim
real(8) vkl_(3),vec_0(3)

! Checks that user has defined bands and projections
if(wann_projlines.eq.-1) then
  write(*,*)
  write(*,'("Error(writew90mmn): No projections specified - please input using the wann_projections block.")')
  write(*,*)
  stop
end if
if(wann_nband.eq.-1) then
  write(*,*)
  write(*,'("Error(writew90mmn): No bands specified - please input using the wann_bands block.")')
  write(*,*)
  stop
end if

! initialise global variables
call init0
! if reducek=1 was used in ground state calculations, need to regenerate the eigenvectors set for the full BZ.
reducek0=reducek
if (reducek0.ne.0) reducek=0
call init1
call init2
call init3

! get the nearest neighbour k-points and projections
wann_natoms = 0
do is = 1,nspecies
  wann_natoms = wann_natoms + natoms(is)
enddo ! ia, loop over species
allocate(wann_atomsymb(wann_natoms))
allocate(wann_atompos(3,wann_natoms))
do is = 1,nspecies
  do ia = 1,natoms(is)
    atomFileName = trim(spfname(is))
    wann_atomsymb(is + ia - 1) = atomFileName(1:(len(trim(atomFileName))-3)) ! erase '.in'
    ! write(*,*) "wann_atomsymb(",is + ia - 1,") = ", atomFileName(1:(len(trim(atomFileName))-3))
    wann_atompos(:,is + ia - 1) = atposc(:,ia,is) * bohr2angstrom
  enddo ! is, loop over atoms of a species
enddo ! ia, loop over species
! allocate(wann_nnkp(nkpt,wann_nntot)) ! wann_nntot - поменять на num_nnmax, переписать логику под новый wann_nnkp
! if(allocated(wann_nnkp))     deallocate(wann_nnkp)
! allocate(wann_nnkp(5,wann_nntot*nkpt))
allocate(nnlist_lib(nkpt,num_nnmax))
allocate(nncell_lib(3,nkpt,num_nnmax))
if(allocated(wann_proj_site_lib))     deallocate(wann_proj_site_lib)
allocate(wann_proj_site_lib(3,wann_nband))
if(allocated(wann_proj_l_lib))     deallocate(wann_proj_l_lib)
allocate(wann_proj_l_lib(wann_nband))
if(allocated(wann_proj_m_lib))     deallocate(wann_proj_m_lib)
allocate(wann_proj_m_lib(wann_nband))
if(allocated(wann_proj_zaxis_lib))     deallocate(wann_proj_zaxis_lib)
allocate(wann_proj_zaxis_lib(3,wann_nband))
if(allocated(wann_proj_xaxis_lib))     deallocate(wann_proj_xaxis_lib)
allocate(wann_proj_xaxis_lib(3,wann_nband))
if(allocated(wann_proj_radial_lib))     deallocate(wann_proj_radial_lib)
allocate(wann_proj_radial_lib(wann_nband))
if(allocated(wann_proj_zona_lib))     deallocate(wann_proj_zona_lib)
allocate(wann_proj_zona_lib(wann_nband))
if(allocated(wann_proj_exclude_bands_lib))     deallocate(wann_proj_exclude_bands_lib)
allocate(wann_proj_exclude_bands_lib(wann_nband))
if(allocated(wann_proj_spin_lib))     deallocate(wann_proj_spin_lib)
allocate(wann_proj_spin_lib(wann_nband))
if(allocated(wann_proj_quantdir_lib))     deallocate(wann_proj_quantdir_lib)
allocate(wann_proj_quantdir_lib(3,wann_nband))

call wannier_setup(trim(wann_seedname),ngridk,nkpt,bohr2angstrom*avec,&
                   (1/bohr2angstrom)*bvec,vkl,wann_nband,wann_natoms,&
                   wann_atomsymb,wann_atompos,.false.,spinpol,wann_nntot,&
                   nnlist_lib,nncell_lib,wann_nband_out,wann_nwf_lib,&
                   wann_proj_site_lib,&
                   wann_proj_l_lib,wann_proj_m_lib,wann_proj_radial_lib,&
                   wann_proj_zaxis_lib,&
                   wann_proj_xaxis_lib,wann_proj_zona_lib,&
                   wann_proj_exclude_bands_lib,&
                   wann_proj_spin_lib,wann_proj_quantdir_lib)

! get the nearest neighbour k-points and projections
! call getw90nnkp
! call getw90proj
if(allocated(wann_proj_isrand))   deallocate(wann_proj_isrand)
allocate(wann_proj_isrand(wann_nband))
wann_proj_isrand = .false.

! write(*,*) 'wann_proj_isrand: '
! do jst = 1,wann_nproj
!   write(*,*) 'jst = ', jst, 'wann_proj_isrand = ', wann_proj_isrand(jst)
! enddo

! Open files for writting
filename = trim(wann_seedname)//'.mmn'
open(500,file=filename,action='WRITE',form='FORMATTED')
call date_and_time(date=dat,time=tim)
write(500,'("Generated by ELK on ",A4,"-",A2,"-",A2," at ",A2,":",A2,":",A2)')&
    dat(1:4),dat(5:6),dat(7:8),tim(1:2),tim(3:4),tim(5:6)
write(500,'(3I8)') wann_nband,nkpt,wann_nntot
filename = trim(wann_seedname)//'.amn'
open(501,file=filename,action='WRITE',form='FORMATTED')
call date_and_time(date=dat,time=tim)
write(501,'("Generated by ELK on ",A4,"-",A2,"-",A2," at ",A2,":",A2,":",A2)')&
    dat(1:4),dat(5:6),dat(7:8),tim(1:2),tim(3:4),tim(5:6)
if (nspinor.eq.1) then
  write(501,'(3I8)') wann_nband,nkpt,wann_nproj
else
  write(501,'(3I8)') wann_nband,nkpt,2*wann_nproj
  filename = trim(wann_seedname)//'.spn'
  open(502,file=filename,action='WRITE',form='FORMATTED')
  call date_and_time(date=dat,time=tim)
  write(502,'("Generated by ELK on ",A4,"-",A2,"-",A2," at ",A2,":",A2,":",A2)')&
      dat(1:4),dat(5:6),dat(7:8),tim(1:2),tim(3:4),tim(5:6)
  write(502,'(3I8)') wann_nband,nkpt
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
! generates the trial wavefunctions from the projection definitions read in
call genw90twf

! creates the trial wavefunction matrix
allocate(twfmt1(npcmtmax))
allocate(twfmt(npcmtmax,natmtot,nspinor,wann_nproj))
allocate(twfir(ngtot,nspinor,wann_nproj))
twfmt = cmplx(0.d0,0.d0,kind=8)
twfir = cmplx(0.d0,0.d0,kind=8)     ! trial wavefunction is zero outside muffin-tin.
do ispn=1,nspinor
  do is=1,nspecies
    nrc=nrcmt(is)
    nrci=nrcmti(is)
    irco=nrci+1
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      do n=1,wann_nproj
        if(.not.wann_proj_haswt(n,ias)) cycle
        ! zero the wavefunction
        twfmt1 = cmplx(0.d0,0.d0,kind=8)
        i=1
        do ir=1,nrci
          lm=0
          do l=0,lmaxi
            do m=-l,l
              lm=lm+1
              z1=wann_projclm(lm,ias,n)
              if(abs(z1).gt.1.d-14) then
                twfmt1(i)=twfmt1(i)+z1*wann_projulr(ir,l,ias,n)
              end if
              i=i+1
            end do
          end do
        end do
        do ir=nrci+1,nrc
          lm=0
          do l=0,lmaxo
            do m=-l,l
              lm=lm+1
              z1=wann_projclm(lm,ias,n)
              if(abs(z1).gt.1.d-14) then
                twfmt1(i)=twfmt1(i)+z1*wann_projulr(ir,l,ias,n)
              end if
              i=i+1
            end do
          end do
        end do
        call zbsht(nrc,nrci,twfmt1,twfmt(:,ias,ispn,n))
      end do
    end do
  end do
end do
deallocate(twfmt1)

! check that EVECSV.OUT has all necessary k-points
allocate(evalsv_(nstsv))
redkfil=0
inquire(iolength=recl) vkl_,nstsv_,evalsv_
do ik=1,nkpt
  exists=.false.
  t1=9.d99
  inquire(file='EVALSV'//trim(filext),exist=exists)
  if(exists) then
    open(70,file='EVALSV'//trim(filext),action='READ',form='UNFORMATTED', &
        access='DIRECT',recl=recl,err=101)
    read(70,rec=ik,err=101) vkl_,nstsv_,evalsv_
    close(70)
    t1=abs(vkl(1,ik)-vkl_(1))+abs(vkl(2,ik)-vkl_(2))+abs(vkl(3,ik)-vkl_(3))
  end if
101 continue
  if (.not.exists.or.t1.gt.epslat.or.nstsv.ne.nstsv_) then
    redkfil=1
    exit
  end if
end do
! If kpoint not found in saved eigen-values/vectors, then need to recompute EVEC*OUT.
if (redkfil.ne.0) then
  write(*,'("Info(writew90mmn): saved k-points do not contain all required k-points. &
           &Recalculating wavefunctions")')
! compute the overlap radial integrals
  call olprad
! compute the Hamiltonian radial integrals
  call hmlrad
! generate the spin-orbit coupling radial functions
  call gensocfr
! generate the first- and second-variational eigenvectors and eigenvalues
  call genevfsv
end if

! write .eig file for wannier90
filename = trim(wann_seedname)//'.eig'
open(50,file=filename,action='WRITE',form='FORMATTED')
do ik=1,nkpt
  call getevalsv(filext,ik,vkl(:,ik),evalsv_)
  do ig=1,wann_nband
    ist=wann_bands(ig)
    write(50,'(2I12,G18.10)') ig,ik,(evalsv_(ist)-efermi)*27.211385  ! convert to eV
  end do
end do
close(50)
deallocate(evalsv_)

if(nspinor.eq.2) then
  allocate(spn_x(nkpt,(wann_nband*(wann_nband+1))/2))
  allocate(spn_y(nkpt,(wann_nband*(wann_nband+1))/2))
  allocate(spn_z(nkpt,(wann_nband*(wann_nband+1))/2))
end if

!Loop over k and k+b points
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ikp,innkp,inn,ik,jk,bqvec,mmn,amn,vkql,wfmt,wfir,wfmtq,wfirq,expmt)

allocate(wfmt(npcmtmax,natmtot,nspinor,wann_nband))
allocate(wfir(ngtot,nspinor,wann_nband))
allocate(expmt(npcmtmax,natmtot))
allocate(wfmtq(npcmtmax,natmtot,nspinor,wann_nband))
allocate(wfirq(ngtot,nspinor,wann_nband))
allocate(mmn(wann_nband,4,wann_nband,4))
allocate(amn(wann_nband,4,wann_nproj,4))

!$OMP DO
do ikp=1,nkpt

  call genwfsvp(.false.,.false.,wann_nband,wann_bands,vkl(:,ikp),wfmt,ngtot,wfir)

  do inn=1,wann_nntot
    innkp=inn+(ikp-1)*wann_nntot
    ! ik = wann_nnkp(1,innkp)
    ik = ikp
    ! jk=wann_nnkp(2,innkp)
    jk = nnlist_lib(ik,inn)
    ! bqvec=wann_nnkp(3:5,innkp)+vkl(:,jk)-vkl(:,ik)
    bqvec=nncell_lib(1:3,ik,inn)+vkl(:,jk)-vkl(:,ik)

    ! b-vector in Cartesian coordinates
    call r3mv(bvec,bqvec,bqc)
    ! generate the phase factor function exp(ib.r) in the muffin-tins
    call genexpmt(bqc,expmt)

    ! k+b-vector in lattice coordinates
    vkql = vkl(:,ik)+bqvec
    call genwfsvp(.false.,.false.,wann_nband,wann_bands,vkql,wfmtq,ngtot,wfirq)

    ! compute Mmn
    mmn = cmplx(0.d0,0.d0,kind=8)
    call genw90overlap(wfmt,wfir,wann_nband,wfmtq,wfirq,mmn,expmt)

!$OMP CRITICAL
    !Write the Mmn matrix elements
    ! write(500,'(5I8)') wann_nnkp(:,innkp)
    write(500,'(5I8)') ikp,nnlist_lib(ikp,inn),nncell_lib(:,ikp,inn)
    do n=1,wann_nband
      do m=1,wann_nband
        if (nspinor.eq.1) then
          write(500,'(2G18.10)') dble(mmn(m,1,n,1)),aimag(mmn(m,1,n,1))
        else
          write(500,'(2G18.10)') dble(mmn(m,1,n,1)+mmn(m,2,n,2)),aimag(mmn(m,1,n,1)+mmn(m,2,n,2))
        end if
      end do
    end do
    write(*,'("Info(writew90mmn): completed ",I5," of ",I5," Mmn(k,k+b) points")')&
        innkp,wann_nntot*nkpt
!$OMP END CRITICAL

  end do !End loop over b points


  ! Calculates the overlap integrals Amn(k)
  amn = cmplx(0.d0,0.d0,kind=8)
  call genw90overlap(wfmt,wfir,wann_nproj,twfmt,twfir,amn)
!$OMP CRITICAL
  if (nspinor.eq.1) then
    do n=1,wann_nproj
      do m=1,wann_nband
          write(501,'(3I8,2G18.10)') m,n,ikp,dble(amn(m,1,n,1)),aimag(amn(m,1,n,1))
      end do
    end do
  else
    ! do n=1,2*wann_nproj,2
    do n=1,wann_nproj
      do m=1,wann_nband
          write(501,'(3I8,2G18.10)') m,2*n-1,ikp,dble(amn(m,1,n,1)),aimag(amn(m,1,n,1))
          write(501,'(3I8,2G18.10)') m,2*n,ikp,dble(amn(m,2,n,1)),aimag(amn(m,2,n,1))
      end do
    end do
  end if
  write(*,'("Info(writew90amn): completed ",I5," of ",I5," Amn(k) points")')&
      ik,nkpt
!$OMP END CRITICAL

  !Write the .spn file
  if(nspinor.eq.2) then
    vec_0(1)=0.d0
    vec_0(2)=0.d0
    vec_0(3)=0.d0
    call genexpmt(vec_0,expmt)
    mmn = cmplx(0.d0,0.d0,kind=8)
    call genw90overlap(wfmt,wfir,wann_nband,wfmt,wfir,mmn,expmt)
  !$OMP CRITICAL
    is=1
    do m=1,wann_nband
      do n=1,m
        spn_x(ikp,is) = mmn(n,1,m,2)+mmn(n,2,m,1)
        spn_y(ikp,is) = cmplx(0.d0,1.d0,kind=8)*(mmn(n,2,m,1)-mmn(n,1,m,2))
        spn_z(ikp,is) = mmn(n,1,m,1)-mmn(n,2,m,2)
        is=is+1
      end do
    end do

  !$OMP END CRITICAL

  end if

end do !End loop over k points
!$OMP END DO

deallocate(wfmt,wfir)
deallocate(expmt)
deallocate(wfmtq,wfirq)
deallocate(mmn)
deallocate(amn)

!$OMP END PARALLEL

close(500)
close(501)

if(nspinor.eq.2) then
  do ikp=1,nkpt
    do is=1,(wann_nband*(wann_nband+1))/2
      write(502,'(2G18.10)') spn_x(ikp,is)
      write(502,'(2G18.10)') spn_y(ikp,is)
      write(502,'(2G18.10)') spn_z(ikp,is)
    end do
  end do
  close(502)
  deallocate(spn_x,spn_y,spn_z)
end if

reducek=reducek0

end subroutine genw90input
!EOC
