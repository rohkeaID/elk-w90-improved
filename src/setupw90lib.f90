! Copyright (C) 2018 Arsenii Gerasimov
! This file is distributed under the terms of the GNU General Public
! License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genw90amn
! !INTERFACE:
subroutine setupw90lib
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
character(20)  :: atomFileName
integer ia,is
logical        :: spinors_lib = .false.
 
! check that user has defined bands and projections
if(wann_projlines.eq.-1) then
  write(*,*)
  write(*,*) "Error(setupw90lib): No projections specified - please input &
                                        &using the wann_projections block."
  write(*,*)
  stop
end if
if(wann_nband.eq.-1) then
  write(*,*)
  write(*,*) "Error(setupw90lib): No bands specified - please input using &
                                                    &the wann_bands block."
  write(*,*)
  stop
end if

! initialise global variables
call init0
ngridk = wann_ngridk ! AG: should solve
call init1
call init2
call init3

! initialise variables for Wannierization
if(allocated(wann_atomsymb)) deallocate(wann_atomsymb,wann_atompos,&
                                        nnlist,nncell,wann_proj_site,&
                                        wann_proj_l,wann_proj_m,&
                                        wann_proj_radial,wann_proj_zaxis,&
                                        wann_proj_xaxis,wann_proj_zona,&
                                        wann_proj_exclude_bands_lib,&
                                        wann_proj_spin,wann_proj_quantdir,&
                                        wann_proj_isrand)
allocate(wann_atomsymb(natmtot),wann_atompos(3,natmtot))
allocate(nnlist(nkpt,num_nnmax),nncell(3,nkpt,num_nnmax))
allocate(wann_proj_site(3,nstsv))
allocate(wann_proj_l(nstsv),wann_proj_m(nstsv),wann_proj_radial(nstsv))
allocate(wann_proj_zaxis(3,nstsv),wann_proj_xaxis(3,nstsv))
allocate(wann_proj_zona(nstsv))
allocate(wann_proj_exclude_bands_lib(nstsv))
allocate(wann_proj_spin(nstsv))
allocate(wann_proj_quantdir(3,nstsv))
allocate(wann_proj_isrand(nstsv))

! prepare variables for calling lib of wannier90
do is = 1,nspecies
  do ia = 1,natoms(is)
    atomFileName = trim(spfname(is))
    wann_atomsymb(is + ia - 1) = atomFileName(1:(len(trim(atomFileName))-3)) !erase '.in'
    wann_atompos(:,is + ia - 1) = atposc(:,ia,is) * au2angstrom
  enddo ! is, loop over atoms of a species
enddo ! ia, loop over species

if ( nspinor .eq. 2 ) then
  spinors_lib = .true.
end if

! 
call wannier_setup(trim(wann_seedname),ngridk,nkpt,au2angstrom*transpose(avec),& !in
                   (1/au2angstrom)*transpose(bvec),vkl,nstsv,natmtot,&           !in
                   wann_atomsymb,wann_atompos,.false.,spinors_lib,&              !in
                   wann_nntot,nnlist,nncell,wann_nband_total,wann_nwf,&          !out
                   wann_proj_site,wann_proj_l,wann_proj_m,&                      !out
                   wann_proj_radial,wann_proj_zaxis,wann_proj_xaxis,&            !out
                   wann_proj_zona,wann_proj_exclude_bands_lib,&                  !out
                   wann_proj_spin,wann_proj_quantdir)                            !out

wann_proj_radial = 0
wann_proj_isrand = .false. ! AG: solve later
if (nspinor .eq. 2)Then
  wann_nproj = wann_nwf / 2
else
  wann_nproj = wann_nwf
endif

end subroutine setupw90lib

