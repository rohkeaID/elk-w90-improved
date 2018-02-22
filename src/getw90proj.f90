! Copyright (C) 2015 Manh Duc Le
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: getw90proj
! !INTERFACE:
subroutine getw90proj
! !USES:
use modmain
use modw90
! !DESCRIPTION:
!   Parses the wann\_projections input block to determine the projection
!   orbitals for the calculation of the Amn matrix required by Wannier90
!
! !REVISION HISTORY:
!   Created July 2015 (Manh Duc Le)
!EOP
!BOC
implicit none
! local variables
integer il,l1,l2,isc,iob,ipr,is,ia,ias,i
logical lastpass
character(256) element, orbitals, orb, options, opt, lv, rv
character(10) lstr(44)
integer lmul(44),lval(44),mval(44),l,mr
data lstr /'s','pz','px','py','dz2','dxz','dyz','dx2-y2','dxy', &
           'fz3','fxz2','fyz2','fz(x2-y2)','fxyz','fx(x2-3y2)','fy(3x2-y2)', &
           'sp-1','sp-2','sp2-1','sp2-2','sp2-3','sp3-1','sp3-2','sp3-3','sp3-4', &
           'sp3d-1','sp3d-2','sp3d-3','sp3d-4','sp3d-5', &
           'sp3d2-1','sp3d2-2','sp3d2-3','sp3d2-4','sp3d2-5','sp3d2-6', &
           'p','d','f','sp','sp2','sp3','sp3d','sp3d2'/
data lmul /1, 1,1,1, 1,1,1,1,1, 1,1,1,1,1,1,1, 1,1, 1,1,1, 1,1,1,1, &
           1,1,1,1,1, 1,1,1,1,1,1, 3,5,7,2,3,4,5,6/
data lval /0, 1,1,1, 2,2,2,2,2, 3,3,3,3,3,3,3, -1,-1, -2,-2,-2, -3,-3,-3,-3, &
           -4,-4,-4,-4,-4, -5,-5,-5,-5,-5,-5, 1,2,3,-1,-2,-3,-4,-5/
data mval /1, 1,2,3, 1,2,3,4,5, 1,2,3,4,5,6,7, 1,2, 1,2,3, 1,2,3,4, &
           1,2,3,4,5, 1,2,3,4,5,6, 0,0,0,0,0,0,0,0/
integer lmr(3,100),rval,ncen
real(8) zx(6,100),zona(100),zonaval,xval(3),zval(3),cen(3,100),cen0(3,100)
logical lrand(100)

! Parse wann_projstr - must be of form (Element):(orbital):(optional_parts)
! First determines the number of wann_nproj
wann_nproj = 0
lrand = .false.
do il=1,wann_projlines
! Checks if this should be a random projection
  if(trim(adjustl(wann_projstr(il))).eq.'random') then
    wann_nproj = wann_nproj + 1
    lrand(wann_nproj) = .true.
    cycle
  end if
! Reads the element/site name
  l1 = index(wann_projstr(il),":")+1
  read(wann_projstr(il)(1:l1-2),'(A)') element
  isc = index(element,"=")
  if(isc.ne.0) then
! Gives coordinates
    lv = element(:isc-1)
    rv = element(isc+1:)
    if(trim(adjustl(lv)).eq.'f'.or.trim(adjustl(lv)).eq.'c') then
      read(rv,*) cen0(:,1)
      ncen = 1
    else
      write(*,*)
      write(*,'("Error(getw90proj): Unrecognised site notation: ",A)') trim(adjustl(element))
      write(*,*)
      stop
    end if
  else
    do is=1,nspecies
      if(trim(spsymb(is)).eq.trim(adjustl(element))) then
        ncen = natoms(is)
        do ia=1,natoms(is)
          ias=idxas(ia,is)
! code taken from writegeom.f90
          if (molecule) then
! map lattice coordinates to [-0.5,0.5)
            cen0(:,ia)=atposl(:,ia,is)
            do i=1,3
              if (cen0(i,ia).gt.0.5d0) cen0(i,ia)=cen0(i,ia)-1.d0
            end do
          else
! otherwise write lattice coordinates
            cen0(:,ia)=atposl(:,ia,is)
          end if
        end do
        exit
      end if
    end do
    if(is.gt.nspecies) then
      write(*,*)
      write(*,'("Error(getw90proj): Unrecognised site notation: ",A)') trim(adjustl(element))
      write(*,*)
      stop
    end if
! Projection centered on a named site
  end if
! Reads the angular momentum of the orbital
  l2 = index(wann_projstr(il)(l1:),":")+l1
  if (l2.eq.l1) l2 = len(wann_projstr(il))
  read(wann_projstr(il)(l1:l2-2),'(A256)') orbitals
  orbitals = trim(orbitals)
! Reads the optional parts (zaxis,xaxis,zona,radial) - split by colons
  rval = 0
  zonaval = 1.d0
  xval = (/1.d0,0.d0,0.d0/)
  zval = (/0.d0,0.d0,1.d0/)
  read(wann_projstr(il)(l2:),'(A256)') options
  lastpass = .false.
  if(trim(options).ne.'') then
    do
      isc = index(options,":")
      opt = options(:isc-1)
      options = options(isc+1:)
      if(isc.eq.0) then
        opt = options
        lastpass = .true.
      end if
      isc = index(opt,"=")
      lv = opt(:isc-1)
      rv = opt(isc+1:)
      if(lv.eq.'zona') then
        read(rv,'(F9.0)') zonaval
      else if(lv.eq.'r') then
        read(rv,'(I9)') rval
      elseif(lv.eq.'x') then
        read(rv,'(3F9.0)') xval
      elseif(lv.eq.'z') then
        read(rv,'(3F9.0)') zval
      else
        write(*,*)
        write(*,'("Error(getw90proj): Unrecognised option: ",A)') trim(lv)
        write(*,*)
        stop
      end if
      if(lastpass) exit
    end do
  end if
! Splits orbital by semicolons
  lastpass = .false.
  do
    isc = index(orbitals,";")
    orb = orbitals(:isc-1)
    orbitals = orbitals(isc+1:)
    if(isc.eq.0) then
       orb = orbitals
       lastpass = .true.
    end if
    isc = index(orb,"=")
    if (isc.eq.0) then         ! String type (e.g. sp3, f)
      do iob=1,44
        if(trim(orb).eq.lstr(iob)) then
          do is=1,ncen
            if(lmul(iob).eq.1) then
              wann_nproj = wann_nproj + 1
              cen(:,wann_nproj) = cen0(:,is)
              lmr(:,wann_nproj) = (/lval(iob),mval(iob),rval/)
              zx(:,wann_nproj) = (/xval,zval/)
              zona(wann_nproj) = zonaval
            else
              do ipr=1,lmul(iob)
                wann_nproj = wann_nproj + 1
                cen(:,wann_nproj) = cen0(:,is)
                lmr(:,wann_nproj) = (/lval(iob),ipr,rval/)
                zx(:,wann_nproj) = (/xval,zval/)
                zona(wann_nproj) = zonaval
              end do
            end if
          end do
          exit
        end if
      end do
      if(iob.gt.44) then
        write(*,*)
        write(*,'("Error(getw90proj): Unrecognised projection: ",A)') trim(orb)
        write(*,*)
        stop
      end if
    else                       ! list of l and/or mr. (e.g. l=2)
      isc = index(orb,",")
      if (isc.eq.0) then       ! just "l="
        isc = index(orb,"l=")
        if(isc.eq.0) then
          write(*,*)
          write(*,'("Error(getw90proj): No ang_mtm given in: ",A)') trim(wann_projstr(il))
          write(*,*)
          stop
        end if
        read(orb(isc+2:),'(I9)') l
        do is=1,ncen
          do ipr=1,2*l+1
            wann_nproj = wann_nproj + 1
            cen(:,wann_nproj) = cen0(:,is)
            lmr(:,wann_nproj) = (/l,ipr,rval/)
            zx(:,wann_nproj) = (/xval,zval/)
            zona(wann_nproj) = zonaval
          end do
        end do
      else                     ! "l=?,mr=?"
        if(index(orb,"l=").eq.0.or.index(orb,"mr=").eq.0) then
          write(*,*)
          write(*,'("Error(getw90proj): No ang_mtm given in: ",A)') trim(wann_projstr(il))
          write(*,*)
          stop
        end if
        read(orb((index(orb,"l=")+2):index(orb,",")-1),'(I9)') l
        read(orb((index(orb,"mr=")+3):),'(I9)') mr
        do is=1,ncen
          wann_nproj = wann_nproj + 1
          cen(:,wann_nproj) = cen0(:,is)
          lmr(:,wann_nproj) = (/l,mr,rval/)
          zx(:,wann_nproj) = (/xval,zval/)
          zona(wann_nproj) = zonaval
        end do
      end if
    end if
    if(lastpass) exit
  end do
end do
! checks that the number of wannier functions given is equal to the number of projections
if(nspinor.eq.2) wann_nproj=2*wann_nproj
if(wann_nwf.ne.wann_nproj) then
  if(wann_nwf.gt.wann_nproj) then
    write(*,*)
    write(*,'("Warning(getw90proj): wann_nwf>wann_nproj - making up the &
          &difference by adding ",I3," random projections")') wann_nwf-wann_nproj
    do i=wann_nproj+1,wann_nwf
      lrand(i) = .true.
    end do
    wann_nproj = wann_nwf
  else
    write(*,*)
    write(*,'("Warning(getw90proj): wann_nwf<wann_nproj - setting wann_nwf=wann_nproj")')
    wann_nwf = wann_nproj
  end if
end if
if(nspinor.eq.2) wann_nproj=wann_nproj/2
! Actually allocated required arrays
if(allocated(wann_proj_site))     deallocate(wann_proj_site)
if(allocated(wann_proj_l))        deallocate(wann_proj_l)
if(allocated(wann_proj_m))        deallocate(wann_proj_m)
if(allocated(wann_proj_zaxis))    deallocate(wann_proj_zaxis)
if(allocated(wann_proj_xaxis))    deallocate(wann_proj_xaxis)
if(allocated(wann_proj_radial))   deallocate(wann_proj_radial)
if(allocated(wann_proj_zona))     deallocate(wann_proj_zona)
if(allocated(wann_proj_quantdir)) deallocate(wann_proj_quantdir)
if(allocated(wann_proj_spin))     deallocate(wann_proj_spin)
if(allocated(wann_proj_isrand))   deallocate(wann_proj_isrand)
allocate(wann_proj_site(3,wann_nproj))
allocate(wann_proj_l(wann_nproj))
allocate(wann_proj_m(wann_nproj))
allocate(wann_proj_zaxis(3,wann_nproj))
allocate(wann_proj_xaxis(3,wann_nproj))
allocate(wann_proj_radial(wann_nproj))
allocate(wann_proj_zona(wann_nproj))
allocate(wann_proj_isrand(wann_nproj))
wann_proj_isrand = .false.
do i=1,wann_nproj
  if(lrand(i)) then
    wann_proj_isrand(i) = .true.
    cycle
  end if
  !попробуй закомментить сначала это, прежде чем убирать остальной код
  wann_proj_site(:,i) = cen(:,i)
  wann_proj_l(i) = lmr(1,i)
  wann_proj_m(i) = lmr(2,i)
  wann_proj_radial(i) = lmr(3,i)
  wann_proj_xaxis(:,i) = zx(1:3,i)
  wann_proj_zaxis(:,i) = zx(4:6,i)
  wann_proj_zona(i) = zona(i)
end do

end subroutine getw90proj
!EOC
