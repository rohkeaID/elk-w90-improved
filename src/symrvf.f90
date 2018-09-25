
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: symrvf
! !INTERFACE:
subroutine symrvf(tspin,tnc,nr,nri,np,ld,rvfmt,rvfir)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   tspin : .true. if spin rotations should be used (in,logical)
!   tnc   : .true. if the vector field is non-collinear, otherwise it is
!           collinear along the z-axis (in,logical)
!   nr    : number of radial points for each species (in,integer(nspecies))
!   nri   : number of radial points on the inner part (in,integer(nspecies))
!   np    : total number of points in each muffin-tin (in,integer(nspecies))
!   ld    : leading dimension (in,integer)
!   rvfmt : real muffin-tin vector field (in,real(ld,natmtot,*))
!   rvfir : real interstitial vector field (in,real(ngtot,*))
! !DESCRIPTION:
!   Symmetrises a vector field defined over the entire unit cell using the full
!   set of crystal symmetries. If a particular symmetry involves rotating atom
!   1 into atom 2, then the spatial and spin rotations of that symmetry are
!   applied to the vector field in atom 2 (expressed in spherical harmonic
!   coefficients), which is then added to the field in atom 1. This is repeated
!   for all symmetry operations. The fully symmetrised field in atom 1 is then
!   rotated and copied to atom 2. Symmetrisation of the interstitial part of the
!   field is performed by {\tt symrvfir}. See also {\tt symrfmt} and
!   {\tt findsym}.
!
! !REVISION HISTORY:
!   Created May 2007 (JKD)
!   Fixed problem with improper rotations, February 2008 (L. Nordstrom,
!    F. Bultmark and F. Cricchio)
!EOP
!BOC
implicit none
! arguments
logical, intent(in) :: tspin,tnc
integer, intent(in) :: nr(nspecies),nri(nspecies),np(nspecies)
integer, intent(in) :: ld
real(8), intent(inout) :: rvfmt(ld,natmtot,*),rvfir(ngtot,*)
! local variables
integer is,ia,ja,ias,jas
integer nd,isym,lspl,lspn,i
real(8) sc(3,3),v1(3),v2(3),t0,t1
! automatic arrays
logical done(natmmax)
! allocatable arrays
real(8), allocatable :: rvfmt1(:,:,:),rvfmt2(:,:)
! dimension of the vector field
if (tnc) then
  nd=3
else
  nd=1
end if
!-------------------------!
!     muffin-tin part     !
!-------------------------!
allocate(rvfmt1(npmtmax,natmmax,nd),rvfmt2(npmtmax,nd))
t0=1.d0/dble(nsymcrys)
do is=1,nspecies
! make copy of vector field for all atoms of current species
  do i=1,nd
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      call dcopy(np(is),rvfmt(:,ias,i),1,rvfmt1(:,ia,i),1)
    end do
  end do
  done(:)=.false.
  do ia=1,natoms(is)
    if (done(ia)) cycle
    ias=idxas(ia,is)
    rvfmt(1:np(is),ias,1:nd)=0.d0
! begin loop over crystal symmetries
    do isym=1,nsymcrys
! equivalent atom
      ja=ieqatom(ia,is,isym)
! parallel transport of vector field
      lspl=lsplsymc(isym)
      do i=1,nd
        call rotrfmt(symlatc(:,:,lspl),nr(is),nri(is),rvfmt1(:,ja,i), &
         rvfmt2(:,i))
      end do
      if (tspin) then
! global spin proper rotation matrix in Cartesian coordinates
        lspn=lspnsymc(isym)
        sc(:,:)=symlatd(lspn)*symlatc(:,:,lspn)
      else
! set spin rotation equal to spatial rotation
        lspn=lspl
        sc(:,:)=symlatc(:,:,lspl)
      end if
! global spin rotation of vector field
      if (tnc) then
! non-collinear case
        do i=1,np(is)
          v1(:)=rvfmt2(i,:)
          v2(1)=sc(1,1)*v1(1)+sc(1,2)*v1(2)+sc(1,3)*v1(3)
          v2(2)=sc(2,1)*v1(1)+sc(2,2)*v1(2)+sc(2,3)*v1(3)
          v2(3)=sc(3,1)*v1(1)+sc(3,2)*v1(2)+sc(3,3)*v1(3)
          rvfmt(i,ias,1:3)=rvfmt(i,ias,1:3)+v2(1:3)
        end do
      else
! collinear case
        t1=sc(3,3)
        do i=1,np(is)
          rvfmt(i,ias,1)=rvfmt(i,ias,1)+t1*rvfmt2(i,1)
        end do
      end if
! end loop over crystal symmetries
    end do
! normalise
    do i=1,nd
      call dscal(np(is),t0,rvfmt(:,ias,i),1)
    end do
! mark atom as done
    done(ia)=.true.
! rotate into equivalent atoms
    do isym=1,nsymcrys
      ja=ieqatom(ia,is,isym)
      if (done(ja)) cycle
      jas=idxas(ja,is)
! parallel transport of vector field (using operation inverse)
      lspl=isymlat(lsplsymc(isym))
      do i=1,nd
        call rotrfmt(symlatc(:,:,lspl),nr(is),nri(is),rvfmt(:,ias,i), &
         rvfmt(:,jas,i))
      end do
      if (tspin) then
! inverse of global proper rotation matrix in Cartesian coordinates
        lspn=isymlat(lspnsymc(isym))
        sc(:,:)=symlatd(lspn)*symlatc(:,:,lspn)
      else
! set spin rotation equal to spatial rotation
        lspn=lspl
        sc(:,:)=symlatc(:,:,lspl)
      end if
! global spin rotation of vector field
      if (tnc) then
! non-collinear case
        do i=1,np(is)
          v1(1:3)=rvfmt(i,jas,1:3)
          v2(1)=sc(1,1)*v1(1)+sc(1,2)*v1(2)+sc(1,3)*v1(3)
          v2(2)=sc(2,1)*v1(1)+sc(2,2)*v1(2)+sc(2,3)*v1(3)
          v2(3)=sc(3,1)*v1(1)+sc(3,2)*v1(2)+sc(3,3)*v1(3)
          rvfmt(i,jas,1:3)=v2(1:3)
        end do
      else
! collinear case
        t1=sc(3,3)
        do i=1,np(is)
          rvfmt(i,jas,1)=t1*rvfmt(i,jas,1)
        end do
      end if
! mark atom as done
      done(ja)=.true.
    end do
! end loop over atoms and species
  end do
end do
deallocate(rvfmt1,rvfmt2)
!---------------------------!
!     interstitial part     !
!---------------------------!
call symrvfir(tspin,tnc,rvfir)
return
end subroutine
!EOC

