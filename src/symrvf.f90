
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: symrvf
! !INTERFACE:
subroutine symrvf(nr,nri,np,ld,rvfmt,rvfir)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   nr    : number of radial points for each species (in,integer(nspecies))
!   nri   : number of radial points on the inner part (in,integer(nspecies))
!   np    : total number of points in each muffin-tin (in,integer(nspecies))
!   ld    : leading dimension (in,integer)
!   rvfmt : real muffin-tin vector field (in,real(ld,natmtot,ndmag))
!   rvfir : real interstitial vector field (in,real(ngtot,ndmag))
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
integer, intent(in) :: nr(nspecies),nri(nspecies),np(nspecies)
integer, intent(in) :: ld
real(8), intent(inout) :: rvfmt(ld,natmtot,ndmag),rvfir(ngtot,ndmag)
! local variables
integer is,ia,ja,ias,jas
integer isym,lspl,lspn,md,i
real(8) sc(3,3),v(3),t0,t1
! automatic arrays
logical done(natmmax)
! allocatable arrays
real(8), allocatable :: rvfmt1(:,:,:),rvfmt2(:,:)
!-------------------------!
!     muffin-tin part     !
!-------------------------!
allocate(rvfmt1(npmtmax,natmmax,ndmag),rvfmt2(npmtmax,ndmag))
t0=1.d0/dble(nsymcrys)
do is=1,nspecies
! make copy of vector field for all atoms of current species
  do i=1,ndmag
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      rvfmt1(1:np(is),ia,i)=rvfmt(1:np(is),ias,i)
    end do
  end do
  done(:)=.false.
  do ia=1,natoms(is)
    if (done(ia)) cycle
    ias=idxas(ia,is)
    rvfmt(1:np(is),ias,1:ndmag)=0.d0
! begin loop over crystal symmetries
    do isym=1,nsymcrys
! equivalent atom
      ja=ieqatom(ia,is,isym)
! parallel transport of vector field
      lspl=lsplsymc(isym)
      do i=1,ndmag
        call rotrfmt(symlatc(:,:,lspl),nr(is),nri(is),rvfmt1(:,ja,i), &
         rvfmt2(:,i))
      end do
! global spin proper rotation matrix in Cartesian coordinates
      lspn=lspnsymc(isym)
      md=symlatd(lspn)
      sc(:,:)=dble(md)*symlatc(:,:,lspn)
! global spin rotation of vector field
      if (ncmag) then
! non-collinear case
        do i=1,np(is)
          v(:)=sc(:,1)*rvfmt2(i,1) &
              +sc(:,2)*rvfmt2(i,2) &
              +sc(:,3)*rvfmt2(i,3)
          rvfmt(i,ias,:)=rvfmt(i,ias,:)+v(:)
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
    rvfmt(1:np(is),ias,1:ndmag)=t0*rvfmt(1:np(is),ias,1:ndmag)
! mark atom as done
    done(ia)=.true.
! rotate into equivalent atoms
    do isym=1,nsymcrys
      ja=ieqatom(ia,is,isym)
      if (.not.done(ja)) then
        jas=idxas(ja,is)
! parallel transport of vector field (using operation inverse)
        lspl=isymlat(lsplsymc(isym))
        do i=1,ndmag
          call rotrfmt(symlatc(:,:,lspl),nr(is),nri(is),rvfmt(:,ias,i), &
           rvfmt(:,jas,i))
        end do
! inverse of global proper rotation matrix in Cartesian coordinates
        lspn=isymlat(lspnsymc(isym))
        md=symlatd(lspn)
        sc(:,:)=dble(md)*symlatc(:,:,lspn)
! global spin rotation of vector field
        if (ncmag) then
! non-collinear case
          do i=1,np(is)
            v(:)=sc(:,1)*rvfmt(i,jas,1) &
                +sc(:,2)*rvfmt(i,jas,2) &
                +sc(:,3)*rvfmt(i,jas,3)
            rvfmt(i,jas,:)=v(:)
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
      end if
    end do
! end loop over atoms and species
  end do
end do
deallocate(rvfmt1,rvfmt2)
!---------------------------!
!     interstitial part     !
!---------------------------!
call symrvfir(ngvec,rvfir)
return
end subroutine
!EOC

