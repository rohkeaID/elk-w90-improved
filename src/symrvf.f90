
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: symrvf
! !INTERFACE:
subroutine symrvf(lrstp,rvfmt,rvfir)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   lrstp : radial step length (in,integer)
!   rvfmt : real muffin-tin vector field
!           (in,real(lmmaxvr,nrmtmax,natmtot,ndmag))
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
integer, intent(in) :: lrstp
real(8), intent(inout) :: rvfmt(lmmaxvr,nrmtmax,natmtot,ndmag)
real(8), intent(inout) :: rvfir(ngtot,ndmag)
! local variables
integer is,ia,ja,ias,jas
integer nr,nri,ir,lmmax,lm
integer isym,lspl,ilspl
integer lspn,ilspn,md,i
real(8) sc(3,3),v(3),t0,t1
! automatic arrays
logical done(natmmax)
! allocatable arrays
real(8), allocatable :: rvfmt1(:,:,:,:),rvfmt2(:,:,:)
!-------------------------!
!     muffin-tin part     !
!-------------------------!
allocate(rvfmt1(lmmaxvr,nrmtmax,natmmax,ndmag))
allocate(rvfmt2(lmmaxvr,nrmtmax,ndmag))
t0=1.d0/dble(nsymcrys)
do is=1,nspecies
  nr=nrmt(is)
  nri=nrmtinr(is)
! make copy of vector field for all atoms of current species
  do i=1,ndmag
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      call rfmtcopy(nr,nri,lrstp,rvfmt(:,:,ias,i),rvfmt1(:,:,ia,i))
    end do
  end do
  done(:)=.false.
  do ia=1,natoms(is)
    if (done(ia)) cycle
    ias=idxas(ia,is)
    do i=1,ndmag
      call rfmtzero(nr,nri,lrstp,rvfmt(:,:,ias,i))
    end do
! begin loop over crystal symmetries
    do isym=1,nsymcrys
! equivalent atom
      ja=ieqatom(ia,is,isym)
! parallel transport of vector field
      lspl=lsplsymc(isym)
      do i=1,ndmag
        call rotrfmt(symlatc(:,:,lspl),nr,nri,lrstp,rvfmt1(:,:,ja,i), &
         rvfmt2(:,:,i))
      end do
! global spin proper rotation matrix in Cartesian coordinates
      lspn=lspnsymc(isym)
      md=symlatd(lspn)
      sc(:,:)=dble(md)*symlatc(:,:,lspn)
! global spin rotation of vector field
      if (ncmag) then
! non-collinear case
        lmmax=lmmaxinr
        do ir=1,nr,lrstp
          do lm=1,lmmax
            v(:)=sc(:,1)*rvfmt2(lm,ir,1) &
                +sc(:,2)*rvfmt2(lm,ir,2) &
                +sc(:,3)*rvfmt2(lm,ir,3)
            rvfmt(lm,ir,ias,:)=rvfmt(lm,ir,ias,:)+v(:)
          end do
          if (ir.eq.nri) lmmax=lmmaxvr
        end do
      else
! collinear case
        t1=sc(3,3)
        lmmax=lmmaxinr
        do ir=1,nr,lrstp
          rvfmt(1:lmmax,ir,ias,1)=rvfmt(1:lmmax,ir,ias,1) &
           +t1*rvfmt2(1:lmmax,ir,1)
          if (ir.eq.nri) lmmax=lmmaxvr
        end do
      end if
! end loop over crystal symmetries
    end do
! normalise
    do i=1,ndmag
      call rfmtscal(nr,nri,lrstp,t0,rvfmt(:,:,ias,i))
    end do
! mark atom as done
    done(ia)=.true.
! rotate into equivalent atoms
    do isym=1,nsymcrys
      ja=ieqatom(ia,is,isym)
      if (.not.done(ja)) then
        jas=idxas(ja,is)
! parallel transport of vector field (using operation inverse)
        lspl=lsplsymc(isym)
        ilspl=isymlat(lspl)
        do i=1,ndmag
          call rotrfmt(symlatc(:,:,ilspl),nr,nri,lrstp,rvfmt(:,:,ias,i), &
           rvfmt(:,:,jas,i))
        end do
! inverse of global proper rotation matrix in Cartesian coordinates
        lspn=lspnsymc(isym)
        ilspn=isymlat(lspn)
        md=symlatd(ilspn)
        sc(:,:)=dble(md)*symlatc(:,:,ilspn)
! global spin rotation of vector field
        if (ncmag) then
! non-collinear case
          lmmax=lmmaxinr
          do ir=1,nr,lrstp
            do lm=1,lmmax
              v(:)=rvfmt(lm,ir,jas,:)
              rvfmt(lm,ir,jas,:)=sc(:,1)*v(1)+sc(:,2)*v(2)+sc(:,3)*v(3)
            end do
            if (ir.eq.nri) lmmax=lmmaxvr
          end do
        else
! collinear case
          t1=sc(3,3)
          lmmax=lmmaxinr
          do ir=1,nr,lrstp
            rvfmt(1:lmmax,ir,jas,1)=t1*rvfmt(1:lmmax,ir,jas,1)
            if (ir.eq.nri) lmmax=lmmaxvr
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

