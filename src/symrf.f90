
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: symrf
! !INTERFACE:
subroutine symrf(nr,nri,np,ld,rfmt,rfir)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   nr   : number of radial points for each species (in,integer(nspecies))
!   nri  : number of radial points on the inner part (in,integer(nspecies))
!   np   : total number of points in each muffin-tin (in,integer(nspecies))
!   ld   : leading dimension (in,integer)
!   rfmt : real muffin-tin function (inout,real(ld,natmtot))
!   rfir : real intersitial function (inout,real(ngtot))
! !DESCRIPTION:
!   Symmetrises a real scalar function defined over the entire unit cell using
!   the full set of crystal symmetries. In the muffin-tin of a particular atom
!   the spherical harmonic coefficients of every equivlent atom are rotated and
!   averaged. The interstitial part of the function is first Fourier transformed
!   to $G$-space, and then averaged over each symmetry by rotating the Fourier
!   coefficients and multiplying them by a phase factor corresponding to the
!   symmetry translation.
!
! !REVISION HISTORY:
!   Created May 2007 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: nr(nspecies),nri(nspecies),np(nspecies)
integer, intent(in) :: ld
real(8), intent(inout) :: rfmt(ld,natmtot),rfir(ngtot)
! local variables
integer is,ia,ja,ias,jas
integer isym,lspl
real(8) t0
! automatic arrays
logical done(natmmax)
! allocatable arrays
real(8), allocatable :: rfmt1(:,:),rfmt2(:)
!-------------------------!
!     muffin-tin part     !
!-------------------------!
allocate(rfmt1(ld,natmmax),rfmt2(ld))
t0=1.d0/dble(nsymcrys)
do is=1,nspecies
! make a copy of the input function
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    rfmt1(1:np(is),ia)=rfmt(1:np(is),ias)
  end do
  done(:)=.false.
! loop over atoms
  do ia=1,natoms(is)
    if (done(ia)) cycle
    ias=idxas(ia,is)
    rfmt(1:np(is),ias)=0.d0
! loop over crystal symmetries
    do isym=1,nsymcrys
! index to spatial rotation lattice symmetry
      lspl=lsplsymc(isym)
! equivalent atom index (symmetry rotates atom ja into atom ia)
      ja=ieqatom(ia,is,isym)
! apply the rotation to the muffin-tin function
      call rotrfmt(symlatc(:,:,lspl),nr(is),nri(is),rfmt1(:,ja),rfmt2)
! accumulate in original function array
      rfmt(1:np(is),ias)=rfmt(1:np(is),ias)+rfmt2(1:np(is))
    end do
! normalise
    rfmt(1:np(is),ias)=t0*rfmt(1:np(is),ias)
    done(ia)=.true.
! rotate into equivalent atoms
    do isym=1,nsymcrys
      ja=ieqatom(ia,is,isym)
      if (.not.done(ja)) then
        jas=idxas(ja,is)
! inverse symmetry (which rotates atom ia into atom ja)
        lspl=isymlat(lsplsymc(isym))
! rotate symmetrised function into equivalent muffin-tin
        call rotrfmt(symlatc(:,:,lspl),nr(is),nri(is),rfmt(:,ias),rfmt(:,jas))
        done(ja)=.true.
      end if
    end do
  end do
end do
deallocate(rfmt1,rfmt2)
!---------------------------!
!     interstitial part     !
!---------------------------!
call symrfir(ngvec,rfir)
return
end subroutine
!EOC

