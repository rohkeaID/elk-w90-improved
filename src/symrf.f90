
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: symrf
! !INTERFACE:
subroutine symrf(lrstp,rfmt,rfir)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   lrstp : radial step length (in,integer)
!   rfmt  : real muffin-tin function (inout,real(lmmaxvr,nrmtmax,natmtot))
!   rfir  : real intersitial function (inout,real(ngtot))
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
integer, intent(in) :: lrstp
real(8), intent(inout) :: rfmt(lmmaxvr,nrmtmax,natmtot),rfir(ngtot)
! local variables
integer is,ia,ja,ias,jas
integer nr,nri
integer isym,lspl,ilspl
real(8) t1
! automatic arrays
logical done(natmmax)
! allocatable arrays
real(8), allocatable :: rfmt1(:,:,:),rfmt2(:,:)
!-------------------------!
!     muffin-tin part     !
!-------------------------!
allocate(rfmt1(lmmaxvr,nrmtmax,natmmax))
allocate(rfmt2(lmmaxvr,nrmtmax))
t1=1.d0/dble(nsymcrys)
do is=1,nspecies
  nr=nrmt(is)
  nri=nrmtinr(is)
! make a copy of the input function
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    call rfmtcopy(nr,nri,lrstp,rfmt(:,:,ias),rfmt1(:,:,ia))
  end do
  done(:)=.false.
! loop over atoms
  do ia=1,natoms(is)
    if (done(ia)) cycle
    ias=idxas(ia,is)
    call rfmtzero(nr,nri,lrstp,rfmt(:,:,ias))
! loop over crystal symmetries
    do isym=1,nsymcrys
! index to spatial rotation lattice symmetry
      lspl=lsplsymc(isym)
! equivalent atom index (symmetry rotates atom ja into atom ia)
      ja=ieqatom(ia,is,isym)
! apply the rotation to the muffin-tin function
      call rotrfmt(symlatc(:,:,lspl),nr,nri,lrstp,rfmt1(:,:,ja),rfmt2)
! accumulate in original function array
      call rfmtadd(nr,nri,lrstp,rfmt2,rfmt(:,:,ias))
    end do
! normalise
    call rfmtscal(nr,nri,lrstp,t1,rfmt(:,:,ias))
    done(ia)=.true.
! rotate into equivalent atoms
    do isym=1,nsymcrys
      ja=ieqatom(ia,is,isym)
      if (.not.done(ja)) then
        jas=idxas(ja,is)
        lspl=lsplsymc(isym)
! inverse symmetry (which rotates atom ia into atom ja)
        ilspl=isymlat(lspl)
! rotate symmetrised function into equivalent muffin-tin
        call rotrfmt(symlatc(:,:,ilspl),nr,nri,lrstp,rfmt(:,:,ias), &
         rfmt(:,:,jas))
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

