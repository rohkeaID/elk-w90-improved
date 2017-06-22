
! Copyright (C) 2007 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine symdmat(lmax,ld,dmat)
use modmain
implicit none
! arguments
integer, intent(in) :: lmax
integer, intent(in) :: ld
complex(8), intent(inout) :: dmat(ld,nspinor,ld,nspinor,natmtot)
! local variables
integer isym,lspl,lspn,lmmax
integer is,ia,ja,ias,jas
real(8) t1
! automatic arrays
logical done(natmmax)
! allocatable arrays
complex(8), allocatable :: dm(:,:,:,:,:)
lmmax=(lmax+1)**2
! allocate local arrays
allocate(dm(ld,nspinor,ld,nspinor,natmmax))
t1=1.d0/dble(nsymcrys)
do is=1,nspecies
! make copy of the density matrices
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    dm(1:lmmax,:,1:lmmax,:,ia)=dmat(1:lmmax,:,1:lmmax,:,ias)
  end do
  done(:)=.false.
  do ia=1,natoms(is)
    if (done(ia)) cycle
    ias=idxas(ia,is)
    dmat(:,:,:,:,ias)=0.d0
    do isym=1,nsymcrys
      lspl=lsplsymc(isym)
      lspn=lspnsymc(isym)
! equivalent atom index (symmetry rotates atom ja into atom ia)
      ja=ieqatom(ia,is,isym)
      call rotdmat(symlatc(:,:,lspl),symlatc(:,:,lspn),lmax,nspinor,ld, &
       dm(:,:,:,:,ja),dmat(:,:,:,:,ias))
! end loop over crystal symmetries
    end do
! normalise
    dmat(:,:,:,:,ias)=t1*dmat(:,:,:,:,ias)
    done(ia)=.true.
! rotate into equivalent atoms
    do isym=1,nsymcrys
      ja=ieqatom(ia,is,isym)
      if (.not.done(ja)) then
        jas=idxas(ja,is)
        lspl=lsplsymc(isym)
        lspn=lspnsymc(isym)
        dmat(:,:,:,:,jas)=0.d0
        call rotdmat(symlatc(:,:,lspl),symlatc(:,:,lspn),lmax,nspinor,ld, &
         dmat(:,:,:,:,ias),dmat(:,:,:,:,jas))
        done(ja)=.true.
      end if
    end do
! end loop over atoms and species
  end do
end do
deallocate(dm)
return
end subroutine

