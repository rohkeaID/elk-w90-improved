
! Copyright (C) 2005-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dynsymapp(isym,vpl,dyn,dyns)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: isym
real(8), intent(in) :: vpl(3)
complex(8), intent(in) :: dyn(nbph,nbph)
complex(8), intent(inout) :: dyns(nbph,nbph)
! local variables
integer is,ia,ja,ias,jas
integer lspl,ilspl,i,j,k,l,m,n
real(8) sl(3,3),sic(3,3),v(3),t1
real(8) a(3,3),b(3,3),c(3,3)
complex(8) z1
! automatic arrays
integer map(natmtot)
complex(8) zph(natmtot)
! index to spatial rotation in lattice point group
lspl=lsplsymc(isym)
! the inverse of the spatial symmetry
ilspl=isymlat(lspl)
! symmetry matrix in lattice coordinates
sl(:,:)=dble(symlat(:,:,lspl))
! inverse symmetry matrix in Cartesian coordinates
sic(:,:)=symlatc(:,:,ilspl)
! operate with symmetry matrix on vpl
call r3mtv(sl,vpl,v)
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! equivalent atom with this symmetry
    ja=ieqatom(ia,is,isym)
    jas=idxas(ja,is)
    map(ias)=jas
! phase factor
    t1=twopi*dot_product(vpl(:),atposl(:,ia,is))
    z1=cmplx(cos(t1),sin(t1),8)
    zph(ias)=z1
    t1=-twopi*dot_product(v(:),atposl(:,ja,is))
    zph(ias)=z1*cmplx(cos(t1),sin(t1),8)
  end do
end do
! rotate and phase-shift dynamical matrix with symmetry
do ias=1,natmtot
  i=3*(ias-1)
  k=3*(map(ias)-1)
  do jas=1,natmtot
    j=3*(jas-1)
    l=3*(map(jas)-1)
    do m=1,3
      do n=1,3
        a(m,n)=dble(dyn(i+m,j+n))
        b(m,n)=aimag(dyn(i+m,j+n))
      end do
    end do
    call r3mm(sic,a,c)
    call r3mmt(c,sic,a)
    call r3mm(sic,b,c)
    call r3mmt(c,sic,b)
    z1=zph(ias)*conjg(zph(jas))
    do m=1,3
      do n=1,3
        dyns(k+m,l+n)=dyns(k+m,l+n)+z1*cmplx(a(m,n),b(m,n),8)
      end do
    end do
  end do
end do
return
end subroutine

