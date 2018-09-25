
! Copyright (C) 2010 Alexey I. Baranov.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genhvec
use modmain
use modpw
implicit none
! local variables
logical lsym(48)
integer ih,jh,kh,lh,k
integer i1,i2,i3,iv(3)
integer nsym,isym,sym(3,3,48)
real(8) v1(3),v2(3),v3(3)
! allocatable arrays
integer, allocatable :: idx(:),ivh0(:,:)
real(8), allocatable :: vhc0(:,:),hc0(:)
! find the H-vector grid sizes
call gridsize(avec,hmaxvr,ngridh,nhtot,inthv)
! allocate global H-vector arrays
if (allocated(ivh)) deallocate(ivh)
allocate(ivh(3,nhtot))
if (allocated(mulh)) deallocate(mulh)
allocate(mulh(nhtot))
if (allocated(vhc)) deallocate(vhc)
allocate(vhc(3,nhtot))
if (allocated(hc)) deallocate(hc)
allocate(hc(nhtot))
! allocate local arrays
allocate(idx(nhtot),ivh0(3,nhtot))
allocate(vhc0(3,nhtot),hc0(nhtot))
ih=0
do i1=inthv(1,1),inthv(2,1)
  v1(:)=dble(i1)*bvec(:,1)
  do i2=inthv(1,2),inthv(2,2)
    v2(:)=v1(:)+dble(i2)*bvec(:,2)
    do i3=inthv(1,3),inthv(2,3)
      v3(:)=v2(:)+dble(i3)*bvec(:,3)
      ih=ih+1
! map from H-vector to (i1,i2,i3) index
      ivh0(1,ih)=i1
      ivh0(2,ih)=i2
      ivh0(3,ih)=i3
! H-vector in Cartesian coordinates
      vhc0(:,ih)=v3(:)
! length of each H-vector
      hc0(ih)=sqrt(v3(1)**2+v3(2)**2+v3(3)**2)
    end do
  end do
end do
! sort by vector length
call sortidx(nhtot,hc0,idx)
! reorder arrays
do ih=1,nhtot
  jh=idx(ih)
  ivh(:,ih)=ivh0(:,jh)
  hc(ih)=hc0(jh)
  vhc(:,ih)=vhc0(:,jh)
end do
! find the number of vectors with H < hmaxvr
nhvec=1
do ih=nhtot,1,-1
  if (hc(ih).lt.hmaxvr) then
    nhvec=ih
    exit
  end if
end do
! find the subgroup of symmorphic, non-magnetic symmetries
lsym(:)=.false.
do isym=1,nsymcrys
  if (tv0symc(isym).and.(lspnsymc(isym).eq.1)) lsym(lsplsymc(isym))=.true.
end do
nsym=0
do isym=1,nsymlat
  if (lsym(isym)) then
    nsym=nsym+1
    sym(:,:,nsym)=symlat(:,:,isym)
  end if
end do
if (reduceh) then
! find the subgroup of symmorphic, non-magnetic symmetries
  lsym(:)=.false.
  do isym=1,nsymcrys
    if (tv0symc(isym).and.(lspnsymc(isym).eq.1)) lsym(lsplsymc(isym))=.true.
  end do
  nsym=0
  do isym=1,nsymlat
    if (lsym(isym)) then
      nsym=nsym+1
      sym(:,:,nsym)=symlat(:,:,isym)
    end if
  end do
else
! use only the identity element if no reduction is required
  nsym=1
end if
! reduce the H-vector set with the symmetries if required
if (nsym.gt.1) then
  ivh0(:,1:nhvec)=ivh(:,1:nhvec)
  hc0(1:nhvec)=hc(1:nhvec)
  vhc0(:,1:nhvec)=vhc(:,1:nhvec)
  kh=0
  lh=nhvec
  do ih=1,nhvec
    do isym=1,nsym
      call i3mtv(sym(:,:,isym),ivh0(:,ih),iv(:))
      do jh=1,kh
        k=abs(ivh(1,jh)-iv(1))+abs(ivh(2,jh)-iv(2))+abs(ivh(3,jh)-iv(3))
        if (k.eq.0) then
          ivh(:,lh)=ivh0(:,ih)
          hc(lh)=hc0(ih)
          vhc(:,lh)=vhc0(:,ih)
          lh=lh-1
          mulh(jh)=mulh(jh)+1
          goto 10
        end if
      end do
    end do
    kh=kh+1
    ivh(:,kh)=ivh0(:,ih)
    hc(kh)=hc0(ih)
    vhc(:,kh)=vhc0(:,ih)
    mulh(kh)=1
10 continue
  end do
  nhvec=kh
else
  mulh(:)=1
end if
deallocate(idx,ivh0,vhc0,hc0)
return
end subroutine

