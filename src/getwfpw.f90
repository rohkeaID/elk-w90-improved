
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine getwfpw(vpl,vhpl,wfpw)
use modmain
use modpw
implicit none
! arguments
real(8), intent(in) :: vpl(3)
real(8), intent(in) :: vhpl(3,nhkmax,nspnfv)
complex(8), intent(out) :: wfpw(nhkmax,nspinor,nstsv)
! local variables
integer isym,lspl,ilspl,lspn
integer ik,ist,ihk,ihp,jhp,ig
integer ispn0,ispn1,jspn,i
integer recl,nhkmax_,nspinor_,nstsv_
real(8) vkl_(3),si(3,3)
real(8) v(3),det,th,t1
complex(8) su2(2,2),z1,z2
! allocatable arrays
complex(8), allocatable :: wfpw_(:,:,:)
! find the equivalent k-point number and crystal symmetry element
call findkpt(vpl,isym,ik)
! index to spatial rotation in lattice point group
lspl=lsplsymc(isym)
! find the record length
inquire(iolength=recl) vkl_,nhkmax_,nspinor_,nstsv_,wfpw
!$OMP CRITICAL
open(105,file='WFPW.OUT',action='READ',form='UNFORMATTED',access='DIRECT', &
 recl=recl)
read(105,rec=ik) vkl_,nhkmax_,nspinor_,nstsv_,wfpw
close(105)
!$OMP END CRITICAL
t1=abs(vkl(1,ik)-vkl_(1))+abs(vkl(2,ik)-vkl_(2))+abs(vkl(3,ik)-vkl_(3))
if (t1.gt.epslat) then
  write(*,*)
  write(*,'("Error(getwfpw): differing vectors for k-point ",I8)') ik
  write(*,'(" current  : ",3G18.10)') vkl(:,ik)
  write(*,'(" WFPW.OUT : ",3G18.10)') vkl_
  write(*,*)
  stop
end if
if (nhkmax.ne.nhkmax_) then
  write(*,*)
  write(*,'("Error(getwfpw): differing nhkmax for k-point ",I8)') ik
  write(*,'(" current  : ",I8)') nhkmax
  write(*,'(" WFPW.OUT : ",I8)') nhkmax_
  write(*,*)
  stop
end if
if (nspinor.ne.nspinor_) then
  write(*,*)
  write(*,'("Error(getwfpw): differing nspinor for k-point ",I8)') ik
  write(*,'(" current  : ",I8)') nspinor
  write(*,'(" WFPW.OUT : ",I8)') nspinor_
  write(*,*)
  stop
end if
if (nstsv.ne.nstsv_) then
  write(*,*)
  write(*,'("Error(getwfpw): differing nstsv for k-point ",I8)') ik
  write(*,'(" current  : ",I8)') nstsv
  write(*,'(" WFPW.OUT : ",I8)') nstsv_
  write(*,*)
  stop
end if
! if p = k then return
t1=abs(vpl(1)-vkl(1,ik))+abs(vpl(2)-vkl(2,ik))+abs(vpl(3)-vkl(3,ik))
if (t1.lt.epslat) return
!--------------------------------------------------------!
!     translate and rotate wavefunction coefficients     !
!--------------------------------------------------------!
! allocate temporary copy of wavefunction
allocate(wfpw_(nhkmax,nspinor,nstsv))
! the inverse of the spatial symmetry rotates k into p
ilspl=isymlat(lspl)
si(:,:)=dble(symlat(:,:,ilspl))
! loop over first-variational spins
do jspn=1,nspnfv
  if (spinsprl) then
    ispn0=jspn; ispn1=jspn
  else
    ispn0=1; ispn1=nspinor
  end if
! apply translation operation if required
  if (tvzsymc(isym)) then
! translation vector is zero
    do ihk=1,nhk(jspn,ik)
      wfpw_(ihk,ispn0:ispn1,:)=wfpw(ihk,ispn0:ispn1,:)
    end do
  else
! non-zero translation vector gives a phase factor
    do ihk=1,nhk(jspn,ik)
      ig=ihkig(ihk,jspn,ik)
      t1=-twopi*dot_product(dble(ivg(:,ig)),vtlsymc(:,isym))
      z1=cmplx(cos(t1),sin(t1),8)
      wfpw_(ihk,ispn0:ispn1,:)=z1*wfpw(ihk,ispn0:ispn1,:)
    end do
  end if
! apply spatial rotation operation (passive transformation)
  i=1
  do ihk=1,nhk(jspn,ik)
    call r3mtv(si,vhkl(:,ihk,jspn,ik),v)
    do ihp=i,nhk(jspn,ik)
      t1=abs(v(1)-vhpl(1,ihp,jspn)) &
        +abs(v(2)-vhpl(2,ihp,jspn)) &
        +abs(v(3)-vhpl(3,ihp,jspn))
      if (t1.lt.epslat) then
        wfpw(ihp,ispn0:ispn1,:)=wfpw_(ihk,ispn0:ispn1,:)
        if (ihp.eq.i) i=i+1
        exit
      end if
    end do
  end do
end do
! apply spin rotation if required
if (spinpol) then
! index to global spin rotation in lattice point group
  lspn=lspnsymc(isym)
! if symmetry element is the identity return
  if (lspn.eq.1) return
! find the SU(2) representation of the spin rotation matrix
  call rotaxang(epslat,symlatc(:,:,lspn),det,v,th)
  call axangsu2(v,th,su2)
! apply SU(2) matrix to spinor wavefunctions (active transformation)
  if (spinsprl) then
! spin-spiral case
    wfpw(:,2,:)=0.d0
    i=1
    do ihp=1,nhk(1,ik)
      v(:)=vhpl(:,ihp,1)-vqlss(:)
      do jhp=i,nhk(2,ik)
        t1=abs(v(1)-vhpl(1,jhp,2)) &
          +abs(v(2)-vhpl(2,jhp,2)) &
          +abs(v(3)-vhpl(3,jhp,2))
        if (t1.lt.epslat) then
          do ist=1,nstsv
            z1=wfpw(ihp,1,ist)
            z2=wfpw(jhp,2,ist)
            wfpw(ihp,1,ist)=su2(1,1)*z1+su2(1,2)*z2
            wfpw(jhp,2,ist)=su2(2,1)*z1+su2(2,2)*z2
          end do
          if (jhp.eq.i) i=i+1
          goto 10
        end if
      end do
      wfpw(ihp,1,:)=0.d0
10 continue
    end do
  else
! normal spin case
    do ist=1,nstsv
      do ihp=1,nhk(1,ik)
        z1=wfpw(ihp,1,ist)
        z2=wfpw(ihp,2,ist)
        wfpw(ihp,1,ist)=su2(1,1)*z1+su2(1,2)*z2
        wfpw(ihp,2,ist)=su2(2,1)*z1+su2(2,2)*z2
      end do
    end do
  end if
end if
deallocate(wfpw_)
return
end subroutine

