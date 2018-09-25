
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine getpmat(tfv,vpl,pmat)
use modmain
implicit none
! arguments
logical, intent(in) :: tfv
real(8), intent(in) :: vpl(3)
complex(8), intent(out) :: pmat(nstsv,nstsv,3)
! local variables
integer isym,ik,ist,jst
integer recl,nstsv_,i
real(8) vkl_(3),sc(3,3),t1
real(8) v1(3),v2(3),v3(3)
! find the k-point number
call findkpt(vpl,isym,ik)
! find the record length
inquire(iolength=recl) vkl_,nstsv_,pmat
if (tfv) then
!$OMP CRITICAL(u152)
  do i=1,2
    open(152,file='PMATFV.OUT',form='UNFORMATTED',access='DIRECT',recl=recl, &
     err=10)
    read(152,rec=ik,err=10) vkl_,nstsv_,pmat
    exit
10 continue
    if (i.eq.2) then
      write(*,*)
      write(*,'("Error(getpmat): unable to read from PMATFV.OUT")')
      write(*,*)
      stop
    end if
    close(152)
  end do
!$OMP END CRITICAL(u152)
else
!$OMP CRITICAL(u150)
  do i=1,2
    open(150,file='PMAT.OUT',form='UNFORMATTED',access='DIRECT',recl=recl, &
     err=20)
    read(150,rec=ik,err=20) vkl_,nstsv_,pmat
    exit
20 continue
    if (i.eq.2) then
      write(*,*)
      write(*,'("Error(getpmat): unable to read from PMAT.OUT")')
      write(*,*)
      stop
    end if
    close(150)
  end do
!$OMP END CRITICAL(u150)
end if
t1=abs(vkl(1,ik)-vkl_(1))+abs(vkl(2,ik)-vkl_(2))+abs(vkl(3,ik)-vkl_(3))
if (t1.gt.epslat) then
  write(*,*)
  write(*,'("Error(getpmat): differing vectors for k-point ",I8)') ik
  write(*,'(" current  : ",3G18.10)') vkl(:,ik)
  write(*,'(" PMAT.OUT : ",3G18.10)') vkl_
  write(*,*)
  stop
end if
if (nstsv.ne.nstsv_) then
  write(*,*)
  write(*,'("Error(getpmat): differing nstsv for k-point ",I8)') ik
  write(*,'(" current  : ",I8)') nstsv
  write(*,'(" PMAT.OUT : ",I8)') nstsv_
  write(*,*)
  stop
end if
! if p = k then return
t1=abs(vpl(1)-vkl(1,ik))+abs(vpl(2)-vkl(2,ik))+abs(vpl(3)-vkl(3,ik))
if (t1.lt.epslat) return
! rotate the matrix elements from the reduced to non-reduced k-point
sc(:,:)=symlatc(:,:,lsplsymc(isym))
do ist=1,nstsv
  do jst=1,nstsv
    v1(:)=dble(pmat(ist,jst,:))
    call r3mv(sc,v1,v2)
    v1(:)=aimag(pmat(ist,jst,:))
    call r3mv(sc,v1,v3)
    pmat(ist,jst,:)=cmplx(v2(:),v3(:),8)
  end do
end do
return
end subroutine

