
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine getevecsv(fext,vpl,evecsv)
use modmain
implicit none
! arguments
character(*), intent(in) :: fext
real(8), intent(in) :: vpl(3)
complex(8), intent(out) :: evecsv(nstsv,nstsv)
! local variables
integer isym,lspn,ik,ist,i
integer recl,nstsv_
real(8) vkl_(3),det,v(3),th,t1
complex(8) su2(2,2),z1,z2
! find the k-point number
call findkpt(vpl,isym,ik)
! find the record length
inquire(iolength=recl) vkl_,nstsv_,evecsv
!$OMP CRITICAL
open(70,file=trim(scrpath)//'EVECSV'//trim(fext),action='READ', &
 form='UNFORMATTED',access='DIRECT',recl=recl)
read(70,rec=ik) vkl_,nstsv_,evecsv
close(70)
!$OMP END CRITICAL
t1=abs(vkl(1,ik)-vkl_(1))+abs(vkl(2,ik)-vkl_(2))+abs(vkl(3,ik)-vkl_(3))
if (t1.gt.epslat) then
  write(*,*)
  write(*,'("Error(getevecsv): differing vectors for k-point ",I8)') ik
  write(*,'(" current    : ",3G18.10)') vkl(:,ik)
  write(*,'(" EVECSV.OUT : ",3G18.10)') vkl_
  write(*,*)
  stop
end if
if (nstsv.ne.nstsv_) then
  write(*,*)
  write(*,'("Error(getevecsv): differing nstsv for k-point ",I8)') ik
  write(*,'(" current    : ",I8)') nstsv
  write(*,'(" EVECSV.OUT : ",I8)') nstsv_
  write(*,*)
  stop
end if
! if eigenvectors are spin-unpolarised return
if (.not.spinpol) return
! index to global spin rotation in lattice point group
lspn=lspnsymc(isym)
! if symmetry element is the identity return
if (lspn.eq.1) return
! find the SU(2) representation of the spin rotation matrix
call rotaxang(epslat,symlatc(:,:,lspn),det,v,th)
call axangsu2(v,th,su2)
! apply SU(2) matrix to second-variational states (active transformation)
do i=1,nstsv
  do ist=1,nstfv
    z1=evecsv(ist,i)
    z2=evecsv(ist+nstfv,i)
    evecsv(ist,i)=su2(1,1)*z1+su2(1,2)*z2
    evecsv(ist+nstfv,i)=su2(2,1)*z1+su2(2,2)*z2
  end do
end do
return
end subroutine

