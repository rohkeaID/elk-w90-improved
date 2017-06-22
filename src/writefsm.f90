
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writefsm(fnum)
use modmain
implicit none
! arguments
integer, intent(in) :: fnum
! local variables
integer is,ia,ias
if (fsmtype.eq.0) return
write(fnum,*)
if ((abs(fsmtype).eq.1).or.(abs(fsmtype).eq.3)) then
  write(fnum,'("FSM global effective field",T30,": ",3G18.10)') bfsmc(1:ndmag)
end if
if ((abs(fsmtype).eq.2).or.(abs(fsmtype).eq.3)) then
  write(fnum,'("FSM local muffin-tin effective fields :")')
  do is=1,nspecies
    write(fnum,'(" species : ",I4," (",A,")")') is,trim(spsymb(is))
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      write(fnum,'("  atom ",I4,T30,": ",3G18.10)') ia,bfsmcmt(1:ndmag,ias)
    end do
  end do
end if
return
end subroutine


