
! Copyright (C) 2017 A. Davydov, A. Sanna, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine putgwsefm(ik,se)
use modmain
use modgw
implicit none
! arguments
integer, intent(in) :: ik
complex(8), intent(in) :: se(nstsv,nstsv,0:nwfm)
! local variables
integer recl
! find the record length
inquire(iolength=recl) vkl(:,ik),nstsv,nwfm,se
!$OMP CRITICAL(u280)
open(280,file='GWSEFM.OUT',form='UNFORMATTED',access='DIRECT',recl=recl)
write(280,rec=ik) vkl(:,ik),nstsv,nwfm,se
close(280)
!$OMP END CRITICAL(u280)
return
end subroutine

