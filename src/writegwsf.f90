
! Copyright (C) 2016 A. Davydov, A. Sanna, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine writegwsf(ik,sf)
use modmain
use modgw
implicit none
! arguments
integer, intent(in) :: ik
real(8), intent(in) :: sf(nwplot)
! local variables
integer iw,ist
real(8) dw,w,e
character(256) fname
write(fname,'("GWSF_K",I6.6,".OUT")') ik
open(50,file=trim(fname),form='FORMATTED')
! write the GW spectral function
dw=(wplot(2)-wplot(1))/dble(nwplot)
do iw=1,nwplot
  w=dw*dble(iw-1)+wplot(1)
  write(50,'(2G18.10)') w,sf(iw)
end do
write(50,'("     ")')
! write the Kohn-Sham eigenvalues for reference
do ist=1,nstsv
  e=evalsv(ist,ik)-efermi
  write(50,'(2G18.10)') e,0.d0
  write(50,'(2G18.10)') e,1.d0/swidth
  write(50,'("     ")')
end do
close(50)
return
end subroutine

