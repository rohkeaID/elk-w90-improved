
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writedvs(fext)
use modmain
use modphonon
implicit none
! arguments
character(*), intent(in) :: fext
! local variables
integer is,ias
! allocatable arrays
complex(8), allocatable :: zfmt(:,:,:)
allocate(zfmt(lmmaxo,nrmtmax,natmtot))
open(50,file='DVS'//trim(fext),form='UNFORMATTED')
write(50) version
write(50) nspecies
write(50) lmmaxo
do is=1,nspecies
  write(50) natoms(is)
  write(50) nrmt(is)
end do
write(50) ngridg
do ias=1,natmtot
  is=idxis(ias)
  call zfmtpack(.false.,nrmt(is),nrmti(is),dvsmt(:,ias),zfmt(:,:,ias))
end do
write(50) zfmt,dvsir
close(50)
deallocate(zfmt)
return
end subroutine

