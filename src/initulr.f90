
! Copyright (C) 2018 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine initulr
use modmain
use modulr
implicit none
! local variables
integer iq

! allocate long-range density and magnetisation arrays
if (allocated(rhormt)) deallocate(rhormt)
allocate(rhormt(npcmtmax,natmtot,nqpt))
if (allocated(rhorir)) deallocate(rhorir)
allocate(rhorir(ngtot,nqpt))
if (allocated(magrmt)) deallocate(magrmt)
if (allocated(magrir)) deallocate(magrir)
if (spinpol) then
  allocate(magrmt(npcmtmax,natmtot,ndmag,nqpt))
  allocate(magrir(ngtot,ndmag,nqpt))
end if
if (allocated(rhoqmt)) deallocate(rhoqmt)
allocate(rhoqmt(npcmtmax,natmtot,nqpt))
if (allocated(rhoqir)) deallocate(rhoqir)
allocate(rhoqir(ngtot,nqpt))
if (allocated(magqmt)) deallocate(magqmt)
if (allocated(magqir)) deallocate(magqir)
if (spinpol) then
  allocate(magqmt(npcmtmax,natmtot,ndmag,nqpt))
  allocate(magqir(ngtot,ndmag,nqpt))
end if

! G+Q-vector arrays
if (allocated(gqc)) deallocate(gqc)
allocate(gqc(ngvec,nqpt))
do iq=1,nqpt

end do

return
end subroutine
