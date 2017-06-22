
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine phdelete
use modmain
use modphonon
implicit none
character(256) fext
! construct the phonon file extension
call phfext(iqph,isph,iaph,ipph,fext)
! delete the eigenvalue/vector files
open(70,file='DEVALFV'//trim(fext))
close(70,status='DELETE')
open(70,file='DEVECFV'//trim(fext))
close(70,status='DELETE')
open(70,file='DEVECSV'//trim(fext))
close(70,status='DELETE')
return
end subroutine

