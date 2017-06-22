
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine phscdelete
use modmain
implicit none
! delete the eigenvector files
call delevec
! delete the eigenvalue files
open(70,file='EVALFV'//trim(filext))
close(70,status='DELETE')
open(70,file='EVALSV'//trim(filext))
close(70,status='DELETE')
! delete the occupancy file
open(70,file='OCCSV'//trim(filext))
close(70,status='DELETE')
return
end subroutine

