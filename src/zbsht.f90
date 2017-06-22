
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: zbsht
! !INTERFACE:
subroutine zbsht(nr,nri,zfmt1,zfmt2)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   nr    : number of radial mesh points (in,integer)
!   nri   : number of points on the inner part of the muffin-tin (in,integer)
!   zfmt1 : input complex muffin-tin function in spherical harmonics
!           (in,complex(lmmaxvr,nr))
!   zfmt2 : output complex muffin-tin function in spherical coordinates
!           (out,complex(lmmaxvr,nr))
! !DESCRIPTION:
!   Performs a backward spherical harmonic transform (SHT) on a complex
!   muffin-tin function expressed in spherical harmonics to obtain a function in
!   spherical coordinates. See also {\tt genshtmat}.
!
! !REVISION HISTORY:
!   Created October 2013 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: nr,nri
complex(8), intent(in) :: zfmt1(lmmaxvr,nr)
complex(8), intent(out) :: zfmt2(lmmaxvr,nr)
! local variables
integer nro,iro
! transform the inner part of the muffin-tin
call zgemm('N','N',lmmaxinr,nri,lmmaxinr,zone,zbshtinr,lmmaxinr,zfmt1,lmmaxvr, &
 zzero,zfmt2,lmmaxvr)
! transform the outer part of the muffin-tin
nro=nr-nri
if (nro.eq.0) return
iro=nri+1
call zgemm('N','N',lmmaxvr,nro,lmmaxvr,zone,zbshtvr,lmmaxvr,zfmt1(:,iro), &
 lmmaxvr,zzero,zfmt2(:,iro),lmmaxvr)
return
end subroutine
!EOC

