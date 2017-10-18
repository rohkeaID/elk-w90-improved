
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: zfsht
! !INTERFACE:
subroutine zfsht(nr,nri,zfmt1,zfmt2)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   nr    : number of radial mesh points (in,integer)
!   nri   : number of points on the inner part of the muffin-tin (in,integer)
!   zfmt1 : input complex muffin-tin function in spherical coordinates
!           (in,complex(*))
!   zfmt2 : output complex muffin-tin function in spherical harmonics
!           (out,complex(*))
! !DESCRIPTION:
!   Performs a forward spherical harmonic transform (SHT) on a complex
!   muffin-tin function in spherical coordinates to obtain a function expressed
!   in spherical harmonics. See also {\tt genshtmat}.

! !REVISION HISTORY:
!   Created October 2013 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: nr,nri
complex(8), intent(in) :: zfmt1(*)
complex(8), intent(out) :: zfmt2(*)
! local variables
integer nro,i
! transform the inner part of the muffin-tin
call zgemm('N','N',lmmaxi,nri,lmmaxi,zone,zfshti,lmmaxi,zfmt1,lmmaxi,zzero, &
 zfmt2,lmmaxi)
! transform the outer part of the muffin-tin
nro=nr-nri
if (nro.eq.0) return
i=lmmaxi*nri+1
call zgemm('N','N',lmmaxo,nro,lmmaxo,zone,zfshto,lmmaxo,zfmt1(i),lmmaxo,zzero, &
 zfmt2(i),lmmaxo)
return
end subroutine
!EOC

