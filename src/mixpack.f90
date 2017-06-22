
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: mixpack
! !INTERFACE:
subroutine mixpack(tpack,n,v)
! !USES:
use modmain
use moddftu
! !INPUT/OUTPUT PARAMETERS:
!   tpack : .true. for packing, .false. for unpacking (in,logical)
!   n     : total number of real values stored (out,integer)
!   v     : packed potential (inout,real(*))
! !DESCRIPTION:
!   Packs/unpacks the muffin-tin and interstitial parts of the Kohn-Sham
!   potential and magnetic field into/from the single array {\tt v}. This array
!   can then be passed directly to the mixing routine. See routine {\tt rfpack}.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
logical, intent(in) :: tpack
integer, intent(out) :: n
real(8), intent(inout) :: v(*)
! local variables
integer idm,ias,lm1,lm2
integer ispn,jspn
n=0
call rfpack(tpack,n,nrmt,nrmtinr,nrmtmax,vsmt,vsir,v)
do idm=1,ndmag
  call rfpack(tpack,n,nrcmt,nrcmtinr,nrcmtmax,bsmt(:,:,:,idm),bsir(:,idm),v)
end do
! pack the DFT+U potential if required
if (tvmatmt) then
  do ias=1,natmtot
    do ispn=1,nspinor
      do jspn=1,nspinor
        do lm1=1,lmmaxdm
          do lm2=1,lmmaxdm
            n=n+1
            if (tpack) then
              v(n)=dble(vmatmt(lm1,ispn,lm2,jspn,ias))
              n=n+1
              v(n)=aimag(vmatmt(lm1,ispn,lm2,jspn,ias))
            else
              vmatmt(lm1,ispn,lm2,jspn,ias)=cmplx(v(n),v(n+1),8)
              n=n+1
            end if
          end do
        end do
      end do
    end do
  end do
end if
return
end subroutine
!EOC

