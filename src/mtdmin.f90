
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: mtdmin
! !INTERFACE:
subroutine mtdmin(is,js,dmin)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   is, js : species numbers (out,integer)
!   dmin   : minimum distance between muffin-tin surfaces (out,real)
! !DESCRIPTION:
!   Finds the atomic species pair for which the distance between the muffin-tin
!   surfaces is a minimum. This distance may be negative if the muffin-tins
!   overlap.
!
! !REVISION HISTORY:
!   Created October 2011 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(out) :: is,js
real(8), intent(out) :: dmin
! local variables
integer i1,i2,i3,ks,ka,ls,la
real(8) v1(3),v2(3),t1,t2,t3
is=1
js=1
dmin=1.d6
do i1=-1,1
  do i2=-1,1
    do i3=-1,1
      v1(:)=dble(i1)*avec(:,1)+dble(i2)*avec(:,2)+dble(i3)*avec(:,3)
      do ks=1,nspecies
        do ka=1,natoms(ks)
          v2(:)=v1(:)+atposc(:,ka,ks)
          do ls=1,nspecies
            t1=rmt(ks)+rmt(ls)
            do la=1,natoms(ls)
              if ((i1.ne.0).or.(i2.ne.0).or.(i3.ne.0).or.(ks.ne.ls).or. &
               (ka.ne.la)) then
                t2=sqrt((v2(1)-atposc(1,la,ls))**2 &
                       +(v2(2)-atposc(2,la,ls))**2 &
                       +(v2(3)-atposc(3,la,ls))**2)
                t3=t2-t1
                if (t3.lt.dmin-epslat) then
                  is=ks
                  js=ls
                  dmin=t3
                end if
              end if
            end do
          end do
        end do
      end do
    end do
  end do
end do
return
end subroutine
!EOC

