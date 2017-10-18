
! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: genppts
! !INTERFACE:
subroutine genppts(tfbz,nsym,sym,ngridp,npptnr,epslat,bvec,boxl,nppt,ipvip, &
 ipvipnr,ivp,vpl,vpc,wppt,wpptnr)
! !INPUT/OUTPUT PARAMETERS:
!   tfbz    : .true. if vpl and vpc should be mapped to the first Brillouin zone
!             (in,logical)
!   nsym    : number of point group symmetries used for reduction, set to 1 for
!             no reduction (in,integer)
!   sym     : symmetry matrices in lattice coordinates (in,integer(3,3,*))
!   ngridp  : p-point grid sizes (in,integer(3))
!   npptnr  : number of non-reduced p-points: ngridp(1)*ngridp(2)*ngridp(3)
!             (in,integer)
!   epslat  : tolerance for determining identical vectors (in,real)
!   bvec    : reciprocal lattice vectors (in,real(3,3))
!   boxl    : corners of box containing p-points in lattice coordinates, the
!             zeroth vector is the origin (in,real(3,0:3))
!   nppt    : total number of p-points (out,integer)
!   ipvip   : map from (i1,i2,i3) to reduced p-point index
!             (out,integer(0:ngridp(1)-1,0:ngridp(2)-1,0:ngridp(3)-1))
!   ipvipnr : map from (i1,i2,i3) to non-reduced p-point index
!             (out,integer(0:ngridp(1)-1,0:ngridp(2)-1,0:ngridp(3)-1))
!   ivp     : integer coordinates of the p-points
!             (out,integer(3,ngridp(1)*ngridp(2)*ngridp(3)))
!   vpl     : lattice coordinates of each p-point
!             (out,real(3,ngridp(1)*ngridp(2)*ngridp(3)))
!   vpc     : Cartesian coordinates of each p-point
!             (out,real(3,ngridp(1)*ngridp(2)*ngridp(3)))
!   wppt    : weights of each reduced p-point
!             (out,real(ngridp(1)*ngridp(2)*ngridp(3)))
!   wpptnr  : weight of each non-reduced p-point (out,real)
! !DESCRIPTION:
!   This routine is used for generating $k$-point or $q$-point sets. Since these
!   are stored in global arrays, the points passed to this and other routines
!   are referred to as $p$-points. In lattice coordinates, the ${\bf p}$ vectors
!   are given by
!   $$ {\bf p}=\left(\begin{matrix} & & \\
!     {\bf B}_2-{\bf B}_1 & {\bf B}_3-{\bf B}_1 & {\bf B}_4-{\bf B}_1 \\
!       & & \end{matrix}\right)
!     \left(\begin{matrix}i_1/n_1 \\ i_2/n_2 \\ i_3/n_3 \end{matrix}\right)
!     +{\bf B}_1 $$
!   where $i_j$ runs from 0 to $n_j-1$, and the ${\bf B}$ vectors define the
!   corners of a box with ${\bf B}_1$ as the origin. If {\tt tfbz} is
!   {\tt .true.} then each {\tt vpl} vector is mapped to the first Brillouin
!   zone. If {\tt tfbz} is {\tt .false.} and {\tt nsym} $>0$, then the
!   coordinates of each {\tt vpl} are mapped to the $[0,1)$ interval. The
!   $p$-point weights are stored in {\tt wppt} and the array {\tt ipvip}
!   contains the map from the integer coordinates to the reduced index.
!
! !REVISION HISTORY:
!   Created August 2002 (JKD)
!   Updated April 2007 (JKD)
!   Added mapping to the first Brillouin zone, September 2008 (JKD)
!   Made independent of modmain, February 2010 (JKD)
!EOP
!BOC
implicit none
! arguments
logical, intent(in) :: tfbz
integer, intent(in) :: nsym,sym(3,3,*)
integer, intent(in) :: ngridp(3),npptnr
real(8), intent(in) :: epslat
real(8), intent(in) :: bvec(3,3),boxl(3,0:3)
integer, intent(out) :: nppt
integer, intent(out) :: ipvip(0:ngridp(1)-1,0:ngridp(2)-1,0:ngridp(3)-1)
integer, intent(out) :: ipvipnr(0:ngridp(1)-1,0:ngridp(2)-1,0:ngridp(3)-1)
integer, intent(out) :: ivp(3,npptnr)
real(8), intent(out) :: vpl(3,npptnr),vpc(3,npptnr)
real(8), intent(out) :: wppt(npptnr),wpptnr
! local variables
integer i1,i2,i3,i
integer isym,ip,jp
real(8) v1(3),v2(3),v3(3),v4(3)
real(8) b(3,3),t1
if ((ngridp(1).le.0).or.(ngridp(2).le.0).or.(ngridp(3).le.0)) then
  write(*,*)
  write(*,'("Error(genppts): invalid ngridp : ",3I8)') ngridp
  write(*,*)
  stop
end if
if (npptnr.ne.ngridp(1)*ngridp(2)*ngridp(3)) then
  write(*,*)
  write(*,'("Error(genppts): mismatched npptnr and ngridp : ",4I8)') npptnr, &
   ngridp
  write(*,*)
  stop
end if
! box vector matrix
b(:,1)=boxl(:,1)-boxl(:,0)
b(:,2)=boxl(:,2)-boxl(:,0)
b(:,3)=boxl(:,3)-boxl(:,0)
! weight of each non-reduced p-point
wpptnr=1.d0/dble(ngridp(1)*ngridp(2)*ngridp(3))
ip=0
jp=npptnr+1
do i3=0,ngridp(3)-1
  v1(3)=dble(i3)/dble(ngridp(3))
  do i2=0,ngridp(2)-1
    v1(2)=dble(i2)/dble(ngridp(2))
    do i1=0,ngridp(1)-1
      v1(1)=dble(i1)/dble(ngridp(1))
      call r3mv(b,v1,v2)
      v2(:)=v2(:)+boxl(:,0)
      if (nsym.gt.1) then
! determine if this point is equivalent to one already in the set
        do isym=1,nsym
          v3(:)=sym(1,:,isym)*v2(1)+sym(2,:,isym)*v2(2)+sym(3,:,isym)*v2(3)
          call r3frac(epslat,v3)
          do i=1,ip
            v4(:)=vpl(:,i)
            call r3frac(epslat,v4)
            t1=abs(v4(1)-v3(1))+abs(v4(2)-v3(2))+abs(v4(3)-v3(3))
            if (t1.lt.epslat) then
! equivalent p-point found so add to existing weight
              ipvip(i1,i2,i3)=i
              wppt(i)=wppt(i)+wpptnr
! add new point to back of set
              jp=jp-1
              ipvipnr(i1,i2,i3)=jp
              ivp(1,jp)=i1; ivp(2,jp)=i2; ivp(3,jp)=i3
              vpl(:,jp)=v2(:)
              wppt(jp)=0.d0
              goto 10
            end if
          end do
        end do
      end if
! add new point to set
      ip=ip+1
      ipvip(i1,i2,i3)=ip
      ipvipnr(i1,i2,i3)=ip
      ivp(1,ip)=i1; ivp(2,ip)=i2; ivp(3,ip)=i3
      vpl(:,ip)=v2(:)
      wppt(ip)=wpptnr
10 continue
    end do
  end do
end do
nppt=ip
do ip=1,npptnr
! map vpl to the first Brillouin zone if required
  if (tfbz) call vecfbz(epslat,bvec,vpl(:,ip))
! determine the Cartesian coordinates of the p-points
  call r3mv(bvec,vpl(:,ip),vpc(:,ip))
end do
return
end subroutine
!EOC

