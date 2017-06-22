
! Copyright (C) 2002-2011 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: xc_pwca
! !INTERFACE:
subroutine xc_pwca(n,rhoup,rhodn,ex,ec,vxup,vxdn,vcup,vcdn)
! !INPUT/OUTPUT PARAMETERS:
!   n     : number of density points (in,integer)
!   rhoup : spin-up charge density (in,real(n))
!   rhodn : spin-down charge density (in,real(n))
!   ex    : exchange energy density (out,real(n))
!   ec    : correlation energy density (out,real(n))
!   vxup  : spin-up exchange potential (out,real(n))
!   vxdn  : spin-down exchange potential (out,real(n))
!   vcup  : spin-up correlation potential (out,real(n))
!   vcdn  : spin-down correlation potential (out,real(n))
! !DESCRIPTION:
!   Spin-polarised exchange-correlation potential and energy of the Perdew-Wang
!   parameterisation of the Ceperley-Alder electron gas: {\it Phys. Rev. B}
!   {\bf 45}, 13244 (1992) and {\it Phys. Rev. Lett.} {\bf 45}, 566 (1980).
!
! !REVISION HISTORY:
!   Created January 2004 (JKD)
!   Rewrote, October 2011 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: rhoup(n),rhodn(n)
real(8), intent(out) :: ex(n),ec(n)
real(8), intent(out) :: vxup(n),vxdn(n)
real(8), intent(out) :: vcup(n),vcdn(n)
! local variables
integer i
real(8), parameter :: pi=3.1415926535897932385d0
real(8), parameter :: thrd=1.d0/3.d0, thrd4=4.d0/3.d0
real(8), parameter :: d2f0=1.709921d0
real(8) a(3),a1(3),b1(3),b2(3),b3(3),b4(3)
data a  / 0.0310907d0, 0.01554535d0, 0.0168869d0 /
data a1 / 0.21370d0,   0.20548d0,    0.11125d0   /
data b1 / 7.5957d0,   14.1189d0,    10.357d0     /
data b2 / 3.5876d0,    6.1977d0,     3.6231d0    /
data b3 / 1.6382d0,    3.3662d0,     0.88026d0   /
data b4 / 0.49294d0,   0.62517d0,    0.49671d0   /
real(8) p1,p2,p3,rup,rdn,r,ri,ri2
real(8) rs,rs2,rs12,rs32,rsi,rs12i
real(8) mz,z,z3,z4,drs,dzu,dzd
real(8) fz,dfz,ders,dez,deu,ded
real(8) a2,ec0,dec0,ec1,dec1,ac,dac
real(8) t1,t2,t3,t4,t5,t6,t7,dt1,dt2
if (n.le.0) then
  write(*,*)
  write(*,'("Error(xc_pwca): invalid n : ",I8)') n
  write(*,*)
  stop
end if
! prefactors
t1=3.d0/(4.d0*pi)
p1=t1**thrd
p2=t1*(9.d0*pi/4.d0)**thrd
p3=1.d0/(2.d0**thrd4-2.d0)
do i=1,n
  rup=rhoup(i); rdn=rhodn(i)
! total density
  r=rup+rdn
  if ((rup.lt.0.d0).or.(rdn.lt.0.d0).or.(r.lt.1.d-20)) then
    ex(i)=0.d0
    ec(i)=0.d0
    vxup(i)=0.d0
    vxdn(i)=0.d0
    vcup(i)=0.d0
    vcdn(i)=0.d0
    cycle
  end if
  ri=1.d0/r
  ri2=ri**2
  rs=p1*ri**thrd
  rs2=rs**2
  rs12=sqrt(rs)
  rs32=rs12*rs
  rsi=1.d0/rs
  rs12i=1.d0/rs12
  mz=rup-rdn
  z=mz/r
  z3=z**3
  z4=z3*z
! drs/drup = drs/drdn = drs/drho
  drs=-thrd*rs*ri
! dz/drup, dz/drdn
  t1=mz*ri2
  dzu=ri-t1
  dzd=-ri-t1
!------------------!
!     exchange     !
!------------------!
  t1=-p2*rsi/2.d0
  t2=1.d0+z
  t3=1.d0-z
  t4=t2**thrd4
  t5=t3**thrd4
  t6=t4+t5
! exchange energy density
  ex(i)=t1*t6
! dex/drs
  ders=-ex(i)*rsi
! f(z)
  fz=p3*(t6-2.d0)
! dex/dz
  t4=t4/t2
  t5=t5/t3
  t6=t4-t5
  t7=thrd4*t6
  dez=t1*t7
! df/dz
  dfz=p3*t7
! dex/drup, dex/drdn
  t1=ders*drs
  deu=t1+dez*dzu
  ded=t1+dez*dzd
! exchange potential
  vxup(i)=ex(i)+r*deu
  vxdn(i)=ex(i)+r*ded
!---------------------!
!     correlation     !
!---------------------!
! ec(rs,0)
  a2=2.d0*a(1)
  t1=a2*(b1(1)*rs12+b2(1)*rs+b3(1)*rs32+b4(1)*rs2)
  dt1=a2*(0.5d0*b1(1)*rs12i+b2(1)+1.5d0*b3(1)*rs12+2.d0*b4(1)*rs)
  t3=1.d0/t1
  t2=1.d0+t3
  dt2=-dt1*t3**2
  t3=1.d0/t2
  t4=1.d0+a1(1)*rs
  t5=log(t2)
  ec0=-a2*t4*t5
  dec0=-a2*(a1(1)*t5+t4*t3*dt2)
! ec(rs,1)
  a2=2.d0*a(2)
  t1=a2*(b1(2)*rs12+b2(2)*rs+b3(2)*rs32+b4(2)*rs2)
  dt1=a2*(0.5d0*b1(2)*rs12i+b2(2)+1.5d0*b3(2)*rs12+2.d0*b4(2)*rs)
  t3=1.d0/t1
  t2=1.d0+t3
  dt2=-dt1*t3**2
  t3=1.d0/t2
  t4=1.d0+a1(2)*rs
  t5=log(t2)
  ec1=-a2*t4*t5
  dec1=-a2*(a1(2)*t5+t4*t3*dt2)
! ac(rs)
  a2=2.d0*a(3)
  t1=a2*(b1(3)*rs12+b2(3)*rs+b3(3)*rs32+b4(3)*rs2)
  dt1=a2*(0.5d0*b1(3)*rs12i+b2(3)+1.5d0*b3(3)*rs12+2.d0*b4(3)*rs)
  t3=1.d0/t1
  t2=1.d0+t3
  dt2=-dt1*t3**2
  t3=1.d0/t2
  t4=1.d0+a1(3)*rs
  t5=log(t2)
  ac=a2*t4*t5
  dac=a2*(a1(3)*t5+t4*t3*dt2)
! correlation energy density
  t1=1.d0-z4
  t2=(fz/d2f0)*t1
  t3=ec1-ec0
  t4=fz*z4
  ec(i)=ec0+ac*t2+t3*t4
! dec/drs
  t5=dec1-dec0
  ders=dec0+dac*t2+t5*t4
! dec/dz
  t6=4.d0*fz*z3
  dez=(ac/d2f0)*(dfz*t1-t6)+t3*(dfz*z4+t6)
! dec/drup, dec/drdn
  t1=ders*drs
  deu=t1+dez*dzu
  ded=t1+dez*dzd
! correlation potential
  vcup(i)=ec(i)+r*deu
  vcdn(i)=ec(i)+r*ded
end do
return
end subroutine
!EOC

