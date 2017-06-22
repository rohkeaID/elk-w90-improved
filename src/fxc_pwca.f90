
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine fxc_pwca(n,rhoup,rhodn,fxcuu,fxcud,fxcdd)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: rhoup(n)
real(8), intent(in) :: rhodn(n)
real(8), intent(out) :: fxcuu(n)
real(8), intent(out) :: fxcud(n)
real(8), intent(out) :: fxcdd(n)
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
real(8) p1,p2,p3,rup,rdn,r,ri,ri2,ri3
real(8) rs,rs2,rs12,rs32,rsi,rs12i,rs32i
real(8) mz,z,z2,z3,z4,fz,dfz,d2fz
real(8) drs,d2rs,dzu,d2zu,dzd,d2zd,d2zud
real(8) ders,d2ers,dez,d2ez,d2ersz
real(8) deu,d2eu,ded,d2ed,d2eud,ex
real(8) ec0,dec0,d2ec0,ec1,dec1,d2ec1
real(8) ac,dac,d2ac,a2,dt1,d2t1,dt2,d2t2
real(8) t1,t2,t3,t4,t5,t6,t7,t8
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
    fxcuu(i)=0.d0
    fxcud(i)=0.d0
    fxcdd(i)=0.d0
    cycle
  end if
  ri=1.d0/r
  ri2=ri**2
  ri3=ri2*ri
  rs=p1*ri**thrd
  rs2=rs**2
  rs12=sqrt(rs)
  rs32=rs12*rs
  rsi=1.d0/rs
  rs12i=1.d0/rs12
  rs32i=1.d0/rs32
  mz=rup-rdn
  z=mz/r
  z2=z**2
  z3=z2*z
  z4=z3*z
! drs/drup = drs/drdn = drs/drho
  drs=-thrd*rs*ri
! d2rs/drup^2 = d^2rs/drn^2 = d^2rs/drho^2
  d2rs=-thrd4*drs*ri
! dz/drup, dz/drdn
  t1=mz*ri2
  dzu=ri-t1
  dzd=-ri-t1
! d^2z/drup^2, d^2z/drdn^2, d^2z/drup*drdn
  t1=2.d0*mz*ri3
  t2=2.d0*ri2
  d2zu=t1-t2
  d2zd=t1+t2
  d2zud=t1
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
  ex=t1*t6
! dex/drs
  ders=-ex*rsi
! d^2ex/drs^2
  d2ers=-2.d0*ders*rsi
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
! d^2ex/dz^2
  t4=t4/t2
  t5=t5/t3
  t6=t4+t5
  t7=thrd4*thrd*t6
  d2ez=t1*t7
! d^2f/dz^2
  d2fz=p3*t7
! d^2ex/drs*dz
  d2ersz=-dez*rsi
! dex/drup, dex/drdn
  t1=ders*drs
  deu=t1+dez*dzu
  ded=t1+dez*dzd
! d^2ex/drup^2
  t1=d2ers*drs
  t2=d2ersz*drs
  t3=ders*d2rs
  t4=(t1+d2ersz*dzu)*drs+t3
  t5=t2+d2ez*dzu
  d2eu=t4+t5*dzu+dez*d2zu
! d^2ex/drdn^2
  d2ed=(t1+d2ersz*dzd)*drs+t3+(t2+d2ez*dzd)*dzd+dez*d2zd
! d^2ex/drup*drdn
  d2eud=t4+t5*dzd+dez*d2zud
! calculate fxc
  fxcuu(i)=2.d0*deu+r*d2eu
  fxcud(i)=deu+ded+r*d2eud
  fxcdd(i)=2.d0*ded+r*d2ed
!---------------------!
!     correlation     !
!---------------------!
! ec(rs,0)
  a2=2.d0*a(1)
  t1=a2*(b1(1)*rs12+b2(1)*rs+b3(1)*rs32+b4(1)*rs2)
  dt1=a2*(0.5d0*b1(1)*rs12i+b2(1)+1.5d0*b3(1)*rs12+2.d0*b4(1)*rs)
  d2t1=a2*(-0.25d0*b1(1)*rs32i+0.75d0*b3(1)*rs12i+2.d0*b4(1))
  t3=1.d0/t1
  t4=t3**2
  t2=1.d0+t3
  dt2=-dt1*t4
  d2t2=t4*(2.d0*t3*dt1**2-d2t1)
  t3=1.d0/t2
  t4=1.d0+a1(1)*rs
  t5=log(t2)
  ec0=-a2*t4*t5
  dec0=-a2*(a1(1)*t5+t4*t3*dt2)
  d2ec0=-a2*(2.d0*a1(1)*t3*dt2+t4*t3*(d2t2-t3*dt2**2))
! ec(rs,1)
  a2=2.d0*a(2)
  t1=a2*(b1(2)*rs12+b2(2)*rs+b3(2)*rs32+b4(2)*rs2)
  dt1=a2*(0.5d0*b1(2)*rs12i+b2(2)+1.5d0*b3(2)*rs12+2.d0*b4(2)*rs)
  d2t1=a2*(-0.25d0*b1(2)*rs32i+0.75d0*b3(2)*rs12i+2.d0*b4(2))
  t3=1.d0/t1
  t4=t3**2
  t2=1.d0+t3
  dt2=-dt1*t4
  d2t2=t4*(2.d0*t3*dt1**2-d2t1)
  t3=1.d0/t2
  t4=1.d0+a1(2)*rs
  t5=log(t2)
  ec1=-a2*t4*t5
  dec1=-a2*(a1(2)*t5+t4*t3*dt2)
  d2ec1=-a2*(2.d0*a1(2)*t3*dt2+t4*t3*(d2t2-t3*dt2**2))
! ac(rs)
  a2=2.d0*a(3)
  t1=a2*(b1(3)*rs12+b2(3)*rs+b3(3)*rs32+b4(3)*rs2)
  dt1=a2*(0.5d0*b1(3)*rs12i+b2(3)+1.5d0*b3(3)*rs12+2.d0*b4(3)*rs)
  d2t1=a2*(-0.25d0*b1(3)*rs32i+0.75d0*b3(3)*rs12i+2.d0*b4(3))
  t3=1.d0/t1
  t4=t3**2
  t2=1.d0+t3
  dt2=-dt1*t4
  d2t2=t4*(2.d0*t3*dt1**2-d2t1)
  t3=1.d0/t2
  t4=1.d0+a1(3)*rs
  t5=log(t2)
  ac=a2*t4*t5
  dac=a2*(a1(3)*t5+t4*t3*dt2)
  d2ac=a2*(2.d0*a1(3)*t3*dt2+t4*t3*(d2t2-t3*dt2**2))
! correlation energy density derivatives
  t1=1.d0-z4
  t2=(fz/d2f0)*t1
  t3=ec1-ec0
  t4=fz*z4
! dec/drs
  t5=dec1-dec0
  ders=dec0+dac*t2+t5*t4
! d^2ec/drs^2
  t6=d2ec1-d2ec0
  d2ers=d2ec0+d2ac*t2+t6*t4
! dec/dz
  t4=ac/d2f0
  t6=4.d0*fz*z3
  t7=dfz*t1-t6
  t8=dfz*z4+t6
  dez=t4*t7+t3*t8
! d^2ec/drs*dz
  d2ersz=(dac/d2f0)*t7+t5*t8
! d^2ec/dz^2
  t7=8.d0*dfz*z3
  t8=12.d0*fz*z2
  d2ez=t4*(d2fz*t1-t7-t8)+t3*(d2fz*z4+t7+t8)
! dec/drup, dec/drdn
  t1=ders*drs
  deu=t1+dez*dzu
  ded=t1+dez*dzd
! d^2ec/drup^2
  t1=d2ers*drs
  t2=d2ersz*drs
  t3=ders*d2rs
  t4=(t1+d2ersz*dzu)*drs+t3
  t5=t2+d2ez*dzu
  d2eu=t4+t5*dzu+dez*d2zu
! d^2ec/drdn^2
  d2ed=(t1+d2ersz*dzd)*drs+t3+(t2+d2ez*dzd)*dzd+dez*d2zd
! d^2ec/drup*drdn
  d2eud=t4+t5*dzd+dez*d2zud
! calculate fxc
  fxcuu(i)=fxcuu(i)+2.d0*deu+r*d2eu
  fxcud(i)=fxcud(i)+deu+ded+r*d2eud
  fxcdd(i)=fxcdd(i)+2.d0*ded+r*d2ed
end do
return
end subroutine

