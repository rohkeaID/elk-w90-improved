
! This routine is based on code written by K. Burke.

subroutine x_pbe(kappa,mu,rho,s,u,v,ex,vx)
implicit none
! arguments
real(8), intent(in) :: kappa,mu
real(8), intent(in) :: rho,s,u,v
real(8), intent(out) :: ex,vx
! local variables
real(8), parameter :: ax=-0.7385587663820224058d0
real(8), parameter :: thrd=1.d0/3.d0
real(8), parameter :: thrd4=4.d0/3.d0
real(8) ul,exu,s2,p0
real(8) fxpbe,fs,fss
ul=mu/kappa
! LDA exchange energy density
exu=ax*rho**thrd
! PBE enhancement factor
s2=s**2
p0=1.d0+ul*s2
fxpbe=1.d0+kappa-kappa/p0
ex=exu*fxpbe
fs=2.d0*kappa*ul/(p0*p0)
fss=-4.d0*ul*s*fs/p0
! exchange potential
vx=exu*(thrd4*fxpbe-(u-thrd4*s2*s)*fss-v*fs)
return
end subroutine

