
real(8) function stheta_lr(x)
implicit none
! arguments
real(8), intent(in) :: x
! local variables
real(8), parameter :: pi=3.1415926535897932385d0
stheta_lr=0.5d0+atan(2*x)/pi
return
end function

