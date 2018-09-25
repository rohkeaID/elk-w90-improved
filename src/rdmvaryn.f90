
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: rdmvaryn
! !INTERFACE:
subroutine rdmvaryn
! !USES:
use modmain
use modrdm
use modmpi
! !DESCRIPTION:
!   Calculates new occupation numbers from old by using the derivatives of the
!   total energy: $n_i^{\rm new} = n_i^{\rm old}-\tau \gamma_i$, where $\tau$ is
!   chosen such that $0 \le n_i \le n_{\rm max}$ with
!   $$ \gamma_i=\begin{cases}
!    g_i(n_{\rm max}-n_i) & g_i > 0 \\
!    g_i n_i & g_i\le 0 \end{cases} $$
!   where $g_i=\partial E/\partial n_i-\kappa$, and $\kappa$ is chosen such that
!   $\sum_i\gamma_i=0$.
!
! !REVISION HISTORY:
!   Created 2009 (JKD,Sharma)
!EOP
!BOC
implicit none
! local variables
integer, parameter :: maxit=10000
integer it,ik,ist
real(8), parameter :: eps=1.d-12
real(8) tau,sum,gs,gsp,dgs
real(8) kapa,dkapa,t1
! allocatable arrays
real(8), allocatable :: dedn(:,:),gamma(:,:)
! add constant to occupancies for charge conservation
sum=0.d0
do ik=1,nkpt
  do ist=1,nstsv
    sum=sum+wkpt(ik)*occsv(ist,ik)
  end do
end do
t1=(chgval-sum)/dble(nstsv)
occsv(:,:)=occsv(:,:)+t1
! redistribute charge so that occupancies are in the interval [0,occmax]
sum=0.d0
do ik=1,nkpt
  do ist=1,nstsv
    if (occsv(ist,ik).gt.occmax) then
      sum=sum+wkpt(ik)*(occsv(ist,ik)-occmax)
      occsv(ist,ik)=occmax
    end if
    if (occsv(ist,ik).lt.0.d0) then
      sum=sum+wkpt(ik)*occsv(ist,ik)
      occsv(ist,ik)=0.d0
    end if
  end do
end do
do ist=1,nstsv
  do ik=1,nkpt
    if (sum.gt.0.d0) then
      t1=wkpt(ik)*(occmax-occsv(ist,ik))
      t1=min(t1,sum)
      occsv(ist,ik)=occsv(ist,ik)+t1/wkpt(ik)
      sum=sum-t1
    else
      t1=wkpt(ik)*occsv(ist,ik)
      t1=min(t1,-sum)
      occsv(ist,ik)=occsv(ist,ik)-t1/wkpt(ik)
      sum=sum+t1
    end if
  end do
end do
allocate(dedn(nstsv,nkpt))
allocate(gamma(nstsv,nkpt))
! get the derivatives
call rdmdedn(dedn)
! find suitable value of kapa such that sum of gamma is 0
gsp=0.d0
kapa=0.d0
dkapa=0.1d0
do it=1,maxit
  gs=0.d0
  sum=0.d0
  do ik=1,nkpt
    do ist=1,nstsv
      t1=dedn(ist,ik)-kapa
      if (t1.gt.0.d0) then
        gamma(ist,ik)=t1*(occmax-occsv(ist,ik))
      else
        gamma(ist,ik)=t1*occsv(ist,ik)
      end if
      gs=gs+wkpt(ik)*gamma(ist,ik)
      sum=sum+wkpt(ik)*gamma(ist,ik)**2
    end do
  end do
  sum=sqrt(sum)
  sum=max(sum,1.d0)
  t1=abs(gs)/sum
  if (t1.lt.eps) goto 10
  if (it.ge.2) then
    dgs=gs-gsp
    if (gs*dgs.gt.0.d0) dkapa=-dkapa
    if (gs*gsp.lt.0.d0) then
      dkapa=0.5d0*dkapa
    else
      dkapa=1.1d0*dkapa
    end if
  end if
  gsp=gs
  kapa=kapa+dkapa
end do
write(*,*)
write(*,'("Error(rdmvaryn): could not find offset")')
write(*,*)
stop
10 continue
! write derivatives and occupancies to file
call rdmwritededn(dedn)
deallocate(dedn)
! normalize gamma if sum of squares is greater than 1
sum=0.d0
do ik=1,nkpt
  do ist=1,nstsv
    sum=sum+wkpt(ik)*gamma(ist,ik)**2
  end do
end do
if (sum.gt.1.d0) then
  t1=1.d0/sqrt(sum)
  gamma(:,:)=t1*gamma(:,:)
end if
! find step size which keeps occupancies in the interval [0,occmax]
tau=taurdmn
20 continue
if (abs(tau).lt.eps) goto 30
do ik=1,nkpt
  do ist=1,nstsv
    t1=occsv(ist,ik)+tau*gamma(ist,ik)
    if (gamma(ist,ik).gt.0.d0) then
      if (t1.gt.occmax+eps) then
        tau=0.75d0*tau
        goto 20
      end if
    end if
    if (gamma(ist,ik).lt.0.d0) then
      if (t1.lt.-eps) then
        tau=0.75d0*tau
        goto 20
      end if
    end if
  end do
end do
30 continue
! update occupancies and write to file
do ik=1,nkpt
  do ist=1,nstsv
    occsv(ist,ik)=occsv(ist,ik)+tau*gamma(ist,ik)
  end do
  call putoccsv(filext,ik,occsv(:,ik))
end do
deallocate(gamma)
return
end subroutine
!EOC

