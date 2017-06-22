
! Copyright (C) 2008 F. Bultmark, F. Cricchio and L. Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writetm2du
! !INTERFACE:
subroutine writetm2du(fnum)
! !USES:
use modmain
use moddftu
use modtest
! !DESCRIPTION:
!   Decompose the density matrix and the DFT+$U$ Hartree-Fock energy in 2-index
!   tensor moments components, see {\it Phys. Rev. B} {\bf 80}, 035121 (2009).
!
! !REVISION HISTORY:
!   Created October 2009 (F. Cricchio and L. Nordstrom)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: fnum
! local variables
integer is,ia,ias,i
integer l,p,k,x,y
real(8) ehf,vh,vx,w2,t1
real(8) tm2ptot,tm2p,tm2p0
! allocatable arrays
real(8), allocatable :: ex(:,:,:),tm22(:,:,:)
complex(8), allocatable :: tm2(:,:)
if (.not.spinpol) then
  write(*,*)
  write(*,'("Error(writetm2du): spin-polarisation must be enabled")')
  write(*,*)
  stop
end if
if (iscl.le.1) then
  write(fnum,*)
  write(fnum,'("Tensor moment decomposition of density matrix and Hartree-Fock &
   &energy")')
  write(fnum,'(" (see Phys. Rev. B. 81, 140403(R) (2010))")')
  write(fnum,'("  W.W = modulus square of tensor moment")')
  write(fnum,'("  Eh = Hartree energy term")')
  write(fnum,'("  Ex = exchange energy term")')
  write(fnum,'("  Pol = polarisation of density matrix")')
end if
if (iscl.ge.1) then
  write(fnum,*)
  write(fnum,'("+--------------------+")')
  write(fnum,'("| Loop number : ",I4," |")') iscl
  write(fnum,'("+--------------------+")')
end if
! allocate arrays
allocate(ex(0:2*lmaxdm,0:1,natmtot))
allocate(tm22(0:2*lmaxdm,0:1,natmtot))
allocate(tm2(-lmmaxdm:lmmaxdm,-1:1))
tm22(:,:,:)=0.d0
ex(:,:,:)=0.d0
tm2p0=0.d0
! loop over DFT+U entries
do i=1,ndftu
  is=idftu(1,i)
  l=idftu(2,i)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! write info to TMDFTU.OUT
    write(fnum,*)
    write(fnum,'("Species : ",I4," (",A,"), atom : ",I4)') is, &
     trim(spsymb(is)),ia
    write(fnum,'(" l = ",I2)') l
! zero total Hartree-Fock energy
    ehf=0.d0
! zero total polarisation
    tm2ptot=0.d0
    if (mp_mpi) write(fnum,*)
    do k=0,2*l
      do p=0,1
! decompose density matrix in 2-index tensor moment components
        call dmtotm2(l,nspinor,k,p,lmmaxdm,dmatmt(:,:,:,:,ias),tm2)
! calculate modulus square of tensor moment
        w2=0.d0
        do x=-k,k
          do y=-p,p
            w2=w2+dble((-1)**(x+y))*dble(tm2(-x,-y)*tm2(x,y))
          end do
        end do
! potential energy components
        call pottm2(i,k,p,vh,vx)
! save square of tensor modulus and exchange energy
        tm22(k,p,ias)=w2
        ex(k,p,ias)=vx*w2
! polarisation terms
        call tm2pol(l,k,w2,tm2p)
! write square of tensor modulus; Hartree and exchange energy; and polarisation
        write(fnum,'("  k = ",I2,", p = ",I2)') k,p
        if ((k+p).eq.0) then
! for k,p = 0 save reference polarisation
          tm2p0=tm2p
! for k,p = 0 do not write out the polarisation
          write(fnum,'("   W.W =",F14.8,", Eh =",F14.8,", Ex =",F14.8)') w2, &
           vh*w2,vx*w2
        else
! relative polarisation
          t1=tm2p/tm2p0
          write(fnum,'("   W.W =",F14.8,", Eh =",F14.8,", Ex =",F14.8,&
           &", Pol =",F14.8)') w2,vh*w2,vx*w2,t1
! total relative polarisation (skipping 000 component)
          tm2ptot=tm2ptot+t1
        end if
        ehf=ehf+(vh+vx)*w2
        write(fnum,*)
        do x=-k,k
          do y=-p,p
! write out single components of tensor moments
            write(fnum,'("    x = ",I2,"    y = ",I2," : ",2F14.8)') x,y, &
             tm2(x,y)
          end do
        end do
        write(fnum,*)
      end do
    end do
    write(fnum,*)
    write(fnum,'("  Total Hartree-Fock energy (without DC correction) : ",&
     &F14.8)') ehf
    write(fnum,'("  Total polarisation of density matrix : ",F14.8)') tm2ptot
    write(fnum,*)
! end loop over atoms and species
  end do
end do
call flushifc(fnum)
! write test files if required
if (test) then
  t1=sqrt(sum(tm22(:,:,:)**2))
  call writetest(820,'Norm of tensor moments',tol=1.d-4,rv=t1)
  t1=sqrt(sum(ex(:,:,:)**2))
  call writetest(830,'Norm of DFT+U Hartree-Fock exchange energies',tol=1.d-4, &
   rv=t1)
end if
deallocate(ex,tm22,tm2)
return
end subroutine
!EOC

