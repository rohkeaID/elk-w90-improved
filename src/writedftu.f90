
! Copyright (C) 2007 F. Bultmark, F. Cricchio and L. Nordstrom.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine writedftu
use modmain
use moddftu
implicit none
! local variables
integer ispn,jspn,is,ia,ias
integer i,k,l,m1,m2,lm1,lm2
if (dftu.eq.0) return
! machine readable density matrix file
open(50,file='DMATMT'//trim(filext),action='WRITE',form='FORMATTED')
do i=1,ndftu
  is=idftu(1,i)
  l=idftu(2,i)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    write(50,*)
    write(50,*)
    write(50,'(3I4," : species, atom, l")') is,ia,l
    do ispn=1,nspinor
      do jspn=1,nspinor
        write(50,*)
        write(50,'(2I4," : ispn, jspn; m1, m2, dmatmt below")') ispn,jspn
        do m1=-l,l
          lm1=idxlm(l,m1)
          do m2=-l,l
            lm2=idxlm(l,m2)
            write(50,'(2I6," ",2G18.10)') m1,m2,dmatmt(lm1,ispn,lm2,jspn,ias)
          end do
        end do
      end do
    end do
  end do
end do
close(50)
! machine readable potential matrix file
open(50,file='VMATMT'//trim(filext),action='WRITE',form='FORMATTED')
do i=1,ndftu
  is=idftu(1,i)
  l=idftu(2,i)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    write(50,*)
    write(50,*)
    write(50,'(3I4," : species, atom, l")') is,ia,l
    do ispn=1,nspinor
      do jspn=1,nspinor
        write(50,*)
        write(50,'(2I4," : ispn, jspn; m1, m2, vmatmt below")') ispn,jspn
        do m1=-l,l
          lm1=idxlm(l,m1)
          do m2=-l,l
            lm2=idxlm(l,m2)
            write(50,'(2I6," ",2G18.10)') m1,m2,vmatmt(lm1,ispn,lm2,jspn,ias)
          end do
        end do
      end do
    end do
  end do
end do
close(50)
! machine readable alpha parameters
if ((dftu.eq.3).and.(.not.readadu)) then
  open(50,file='ALPHADU'//trim(filext),action='WRITE',form='FORMATTED')
  do i=1,ndftu
    is=idftu(1,i)
    l=idftu(2,i)
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      write(50,*)
      write(50,'(3I4," : species, atom, l")') is,ia,l
      write(50,'(G18.10," : alpha")') alphadu(ia,i)
    end do
  end do
  close(50)
end if
! Slater parameters
open(50,file='FDU'//trim(filext),action='WRITE',form='FORMATTED')
do i=1,ndftu
  is=idftu(1,i)
  l=idftu(2,i)
  write(50,*)
  write(50,'(3I4," : species, l")') is,l
  do k=0,2*l,2
    write(50,'(G18.10," : F^(",I1,")")') fdu(k,i),k
  end do
  write(50,'(G18.10," : U")') ujdu(1,i)
  write(50,'(G18.10," : J")') ujdu(2,i)
  if (inpdftu.ge.4) write(50,'(G18.10," : lambdadu")') lambdadu(i)
end do
close(50)
return
end subroutine

