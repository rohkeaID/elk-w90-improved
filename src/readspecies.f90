
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine readspecies
use modmain
implicit none
! local variables
integer is,ist,iostat
integer nlx,ilx,lx,ilo
integer io,jo,ko,l,i,j
e0min=0.d0
do is=1,nspecies
  open(50,file=trim(sppath)//trim(spfname(is)),status='OLD',form='FORMATTED', &
   iostat=iostat)
  if (iostat.ne.0) then
    write(*,*)
    write(*,'("Error(readspecies): error opening species file ",A)') &
     trim(sppath)//trim(spfname(is))
    write(*,*)
    stop
  end if
  read(50,*) spsymb(is)
  read(50,*) spname(is)
  read(50,*) spzn(is)
  read(50,*) spmass(is)
  read(50,*) rminsp(is),rmt(is),rmaxsp(is),nrmt(is)
  if (rminsp(is).le.0.d0) then
    write(*,*)
    write(*,'("Error(readspecies): rminsp <= 0 : ",G18.10)') rminsp(is)
    write(*,'(" for species ",I4)') is
    write(*,*)
    stop
  end if
  if (rmt(is).le.rminsp(is)) then
    write(*,*)
    write(*,'("Error(readspecies): rmt <= rminsp : ",2G18.10)') rmt(is), &
     rminsp(is)
    write(*,'(" for species ",I4)') is
    write(*,*)
    stop
  end if
  if (rmaxsp(is).lt.rmt(is)) then
    write(*,*)
    write(*,'("Error(readspecies): rmaxsp < rmt : ",2G18.10)') rmaxsp(is), &
     rmt(is)
    write(*,*)
    stop
  end if
  if (nrmt(is).lt.20) then
    write(*,*)
    write(*,'("Error(readspecies): nrmt too small : ",I8)') nrmt(is)
    write(*,'(" for species ",I4)') is
    write(*,*)
    stop
  end if
! multiply nrmt by the scale factor
  nrmt(is)=nint(nrmt(is)*nrmtscf)
! reduce the minimum radial mesh point by the same factor
  rminsp(is)=rminsp(is)/nrmtscf
  read(50,*) nstsp(is)
  if ((nstsp(is).le.0).or.(nstsp(is).gt.maxstsp)) then
    write(*,*)
    write(*,'("Error(readspecies): nstsp out of range : ",I8)') nstsp(is)
    write(*,'(" for species ",I4)') is
    write(*,*)
    stop
  end if
  do ist=1,nstsp(is)
    read(50,*) nsp(ist,is),lsp(ist,is),ksp(ist,is),occsp(ist,is),spcore(ist,is)
    if (nsp(ist,is).lt.1) then
      write(*,*)
      write(*,'("Error(readspecies): nsp < 1 : ",I8)') nsp(ist,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and state ",I4)') ist
      write(*,*)
      stop
    end if
    if (lsp(ist,is).lt.0) then
      write(*,*)
      write(*,'("Error(readspecies): lsp < 0 : ",I8)') lsp(ist,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and state ",I4)') ist
      write(*,*)
      stop
    end if
    if (ksp(ist,is).lt.1) then
      write(*,*)
      write(*,'("Error(readspecies): ksp < 1 : ",I8)') ksp(ist,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and state ",I4)') ist
      write(*,*)
      stop
    end if
    if (occsp(ist,is).lt.0.d0) then
      write(*,*)
      write(*,'("Error(readspecies): occsp < 0 : ",G18.10)') occsp(ist,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and state ",I4)') ist
      write(*,*)
      stop
    end if
  end do
  read(50,*) apword(0,is)
  if (apword(0,is).le.0) then
    write(*,*)
    write(*,'("Error(readspecies): apword <= 0 : ",I8)') apword(0,is)
    write(*,'(" for species ",I4)') is
    write(*,*)
    stop
  end if
  if (apword(0,is).gt.maxapword) then
    write(*,*)
    write(*,'("Error(readspecies): apword too large : ",I8)') apword(0,is)
    write(*,'(" for species ",I4)') is
    write(*,'("Adjust maxapword in modmain and recompile code")')
    write(*,*)
    stop
  end if
! set the APW orders for l>0
  apword(1:lmaxapw,is)=apword(0,is)
  do io=1,apword(0,is)
    read(50,*) apwe0(io,0,is),apwdm(io,0,is),apwve(io,0,is)
    if (apwdm(io,0,is).lt.0) then
      write(*,*)
      write(*,'("Error(readspecies): apwdm < 0 : ",I8)') apwdm(io,0,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and order ",I4)') io
      write(*,*)
      stop
    end if
! set the APW linearisation energies, derivative orders and variability for l>0
    apwe0(io,1:lmaxapw,is)=apwe0(io,0,is)
    apwdm(io,1:lmaxapw,is)=apwdm(io,0,is)
    apwve(io,1:lmaxapw,is)=apwve(io,0,is)
    e0min=min(e0min,apwe0(io,0,is))
  end do
  read(50,*) nlx
  if (nlx.lt.0) then
    write(*,*)
    write(*,'("Error(readspecies): nlx < 0 : ",I8)') nlx
    write(*,'(" for species ",I4)') is
    write(*,*)
    stop
  end if
  do ilx=1,nlx
    read(50,*) lx,io
    if (lx.lt.0) then
      write(*,*)
      write(*,'("Error(readspecies): lx < 0 : ",I8)') lx
      write(*,'(" for species ",I4)') is
      write(*,'(" and exception number ",I4)') ilx
      write(*,*)
      stop
    end if
    if (lx.gt.lmaxapw) then
      write(*,*)
      write(*,'("Error(readspecies): lx > lmaxapw : ",I8)') lx
      write(*,'(" for species ",I4)') is
      write(*,'(" and exception number ",I4)') ilx
      write(*,*)
      stop
    end if
    apword(lx,is)=io
    if (apword(lx,is).le.0) then
      write(*,*)
      write(*,'("Error(readspecies): apword <= 0 : ",I8)') apword(lx,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and exception number ",I4)') ilx
      write(*,*)
      stop
    end if
    if (apword(lx,is).gt.maxapword) then
      write(*,*)
      write(*,'("Error(readspecies): apword too large : ",I8)') apword(lx,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and exception number ",I4)') ilx
      write(*,'("Adjust maxapword in modmain and recompile code")')
      write(*,*)
      stop
    end if
    do io=1,apword(lx,is)
      read(50,*) apwe0(io,lx,is),apwdm(io,lx,is),apwve(io,lx,is)
      if (apwdm(io,lx,is).lt.0) then
        write(*,*)
        write(*,'("Error(readspecies): apwdm < 0 : ",I8)') apwdm(io,lx,is)
        write(*,'(" for species ",I4)') is
        write(*,'(" exception number ",I4)') ilx
        write(*,'(" and order ",I4)') io
        write(*,*)
        stop
      end if
      e0min=min(e0min,apwe0(io,lx,is))
    end do
  end do
! add excess order to APW functions if required
  if (nxoapwlo.gt.0) then
    do l=0,lmaxapw
      jo=apword(l,is)
      ko=jo+nxoapwlo
      if (ko.gt.maxapword) ko=maxapword
      i=0
      do io=jo+1,ko
        i=i+1
        apwe0(io,l,is)=apwe0(jo,l,is)
        apwdm(io,l,is)=apwdm(jo,l,is)+i
        apwve(io,l,is)=apwve(jo,l,is)
      end do
      apword(l,is)=ko
    end do
  end if
  read(50,*) nlorb(is)
  if (nlorb(is).lt.0) then
    write(*,*)
    write(*,'("Error(readspecies): nlorb < 0 : ",I8)') nlorb(is)
    write(*,'(" for species ",I4)') is
    write(*,*)
    stop
  end if
  if (nlorb(is).gt.maxlorb) then
    write(*,*)
    write(*,'("Error(readspecies): nlorb too large : ",I8)') nlorb(is)
    write(*,'(" for species ",I4)') is
    write(*,'("Adjust maxlorb in modmain and recompile code")')
    write(*,*)
    stop
  end if
  do ilo=1,nlorb(is)
    read(50,*) lorbl(ilo,is),lorbord(ilo,is)
    if (lorbl(ilo,is).lt.0) then
      write(*,*)
      write(*,'("Error(readspecies): lorbl < 0 : ",I8)') lorbl(ilo,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and local-orbital ",I4)') ilo
      write(*,*)
      stop
    end if
    if (lorbl(ilo,is).gt.lmaxapw) then
      write(*,*)
      write(*,'("Error(readspecies): lorbl > lmaxapw : ",2I8)') lorbl(ilo,is), &
       lmaxapw
      write(*,'(" for species ",I4)') is
      write(*,'(" and local-orbital ",I4)') ilo
      write(*,*)
      stop
    end if
    if (lorbord(ilo,is).lt.2) then
      write(*,*)
      write(*,'("Error(readspecies): lorbord < 2 : ",I8)') lorbord(ilo,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and local-orbital ",I4)') ilo
      write(*,*)
      stop
    end if
    if (lorbord(ilo,is).gt.maxlorbord) then
      write(*,*)
      write(*,'("Error(readspecies): lorbord too large : ",I8)') lorbord(ilo,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and local-orbital ",I4)') ilo
      write(*,'("Adjust maxlorbord in modmain and recompile code")')
      write(*,*)
      stop
    end if
    do io=1,lorbord(ilo,is)
      read(50,*) lorbe0(io,ilo,is),lorbdm(io,ilo,is),lorbve(io,ilo,is)
      if (lorbdm(io,ilo,is).lt.0) then
        write(*,*)
        write(*,'("Error(readspecies): lorbdm < 0 : ",I8)') lorbdm(io,ilo,is)
        write(*,'(" for species ",I4)') is
        write(*,'(" local-orbital ",I4)') ilo
        write(*,'(" and order ",I4)') io
        write(*,*)
        stop
      end if
      e0min=min(e0min,lorbe0(io,ilo,is))
    end do
  end do
! add excess local-orbitals if required
  if (nxlo.gt.0) then
    lx=-1
    do ilo=1,nlorb(is)
      do io=1,lorbord(ilo,is)
        if (lorbe0(io,ilo,is).lt.0.d0) goto 10
      end do
      if (lorbl(ilo,is).gt.lx) lx=lorbl(ilo,is)
10 continue
    end do
    ilo=nlorb(is)
    do i=1,nxlo
      if (ilo.eq.maxlorb) exit
      l=lx+i
      if (l.gt.lmaxapw) exit
      ilo=ilo+1
      lorbl(ilo,is)=l
      lorbord(ilo,is)=apword(l,is)+1
      do io=1,lorbord(ilo,is)
        lorbe0(io,ilo,is)=apwe0(1,l,is)
        lorbdm(io,ilo,is)=io-1
        lorbve(io,ilo,is)=apwve(1,l,is)
      end do
    end do
    nlorb(is)=ilo
  end if
! add excess order to local-orbitals if required
  if (nxoapwlo.gt.0) then
    do ilo=1,nlorb(is)
! find the maximum energy derivative
      jo=1
      j=lorbdm(jo,ilo,is)
      do io=1,lorbord(ilo,is)
        i=lorbdm(io,ilo,is)
        if (i.gt.j) then
          jo=io
          j=i
        end if
      end do
      ko=lorbord(ilo,is)+nxoapwlo
      if (ko.gt.maxlorbord) ko=maxlorbord
      i=0
      do io=lorbord(ilo,is)+1,ko
        i=i+1
        lorbe0(io,ilo,is)=lorbe0(jo,ilo,is)
        lorbdm(io,ilo,is)=lorbdm(jo,ilo,is)+i
        lorbve(io,ilo,is)=lorbve(jo,ilo,is)
      end do
      lorbord(ilo,is)=ko
    end do
  end if
  close(50)
end do
! set all muffin-tin radii to single value if required
if (rmtall.gt.0.d0) rmt(1:nspecies)=rmtall
! add conduction state local-orbitals if required
call addlorbcnd
! subtract 2 Hartree from the minimum energy
e0min=e0min-2.d0
return
end subroutine

