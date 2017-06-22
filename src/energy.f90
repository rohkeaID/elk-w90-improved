
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: energy
! !INTERFACE:
subroutine energy
! !USES:
use modmain
use moddftu
use modtest
! !DESCRIPTION:
!   Computes the total energy and its individual contributions. The kinetic
!   energy is given by
!   $$ T_s=\sum_i n_i\epsilon_i-\int\rho({\bf r})[v_{\rm C}({\bf r})
!    +v_{\rm xc}({\bf r})]d{\bf r}-\int {\bf m}({\bf r})\cdot
!    ({\bf B}_{\rm xc}({\bf r})+{\bf B}_{\rm ext}({\bf r}))d{\bf r}, $$
!   where $n_i$ are the occupancies and $\epsilon_i$ are the eigenvalues of both
!   the core and valence states; $\rho$ is the density; ${\bf m}$ is the
!   magnetisation density; $v_{\rm C}$ is the Coulomb potential; $v_{\rm xc}$
!   and ${\bf B}_{\rm xc}$ are the exchange-correlation potential and magnetic
!   field, respectively; and ${\bf B}_{\rm ext}$ is the external magnetic field.
!   The Hartree, electron-nuclear and nuclear-nuclear electrostatic energies are
!   combined into the Coulomb energy:
!   \begin{align*}
!    E_{\rm C}&=E_{\rm H}+E_{\rm en}+E_{\rm nn} \\
!             &=\frac{1}{2}V_{\rm C}+E_{\rm Mad},
!   \end{align*}
!   where
!   $$ V_{\rm C}=\int\rho({\bf r})v_{\rm C}({\bf r})d{\bf r} $$
!   is the Coulomb potential energy. The Madelung energy is given by
!   $$ E_{\rm Mad}=\frac{1}{2}\sum_{\alpha}z_{\alpha}R_{\alpha}, $$
!   where
!   $$ R_{\alpha}=\lim_{r\rightarrow 0}\left(v^{\rm C}_{\alpha;00}(r)Y_{00}
!    +\frac{z_{\alpha}}{r}\right) $$
!   for atom $\alpha$, with $v^{\rm C}_{\alpha;00}$ being the $l=0$ component of
!   the spherical harmonic expansion of $v_{\rm C}$ in the muffin-tin, and
!   $z_{\alpha}$ is the nuclear charge. Using the nuclear-nuclear energy
!   determined at the start of the calculation, the electron-nuclear and Hartree
!   energies can be isolated with
!   $$ E_{\rm en}=2\left(E_{\rm Mad}-E_{\rm nn}\right) $$
!   and
!   $$ E_{\rm H}=\frac{1}{2}(E_{\rm C}-E_{\rm en}). $$
!   Finally, the total energy is
!   $$ E=T_s+E_{\rm C}+E_{\rm xc}, $$
!   where $E_{\rm xc}$ is obtained either by integrating the
!   exchange-correlation energy density, or in the case of exact exchange, the
!   explicit calculation of the Fock exchange integral. The energy from the
!   external magnetic fields in the muffin-tins, {\tt bfcmt}, is always removed
!   from the total since these fields are non-physical: their field lines do not
!   close. The energy of the physical external field, {\tt bfieldc}, is also not
!   included in the total because this field, like those in the muffin-tins, is
!   used for breaking spin symmetry and taken to be infintesimal. If this field
!   is intended to be finite, then the associated energy, {\tt engybext}, should
!   be added to the total by hand. See {\tt potxc}, {\tt exxengy} and related
!   subroutines.
!
! !REVISION HISTORY:
!   Created May 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer ik,ist,idm,jdm,n2,i
integer is,ia,ias,nrc,nrci
real(8) cb,vn,sum,f
complex(8) z1
! allocatable arrays
real(8), allocatable :: rfmt(:,:,:)
complex(8), allocatable :: evecsv(:,:),kmat(:,:),c(:,:)
! external functions
real(8) rfinp
complex(8) zdotc
external rfinp,zdotc
! coupling constant of the external field (g_e/4c)
cb=gfacte/(4.d0*solsc)
!-----------------------------------------------!
!     exchange-correlation potential energy     !
!-----------------------------------------------!
engyvxc=rfinp(1,rhomt,rhoir,vxcmt,vxcir)
!-----------------------------------------------------!
!     exchange-correlation effective field energy     !
!-----------------------------------------------------!
engybxc=0.d0
do idm=1,ndmag
  engybxc=engybxc+rfinp(1,magmt(:,:,:,idm),magir(:,idm),bxcmt(:,:,:,idm), &
   bxcir(:,idm))
end do
!------------------------------------------!
!     external magnetic field energies     !
!------------------------------------------!
engybext=0.d0
do idm=1,ndmag
  if (ncmag) then
    jdm=idm
  else
    jdm=3
  end if
! energy of physical global field
  engybext=engybext+cb*momtot(idm)*bfieldc(jdm)
end do
!----------------------------------!
!     Coulomb potential energy     !
!----------------------------------!
engyvcl=rfinp(1,rhomt,rhoir,vclmt,vclir)
!-----------------------!
!     Madelung term     !
!-----------------------!
engymad=0.d0
do is=1,nspecies
! compute the bare nucleus potential at the origin
  call potnucl(ptnucl,1,rsp(:,is),spzn(is),vn)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    engymad=engymad+0.5d0*spzn(is)*(vclmt(1,1,ias)*y00-vn)
  end do
end do
!---------------------------------------------!
!     electron-nuclear interaction energy     !
!---------------------------------------------!
engyen=2.d0*(engymad-engynn)
!------------------------!
!     Hartree energy     !
!------------------------!
engyhar=0.5d0*(engyvcl-engyen)
!------------------------!
!     Coulomb energy     !
!------------------------!
engycl=engynn+engyen+engyhar
!-------------------------!
!     exchange energy     !
!-------------------------!
if ((xctype(1).lt.0).or.(task.eq.5)) then
! exact exchange for OEP-EXX or Hartree-Fock on last self-consistent loop
  if (tlast) then
    call exxengy
! mix exact and DFT exchange energies for hybrid functionals
    if (hybrid) then
      engyx=engyx*hybridc
      engyx=engyx+rfinp(1,rhomt,rhoir,exmt,exir)
    end if
  else
    engyx=0.d0
  end if
else
! exchange energy from the density
  engyx=rfinp(1,rhomt,rhoir,exmt,exir)
end if
!----------------------------!
!     correlation energy     !
!----------------------------!
if (task.eq.5) then
  if (hybrid) then
! fraction of DFT correlation energy for hybrid functionals
    engyc=rfinp(1,rhomt,rhoir,ecmt,ecir)
  else
! zero correlation energy for pure Hartree-Fock
    engyc=0.d0
  end if
else
! correlation energy from the density
  engyc=rfinp(1,rhomt,rhoir,ecmt,ecir)
end if
!----------------------!
!     DFT+U energy     !
!----------------------!
engydu=0.d0
if (dftu.ne.0) then
  do i=1,ndftu
    is=idftu(1,i)
    do ia=1,natoms(is)
      engydu=engydu+engyadu(ia,i)
    end do
  end do
end if
!----------------------------!
!     sum of eigenvalues     !
!----------------------------!
! core eigenvalues
evalsum=0.d0
do ias=1,natmtot
  is=idxis(ias)
  do ist=1,nstsp(is)
    if (spcore(ist,is)) evalsum=evalsum+occcr(ist,ias)*evalcr(ist,ias)
  end do
end do
! valence eigenvalues
do ik=1,nkpt
  do ist=1,nstsv
    evalsum=evalsum+wkpt(ik)*occsv(ist,ik)*evalsv(ist,ik)
  end do
end do
!------------------------!
!     kinetic energy     !
!------------------------!
! core electron kinetic energy
call energykncr
! total electron kinetic energy
if (task.eq.5) then
! Hartree-Fock case
  engykn=engykncr
! kinetic energy from valence states
  allocate(evecsv(nstsv,nstsv),kmat(nstsv,nstsv),c(nstsv,nstsv))
  do ik=1,nkpt
    call getevecsv(filext,vkl(:,ik),evecsv)
    call getkmat(ik,kmat)
    call zgemm('N','N',nstsv,nstsv,nstsv,zone,kmat,nstsv,evecsv,nstsv,zzero,c, &
     nstsv)
    do ist=1,nstsv
      z1=zdotc(nstsv,evecsv(:,ist),1,c(:,ist),1)
      engykn=engykn+wkpt(ik)*occsv(ist,ik)*dble(z1)
    end do
  end do
  deallocate(evecsv,kmat,c)
else
! Kohn-Sham case
  allocate(rfmt(lmmaxvr,nrmtmax,natmtot))
  sum=0.d0
  do idm=1,ndmag
    do ias=1,natmtot
      is=idxis(ias)
      nrc=nrcmt(is)
      nrci=nrcmtinr(is)
      call rfsht(nrc,nrci,1,bsmt(:,:,ias,idm),lradstp,rfmt(:,:,ias))
    end do
    call rfmtctof(rfmt)
    sum=sum+rfinp(1,magmt(:,:,:,idm),magir(:,idm),rfmt,bsir(:,idm))
  end do
! remove fixed tensor moment potential matrix contribution
  if (ftmtype.ne.0) then
    n2=(lmmaxdm*nspinor)**2
    do ias=1,natmtot
      z1=zdotc(n2,dmatmt(:,:,:,:,ias),1,vmftm(:,:,:,:,ias),1)
      sum=sum+dble(z1)
    end do
  end if
  engykn=evalsum-engyvcl-engyvxc-sum
  deallocate(rfmt)
end if
!-------------------------------!
!     entropic contribution     !
!-------------------------------!
entrpy=0.d0
engyts=0.d0
! non-zero only for the Fermi-Dirac smearing function
if (stype.eq.3) then
  sum=0.d0
  do ik=1,nkpt
    do ist=1,nstsv
      f=occsv(ist,ik)/occmax
      if ((f.gt.0.d0).and.(f.lt.1.d0)) then
        sum=sum+wkpt(ik)*(f*log(f)+(1.d0-f)*log(1.d0-f))
      end if
    end do
  end do
! entropy
  entrpy=-occmax*kboltz*sum
! contribution to free energy
  engyts=-swidth*entrpy/kboltz
end if
!----------------------!
!     total energy     !
!----------------------!
engytot=engykn+0.5d0*engyvcl+engymad+engyx+engyc+engyts
! add the DFT+U correction if required
if (dftu.ne.0) engytot=engytot+engydu
! write total energy to test file
call writetest(0,'total energy',tol=1.d-4,rv=engytot)
return
end subroutine
!EOC

