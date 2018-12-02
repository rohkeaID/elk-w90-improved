
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: readinput
! !INTERFACE:
subroutine readinput
! !USES:
use modmain
use moddftu
use modrdm
use modphonon
use modtest
use modrandom
use modscdft
use modpw
use modtddft
use modulr
use modvars
use modgw
use modjx
use modw90
use modomp
use modmpi
! !DESCRIPTION:
!   Reads in the input parameters from the file {\tt elk.in}. Also sets default
!   values for the input parameters.
!
! !REVISION HISTORY:
!   Created September 2002 (JKD)
!EOP
!BOC
implicit none
! local variables
logical lv
integer is,ia,ias,iostat
integer i,j,k,l,n,p
real(8) sc,sc1,sc2,sc3
real(8) solscf,zn,a,b
real(8) axang(4),rot(3,3)
real(8) v1(3),v2(3),t1
character(256) block,symb,str

!------------------------!
!     default values     !
!------------------------!
ntasks=1
tasks(1)=-1
avec(:,:)=0.d0
avec(1,1)=1.d0
avec(2,2)=1.d0
avec(3,3)=1.d0
sc=1.d0
sc1=1.d0
sc2=1.d0
sc3=1.d0
epslat=1.d-6
primcell=.false.
tshift=.true.
ngridk(:)=1
vkloff(:)=0.d0
autokpt=.false.
radkpt=30.d0
reducek=1
ngridq(:)=-1
reduceq=1
rgkmax=7.d0
gmaxvr=12.d0
lmaxapw=7
lmaxo=7
lmaxi=1
fracinr=0.01d0
trhonorm=.true.
xctype(1)=3
xctype(2:3)=0
xctsp(1)=3
xctsp(2:3)=0
stype=3
swidth=0.001d0
autoswidth=.false.
mstar=10.d0
epsocc=1.d-8
epschg=1.d-3
nempty0=4.d0
maxscl=200
mixtype=3
beta0=0.05d0
betamax=1.d0
! Broyden parameters recommended by M. Meinert
mixsdb=5
broydpm(1)=0.4d0
broydpm(2)=0.15d0
epspot=1.d-6
epsengy=1.d-4
epsforce=5.d-3
epsstress=1.d-3
molecule=.false.
nspecies=0
natoms(:)=0
atposl(:,:,:)=0.d0
atposc(:,:,:)=0.d0
bfcmt0(:,:,:)=0.d0
sppath=''
scrpath=''
nvp1d=2
if (allocated(vvlp1d)) deallocate(vvlp1d)
allocate(vvlp1d(3,nvp1d))
vvlp1d(:,1)=0.d0
vvlp1d(:,2)=1.d0
npp1d=200
vclp2d(:,:)=0.d0
vclp2d(1,1)=1.d0
vclp2d(2,2)=1.d0
np2d(:)=40
vclp3d(:,:)=0.d0
vclp3d(1,1)=1.d0
vclp3d(2,2)=1.d0
vclp3d(3,3)=1.d0
np3d(:)=20
nwplot=500
ngrkf=100
nswplot=1
wplot(1)=-0.5d0
wplot(2)=0.5d0
dosocc=.false.
dosmsum=.false.
dosssum=.false.
lmirep=.true.
spinpol=.false.
spinorb=.false.
socscf=1.d0
atpopt=1
maxatpstp=200
tau0atp=0.25d0
deltast=0.001d0
latvopt=0
maxlatvstp=30
tau0latv=0.2d0
lradstp=4
chgexs=0.d0
scissor=0.d0
noptcomp=1
optcomp(:,1)=1
intraband=.false.
evaltol=-1.d0
epsband=1.d-12
demaxbnd=2.5d0
autolinengy=.false.
dlefe=-0.1d0
deapwlo=0.05d0
bfieldc0(:)=0.d0
efieldc(:)=0.d0
afieldc(:)=0.d0
fsmtype=0
momfix(:)=0.d0
mommtfix(:,:,:)=0.d0
taufsm=0.01d0
rmtdelta=0.05d0
isgkmax=-1
symtype=1
deltaph=0.005d0
nphwrt=1
if (allocated(vqlwrt)) deallocate(vqlwrt)
allocate(vqlwrt(3,nphwrt))
vqlwrt(:,:)=0.d0
notelns=0
tforce=.false.
maxitoep=200
tauoep(1)=1.d0
tauoep(2)=0.75d0
tauoep(3)=1.25d0
nkstlist=1
kstlist(:,1)=1
vklem(:)=0.d0
deltaem=0.025d0
ndspem=1
nosource=.false.
spinsprl=.false.
ssdph=.true.
vqlss(:)=0.d0
nwrite=0
dftu=0
inpdftu=1
ndftu=0
ujdu(:,:)=0.d0
fdu(:,:)=0.d0
edu(:,:)=0.d0
lambdadu(:)=0.d0
udufix(:)=0.d0
lambdadu0(:)=0.d0
tmwrite=.false.
readadu=.false.
rdmxctype=2
rdmmaxscl=2
maxitn=200
maxitc=0
taurdmn=0.5d0
taurdmc=0.25d0
rdmalpha=0.565d0
rdmtemp=0.d0
reducebf=1.d0
ptnucl=.true.
tefvr=.true.
tefvit=.false.
minitefv=6
maxitefv=4
befvit=0.25d0
epsefvit=1.d-5
vecql(:)=0.d0
mustar=0.15d0
sqados(1:2)=0.d0
sqados(3)=1.d0
test=.false.
spincore=.false.
solscf=1.d0
emaxelnes=-1.2d0
wsfac(1)=-1.d6; wsfac(2)=1.d6
vhmat(:,:)=0.d0
vhmat(1,1)=1.d0
vhmat(2,2)=1.d0
vhmat(3,3)=1.d0
reduceh=.true.
hybrid=.false.
hybridc=1.d0
ecvcut=-3.5d0
esccut=-0.4d0
gmaxrf=3.d0
emaxrf=1.d6
ntemp=40
taubdg=0.1d0
nvbse0=2
ncbse0=3
nvxbse=0
ncxbse=0
bsefull=.false.
hxbse=.true.
hdbse=.true.
fxctype=-1
fxclrc(1)=0.d0
fxclrc(2)=0.d0
rndatposc=0.d0
rndbfcmt=0.d0
rndavec=0.d0
emaxbdg=0.5d0
c_tb09=0.d0
tc_tb09=.false.
rndachi=0.1d0
hmaxvr=20.d0
hkmax=12.d0
lorbcnd=.false.
lorbordc=3
nrmtscf=1.d0
lmaxdos=3
epsdev=0.005d0
msmooth=0
npmae0=-1
wrtvars=.false.
ftmtype=0
ntmfix=0
tauftm=0.5d0
ftmstep=1
cmagz=.false.
axang(:)=0.d0
dncgga=1.d-8
tstime=1000.d0
dtimes=0.1d0
npulse=0
nramp=0
ntswrite=10
nxoapwlo=0
nxlo=0
tdrho1d=.false.
tdrho2d=.false.
tdrho3d=.false.
tdmag2d=.false.
tdmag3d=.false.
tdcd2d=.false.
tdcd3d=.false.
rndevt0=0.d0
ssxc=1.d0
avecu(:,:)=0.d0
avecu(1,1)=1.d0
avecu(2,2)=1.d0
avecu(3,3)=1.d0
kpamax=0.d0
evtype=1
wmaxgw=-5.d0
gwdiag=1
actype=10
npole=3
nspade=100
tfav0=.true.
nqssjx=1
thssjx(1)=0.d0
thssjx(2)=pi/6.d0
rmtall=-1.d0
maxthd=0
maxlvl=4
emaxulr=1.d6
taudft=.false.
t0tdlr=0.d0
wann_nwf=-1
wann_nband=-1
wann_projlines=-1
wann_seedname=''
wann_numiter=500
wann_inputlines=0
wann_ngridk(:)=1

!--------------------------!
!     read from elk.in     !
!--------------------------!
open(50,file='elk.in',status='OLD',form='FORMATTED',iostat=iostat)
if (iostat.ne.0) then
  write(*,*)
  write(*,'("Error(readinput): error opening elk.in")')
  write(*,*)
  stop
end if
10 continue
read(50,*,end=30) block
! check for a comment
if ((scan(trim(block),'!').eq.1).or.(scan(trim(block),'#').eq.1)) goto 10
select case(trim(block))
case('tasks')
  do i=1,maxtasks
    read(50,'(A256)',err=20) str
    if (trim(str).eq.'') then
      if (i.eq.1) then
        write(*,*)
        write(*,'("Error(readinput): no tasks to perform")')
        write(*,*)
        stop
      end if
      ntasks=i-1
      goto 10
    end if
    read(str,*,iostat=iostat) tasks(i)
    if (iostat.ne.0) then
      write(*,*)
      write(*,'("Error(readinput): error reading tasks")')
      write(*,'("(blank line required after tasks block)")')
      write(*,*)
      stop
    end if
  end do
  write(*,*)
  write(*,'("Error(readinput): too many tasks")')
  write(*,'("Adjust maxtasks in modmain and recompile code")')
  write(*,*)
  stop
case('species')
! generate a species file
  call genspecies(50)
case('fspecies')
! generate fractional species files
  do is=1,maxspecies
    read(50,'(A256)',err=20) str
    if (trim(str).eq.'') goto 10
    read(str,*,iostat=iostat) zn,symb
    if (zn.gt.-1.d0+epsocc) then
      write(*,*)
      write(*,'("Error(readinput): fractional nuclear Z > -1 : ",G18.10)') zn
      write(*,*)
      stop
    end if
    call genfspecies(zn,symb)
  end do
  write(*,*)
  write(*,'("Error(readinput): too many fractional nucleus species")')
  write(*,*)
  stop
case('avec')
  read(50,*,err=20) avec(:,1)
  read(50,*,err=20) avec(:,2)
  read(50,*,err=20) avec(:,3)
case('scale')
  read(50,*,err=20) sc
case('scale1')
  read(50,*,err=20) sc1
case('scale2')
  read(50,*,err=20) sc2
case('scale3')
  read(50,*,err=20) sc3
case('epslat')
  read(50,*,err=20) epslat
  if (epslat.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): epslat <= 0 : ",G18.10)') epslat
    write(*,*)
    stop
  end if
case('primcell')
  read(50,*,err=20) primcell
case('tshift')
  read(50,*,err=20) tshift
case('autokpt')
  read(50,*,err=20) autokpt
case('radkpt')
  read(50,*,err=20) radkpt
  if (radkpt.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): radkpt <= 0 : ",G18.10)') radkpt
    write(*,*)
    stop
  end if
case('ngridk')
  read(50,*,err=20) ngridk(:)
  if ((ngridk(1).le.0).or.(ngridk(2).le.0).or.(ngridk(3).le.0)) then
    write(*,*)
    write(*,'("Error(readinput): invalid ngridk : ",3I8)') ngridk
    write(*,*)
    stop
  end if
  autokpt=.false.
case('vkloff')
  read(50,*,err=20) vkloff(:)
case('reducek')
  read(50,*,err=20) reducek
case('ngridq')
  read(50,*,err=20) ngridq(:)
  if ((ngridq(1).le.0).or.(ngridq(2).le.0).or.(ngridq(3).le.0)) then
    write(*,*)
    write(*,'("Error(readinput): invalid ngridq : ",3I8)') ngridq
    write(*,*)
    stop
  end if
case('reduceq')
  read(50,*,err=20) reduceq
case('rgkmax')
  read(50,*,err=20) rgkmax
  if (rgkmax.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): rgkmax <= 0 : ",G18.10)') rgkmax
    write(*,*)
    stop
  end if
case('gmaxvr')
  read(50,*,err=20) gmaxvr
case('lmaxapw')
  read(50,*,err=20) lmaxapw
  if (lmaxapw.lt.0) then
    write(*,*)
    write(*,'("Error(readinput): lmaxapw < 0 : ",I8)') lmaxapw
    write(*,*)
    stop
  end if
  if (lmaxapw.ge.maxlapw) then
    write(*,*)
    write(*,'("Error(readinput): lmaxapw too large : ",I8)') lmaxapw
    write(*,'("Adjust maxlapw in modmain and recompile code")')
    write(*,*)
    stop
  end if
case('lmaxo','lmaxvr')
  read(50,*,err=20) lmaxo
  if (lmaxo.lt.3) then
    write(*,*)
    write(*,'("Error(readinput): lmaxo < 3 : ",I8)') lmaxo
    write(*,*)
    stop
  end if
case('lmaxi','lmaxinr')
  read(50,*,err=20) lmaxi
  if (lmaxi.lt.0) then
    write(*,*)
    write(*,'("Error(readinput): lmaxi < 0 : ",I8)') lmaxi
    write(*,*)
    stop
  end if
case('lmaxmat')
  read(50,*,err=20)
  write(*,*)
  write(*,'("Info(readinput): variable ''lmaxmat'' is no longer used")')
case('fracinr')
  read(50,*,err=20) fracinr
case('trhonorm')
  read(50,*,err=20) trhonorm
case('spinpol')
  read(50,*,err=20) spinpol
case('spinorb')
  read(50,*,err=20) spinorb
case('socscf')
  read(50,*,err=20) socscf
  if (socscf.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): socscf < 0 : ",G18.10)') socscf
    write(*,*)
    stop
  end if
case('xctype')
  read(50,'(A256)',err=20) str
  str=trim(str)//' 0 0'
  read(str,*,err=20) xctype
case('xctsp')
  read(50,'(A256)',err=20) str
  str=trim(str)//' 0 0'
  read(str,*,err=20) xctsp
case('stype')
  read(50,*,err=20) stype
case('swidth')
  read(50,*,err=20) swidth
  if (swidth.lt.1.d-9) then
    write(*,*)
    write(*,'("Error(readinput): swidth too small or negative : ",G18.10)') &
     swidth
    write(*,*)
    stop
  end if
case('autoswidth')
  read(50,*,err=20) autoswidth
case('mstar')
  read(50,*,err=20) mstar
  if (mstar.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): mstar <= 0 : ",G18.10)') mstar
    write(*,*)
    stop
  end if
case('epsocc')
  read(50,*,err=20) epsocc
  if (epsocc.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): epsocc <= 0 : ",G18.10)') epsocc
    write(*,*)
    stop
  end if
case('epschg')
  read(50,*,err=20) epschg
  if (epschg.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): epschg <= 0 : ",G18.10)') epschg
    write(*,*)
    stop
  end if
case('nempty','nempty0')
  read(50,*,err=20) nempty0
  if (nempty0.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): nempty <= 0 : ",G18.10)') nempty0
    write(*,*)
    stop
  end if
case('mixtype')
  read(50,*,err=20) mixtype
case('beta0')
  read(50,*,err=20) beta0
  if (beta0.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): beta0 < 0 : ",G18.10)') beta0
    write(*,*)
    stop
  end if
case('betamax')
  read(50,*,err=20) betamax
  if ((betamax.lt.0.d0).or.(betamax.gt.1.d0)) then
    write(*,*)
    write(*,'("Error(readinput): betmax not in [0,1] : ",G18.10)') betamax
    write(*,*)
    stop
  end if
case('mixsdb')
  read(50,*,err=20) mixsdb
  if (mixsdb.lt.2) then
    write(*,*)
    write(*,'("Error(readinput): mixsdb < 2 : ",I8)') mixsdb
    write(*,*)
    stop
  end if
case('broydpm')
  read(50,*,err=20) broydpm(:)
  if ((broydpm(1).lt.0.d0).or.(broydpm(1).gt.1.d0).or. &
      (broydpm(2).lt.0.d0).or.(broydpm(2).gt.1.d0)) then
    write(*,*)
    write(*,'("Error(readinput): invalid Broyden mixing parameters : ",&
     &2G18.10)') broydpm
    write(*,*)
    stop
  end if
case('maxscl')
  read(50,*,err=20) maxscl
  if (maxscl.lt.0) then
    write(*,*)
    write(*,'("Error(readinput): maxscl < 0 : ",I8)') maxscl
    write(*,*)
    stop
  end if
case('epspot')
  read(50,*,err=20) epspot
case('epsengy')
  read(50,*,err=20) epsengy
case('epsforce')
  read(50,*,err=20) epsforce
case('epsstress')
  read(50,*,err=20) epsstress
case('sppath')
  read(50,*,err=20) sppath
  sppath=adjustl(sppath)
case('scrpath')
  read(50,*,err=20) scrpath
case('molecule')
  read(50,*,err=20) molecule
case('atoms')
  read(50,*,err=20) nspecies
  if (nspecies.le.0) then
    write(*,*)
    write(*,'("Error(readinput): nspecies <= 0 : ",I8)') nspecies
    write(*,*)
    stop
  end if
  if (nspecies.gt.maxspecies) then
    write(*,*)
    write(*,'("Error(readinput): nspecies too large : ",I8)') nspecies
    write(*,'("Adjust maxspecies in modmain and recompile code")')
    write(*,*)
    stop
  end if
  do is=1,nspecies
    read(50,*,err=20) spfname(is)
    spfname(is)=adjustl(spfname(is))
    read(50,*,err=20) natoms(is)
    if (natoms(is).le.0) then
      write(*,*)
      write(*,'("Error(readinput): natoms <= 0 : ",I8)') natoms(is)
      write(*,'(" for species ",I4)') is
      write(*,*)
      stop
    end if
    if (natoms(is).gt.maxatoms) then
      write(*,*)
      write(*,'("Error(readinput): natoms too large : ",I8)') natoms(is)
      write(*,'(" for species ",I4)') is
      write(*,'("Adjust maxatoms in modmain and recompile code")')
      write(*,*)
      stop
    end if
    do ia=1,natoms(is)
      read(50,'(A256)',err=20) str
      str=trim(str)//' 0.0 0.0 0.0'
      read(str,*,err=20) atposl(:,ia,is),bfcmt0(:,ia,is)
    end do
  end do
case('plot1d')
  read(50,*,err=20) nvp1d,npp1d
  if (nvp1d.lt.1) then
    write(*,*)
    write(*,'("Error(readinput): nvp1d < 1 : ",I8)') nvp1d
    write(*,*)
    stop
  end if
  if (npp1d.lt.nvp1d) then
    write(*,*)
    write(*,'("Error(readinput): npp1d < nvp1d : ",2I8)') npp1d,nvp1d
    write(*,*)
    stop
  end if
  if (allocated(vvlp1d)) deallocate(vvlp1d)
  allocate(vvlp1d(3,nvp1d))
  do i=1,nvp1d
    read(50,*,err=20) vvlp1d(:,i)
  end do
case('plot2d')
  read(50,*,err=20) vclp2d(:,0)
  read(50,*,err=20) vclp2d(:,1)
  read(50,*,err=20) vclp2d(:,2)
  read(50,*,err=20) np2d(:)
  if ((np2d(1).lt.1).or.(np2d(2).lt.1)) then
    write(*,*)
    write(*,'("Error(readinput): np2d < 1 : ",2I8)') np2d
    write(*,*)
    stop
  end if
case('plot3d')
  read(50,*,err=20) vclp3d(:,0)
  read(50,*,err=20) vclp3d(:,1)
  read(50,*,err=20) vclp3d(:,2)
  read(50,*,err=20) vclp3d(:,3)
  read(50,*,err=20) np3d(:)
  if ((np3d(1).lt.1).or.(np3d(2).lt.1).or.(np3d(3).lt.1)) then
    write(*,*)
    write(*,'("Error(readinput): np3d < 1 : ",3I8)') np3d
    write(*,*)
    stop
  end if
case('wplot','dos')
  read(50,*,err=20) nwplot,ngrkf,nswplot
  if (nwplot.lt.2) then
    write(*,*)
    write(*,'("Error(readinput): nwplot < 2 : ",I8)') nwplot
    write(*,*)
    stop
  end if
  if (ngrkf.lt.1) then
    write(*,*)
    write(*,'("Error(readinput): ngrkf < 1 : ",I8)') ngrkf
    write(*,*)
    stop
  end if
  if (nswplot.lt.0) then
    write(*,*)
    write(*,'("Error(readinput): nswplot < 0 : ",I8)') nswplot
    write(*,*)
    stop
  end if
  read(50,*,err=20) wplot(:)
  if (wplot(1).gt.wplot(2)) then
    write(*,*)
    write(*,'("Error(readinput): wplot(1) > wplot(2) : ",2G18.10)') wplot
    write(*,*)
    stop
  end if
case('dosocc')
  read(50,*,err=20) dosocc
case('dosmsum')
  read(50,*,err=20) dosmsum
case('dosssum')
  read(50,*,err=20) dosssum
case('lmirep')
  read(50,*,err=20) lmirep
case('atpopt')
  read(50,*,err=20) atpopt
case('maxatpstp','maxatmstp')
  read(50,*,err=20) maxatpstp
  if (maxatpstp.le.0) then
    write(*,*)
    write(*,'("Error(readinput): maxatpstp <= 0 : ",I8)') maxatpstp
    write(*,*)
    stop
  end if
case('tau0atp','tau0atm')
  read(50,*,err=20) tau0atp
case('deltast')
  read(50,*,err=20) deltast
  if (deltast.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): deltast <= 0 : ",G18.10)') deltast
    write(*,*)
    stop
  end if
case('latvopt')
  read(50,*,err=20) latvopt
case('maxlatvstp')
  read(50,*,err=20) maxlatvstp
  if (maxlatvstp.le.0) then
    write(*,*)
    write(*,'("Error(readinput): maxlatvstp <= 0 : ",I8)') maxlatvstp
    write(*,*)
    stop
  end if
case('tau0latv')
  read(50,*,err=20) tau0latv
case('nstfsp')
  read(50,*,err=20)
  write(*,*)
  write(*,'("Info(readinput): variable ''nstfsp'' is no longer used")')
case('lradstp')
  read(50,*,err=20) lradstp
  if (lradstp.le.0) then
    write(*,*)
    write(*,'("Error(readinput): lradstp <= 0 : ",I8)') lradstp
    write(*,*)
    stop
  end if
case('chgexs')
  read(50,*,err=20) chgexs
case('nprad')
  read(50,*,err=20)
  write(*,*)
  write(*,'("Info(readinput): variable ''nprad'' is no longer used")')
case('scissor')
  read(50,*,err=20) scissor
case('optcomp')
  do i=1,27
    read(50,'(A256)',err=20) str
    if (trim(str).eq.'') then
      if (i.eq.1) then
        write(*,*)
        write(*,'("Error(readinput): empty optical component list")')
        write(*,*)
        stop
      end if
      noptcomp=i-1
      goto 10
    end if
    str=trim(str)//' 1 1'
    read(str,*,iostat=iostat) optcomp(:,i)
    if (iostat.ne.0) then
      write(*,*)
      write(*,'("Error(readinput): error reading optical component list")')
      write(*,'("(blank line required after optcomp block)")')
      write(*,*)
      stop
    end if
    if ((optcomp(1,i).lt.1).or.(optcomp(1,i).gt.3).or. &
        (optcomp(2,i).lt.1).or.(optcomp(2,i).gt.3).or. &
        (optcomp(3,i).lt.1).or.(optcomp(3,i).gt.3)) then
      write(*,*)
      write(*,'("Error(readinput): invalid optcomp : ",3I8)') optcomp
      write(*,*)
      stop
    end if
  end do
  write(*,*)
  write(*,'("Error(readinput): optical component list too long")')
  write(*,*)
  stop
case('intraband')
  read(50,*,err=20) intraband
case('evaltol')
  read(50,*,err=20) evaltol
case('deband')
  read(50,*,err=20)
  write(*,*)
  write(*,'("Info(readinput): variable ''deband'' is no longer used")')
case('epsband')
  read(50,*,err=20) epsband
  if (epsband.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): epsband <= 0 : ",G18.10)') epsband
    write(*,*)
    stop
  end if
case('demaxbnd')
  read(50,*,err=20) demaxbnd
  if (demaxbnd.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): demaxbnd <= 0 : ",G18.10)') demaxbnd
    write(*,*)
    stop
  end if
case('autolinengy')
  read(50,*,err=20) autolinengy
case('dlefe')
  read(50,*,err=20) dlefe
case('deapwlo')
  read(50,*,err=20) deapwlo
  if (abs(deapwlo).lt.1.d-8) then
    write(*,*)
    write(*,'("Error(readinput): invalid deapwlo : ",G18.10)') deapwlo
    write(*,*)
    stop
  end if
case('bfieldc')
  read(50,*,err=20) bfieldc0
case('efieldc')
  read(50,*,err=20) efieldc
case('afieldc')
  read(50,*,err=20) afieldc
case('fsmtype','fixspin')
  read(50,*,err=20) fsmtype
case('momfix')
  read(50,*,err=20) momfix
case('mommtfix')
  do ias=1,maxspecies*maxatoms
    read(50,'(A256)',err=20) str
    if (trim(str).eq.'') goto 10
    read(str,*,iostat=iostat) is,ia,mommtfix(:,ia,is)
    if (iostat.ne.0) then
      write(*,*)
      write(*,'("Error(readinput): error reading muffin-tin fixed spin &
       &moments")')
      write(*,'("(blank line required after mommtfix block")')
      write(*,*)
      stop
    end if
  end do
case('taufsm')
  read(50,*,err=20) taufsm
  if (taufsm.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): taufsm < 0 : ",G18.10)') taufsm
    write(*,*)
    stop
  end if
case('autormt')
  read(50,*,err=20)
  write(*,*)
  write(*,'("Info(readinput): variable ''autormt'' is no longer used")')
case('rmtdelta')
  read(50,*,err=20) rmtdelta
  if (rmtdelta.lt.0.d0) then
    write(*,*)
    write(*,'("Warning(readinput): rmtdelta < 0 : ",G18.10)') rmtdelta
  end if
case('isgkmax')
  read(50,*,err=20) isgkmax
case('nosym')
  read(50,*,err=20) lv
  if (lv) symtype=0
case('symtype')
  read(50,*,err=20) symtype
  if ((symtype.lt.0).or.(symtype.gt.2)) then
    write(*,*)
    write(*,'("Error(readinput): symtype not defined : ",I8)') symtype
    write(*,*)
    stop
  end if
case('deltaph')
  read(50,*,err=20) deltaph
  if (deltaph.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): deltaph <= 0 : ",G18.10)') deltaph
    write(*,*)
    stop
  end if
case('phwrite')
  read(50,*,err=20) nphwrt
  if (nphwrt.le.0) then
    write(*,*)
    write(*,'("Error(readinput): nphwrt <= 0 : ",I8)') nphwrt
    write(*,*)
    stop
  end if
  if (allocated(vqlwrt)) deallocate(vqlwrt)
  allocate(vqlwrt(3,nphwrt))
  do i=1,nphwrt
    read(50,*,err=20) vqlwrt(:,i)
  end do
case('notes')
  do i=1,maxnlns
    read(50,'(A80)') notes(i)
    if (trim(notes(i)).eq.'') then
      notelns=i-1
      goto 10
    end if
  end do
  write(*,*)
  write(*,'("Error(readinput): too many note lines")')
  write(*,*)
  stop
case('tforce')
  read(50,*,err=20) tforce
case('tfibs')
  read(50,*,err=20)
  write(*,*)
  write(*,'("Info(readinput): variable ''tfibs'' is no longer used")')
case('maxitoep')
  read(50,*,err=20) maxitoep
  if (maxitoep.lt.1) then
    write(*,*)
    write(*,'("Error(readinput): maxitoep < 1 : ",I8)') maxitoep
    write(*,*)
    stop
  end if
case('tauoep')
  read(50,*,err=20) tauoep(:)
  if ((tauoep(1).lt.0.d0).or.(tauoep(2).lt.0.d0).or.(tauoep(3).lt.0.d0)) then
    write(*,*)
    write(*,'("Error(readinput): tauoep < 0 : ",3G18.10)') tauoep
    write(*,*)
    stop
  end if
case('kstlist')
  do i=1,maxkst
    read(50,'(A256)',err=20) str
    if (trim(str).eq.'') then
      if (i.eq.1) then
        write(*,*)
        write(*,'("Error(readinput): empty k-point and state list")')
        write(*,*)
        stop
      end if
      nkstlist=i-1
      goto 10
    end if
    read(str,*,iostat=iostat) kstlist(:,i)
    if (iostat.ne.0) then
      write(*,*)
      write(*,'("Error(readinput): error reading k-point and state list")')
      write(*,'("(blank line required after kstlist block)")')
      write(*,*)
      stop
    end if
  end do
  write(*,*)
  write(*,'("Error(readinput): k-point and state list too long")')
  write(*,*)
  stop
case('vklem')
  read(50,*,err=20) vklem
case('deltaem')
  read(50,*,err=20) deltaem
case('ndspem')
  read(50,*,err=20) ndspem
  if ((ndspem.lt.1).or.(ndspem.gt.4)) then
    write(*,*)
    write(*,'("Error(readinput): ndspem out of range : ",I8)') ndspem
    write(*,*)
    stop
  end if
case('nosource')
  read(50,*,err=20) nosource
case('spinsprl')
  read(50,*,err=20) spinsprl
case('ssdph')
  read(50,*,err=20) ssdph
case('vqlss')
  read(50,*,err=20) vqlss
case('nwrite')
  read(50,*,err=20) nwrite
case('DFT+U','dft+u','lda+u')
  read(50,*,err=20) dftu,inpdftu
  do i=1,maxdftu
    read(50,'(A256)',err=20) str
    if (trim(str).eq.'') then
      ndftu=i-1
      goto 10
    end if
    select case(inpdftu)
    case(1)
      read(str,*,iostat=iostat) is,l,ujdu(1:2,i)
    case(2)
      read(str,*,iostat=iostat) is,l,(fdu(k,i),k=0,2*l,2)
    case(3)
      read(str,*,iostat=iostat) is,l,(edu(k,i),k=0,l)
    case(4)
      read(str,*,iostat=iostat) is,l,lambdadu(i)
    case(5)
      read(str,*,iostat=iostat) is,l,udufix(i)
    case default
      write(*,*)
      write(*,'("Error(readinput): invalid inpdftu : ",I8)') inpdftu
      write(*,*)
      stop
    end select
    if (iostat.ne.0) then
      write(*,*)
      write(*,'("Error(readinput): error reading DFT+U parameters")')
      write(*,'("(blank line required after dft+u block)")')
      write(*,*)
      stop
    end if
    if ((is.le.0).or.(is.ge.maxspecies)) then
      write(*,*)
      write(*,'("Error(readinput): invalid species number in dft+u block : ", &
       &I8)') is
      write(*,*)
      stop
    end if
    if (l.lt.0) then
      write(*,*)
      write(*,'("Error(readinput): l < 0 in dft+u block : ",I8)') l
      write(*,*)
      stop
    end if
    if (l.gt.lmaxdm) then
      write(*,*)
      write(*,'("Error(readinput): l > lmaxdm in dft+u block : ",2I8)') l,lmaxdm
      write(*,*)
      stop
    end if
! check for repeated entries
    do j=1,i-1
      if ((is.eq.idftu(1,j)).and.(l.eq.idftu(2,j))) then
        write(*,*)
        write(*,'("Error(readinput): repeated entry in DFT+U block")')
        write(*,*)
        stop
      end if
    end do
    idftu(1,i)=is
    idftu(2,i)=l
  end do
  write(*,*)
  write(*,'("Error(readinput): too many DFT+U entries")')
  write(*,'("Adjust maxdftu in modmain and recompile code")')
  write(*,*)
  stop
case('tmwrite','tmomlu')
  read(50,*,err=20) tmwrite
case('readadu','readalu')
  read(50,*,err=20) readadu
case('rdmxctype')
  read(50,*,err=20) rdmxctype
case('rdmmaxscl')
  read(50,*,err=20) rdmmaxscl
  if (rdmmaxscl.lt.0) then
    write(*,*)
    write(*,'("Error(readinput): rdmmaxscl < 0 : ",I8)') rdmmaxscl
    write(*,*)
  end if
case('maxitn')
  read(50,*,err=20) maxitn
  if (maxitn.lt.1) then
    write(*,*)
    write(*,'("Error(readinput): maxitn < 1 : ",I8)') maxitn
    write(*,*)
    stop
  end if
case('maxitc')
  read(50,*,err=20) maxitc
case('taurdmn')
  read(50,*,err=20) taurdmn
  if (taurdmn.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): taurdmn < 0 : ",G18.10)') taurdmn
    write(*,*)
    stop
  end if
case('taurdmc')
  read(50,*,err=20) taurdmc
  if (taurdmc.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): taurdmc < 0 : ",G18.10)') taurdmc
    write(*,*)
    stop
  end if
case('rdmalpha')
  read(50,*,err=20) rdmalpha
  if ((rdmalpha.le.0.d0).or.(rdmalpha.ge.1.d0)) then
    write(*,*)
    write(*,'("Error(readinput): rdmalpha not in (0,1) : ",G18.10)') rdmalpha
    write(*,*)
    stop
  end if
case('rdmtemp')
  read(50,*,err=20) rdmtemp
  if (rdmtemp.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): rdmtemp < 0 : ",G18.10)') rdmtemp
    write(*,*)
    stop
  end if
case('reducebf')
  read(50,*,err=20) reducebf
  if ((reducebf.lt.0.49d0).or.(reducebf.gt.1.d0)) then
    write(*,*)
    write(*,'("Error(readinput): reducebf not in [0.5,1] : ",G18.10)') reducebf
    write(*,*)
    stop
  end if
case('ptnucl')
  read(50,*,err=20) ptnucl
case('tefvr','tseqr')
  read(50,*,err=20) tefvr
case('tefvit','tseqit')
  read(50,*,err=20) tefvit
case('minitefv','minseqit')
  read(50,*,err=20) minitefv
  if (minitefv.lt.1) then
    write(*,*)
    write(*,'("Error(readinput): minitefv < 1 : ",I8)') minitefv
    write(*,*)
    stop
  end if
case('maxitefv','maxseqit')
  read(50,*,err=20) maxitefv
  if (maxitefv.lt.1) then
    write(*,*)
    write(*,'("Error(readinput): maxitefv < 1 : ",I8)') maxitefv
    write(*,*)
    stop
  end if
case('befvit','bseqit')
  read(50,*,err=20) befvit
  if (befvit.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): befvit <= 0 : ",G18.10)') befvit
    write(*,*)
    stop
  end if
case('epsefvit','epsseqit')
  read(50,*,err=20) epsefvit
  if (epsefvit.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): epsefvit < 0 : ",G18.10)') epsefvit
    write(*,*)
    stop
  end if
case('nseqit')
  read(50,*,err=20)
  write(*,*)
  write(*,'("Info(readinput): variable ''nseqit'' is no longer used")')
case('tauseq')
  read(50,*,err=20)
  write(*,*)
  write(*,'("Info(readinput): variable ''tauseq'' is no longer used")')
case('vecql')
  read(50,*,err=20) vecql(:)
case('mustar')
  read(50,*,err=20) mustar
case('sqados')
  read(50,*,err=20) sqados(:)
case('test')
  read(50,*,err=20) test
case('frozencr')
  read(50,*,err=20)
  write(*,*)
  write(*,'("Info(readinput): variable ''frozencr'' is no longer used")')
case('spincore')
  read(50,*,err=20) spincore
case('solscf')
  read(50,*,err=20) solscf
  if (solscf.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): solscf < 0 : ",G18.10)') solscf
    write(*,*)
    stop
  end if
case('emaxelnes')
  read(50,*,err=20) emaxelnes
case('wsfac')
  read(50,*,err=20) wsfac(:)
case('vhmat')
  read(50,*,err=20) vhmat(1,:)
  read(50,*,err=20) vhmat(2,:)
  read(50,*,err=20) vhmat(3,:)
case('reduceh')
  read(50,*,err=20) reduceh
case('hybrid')
  read(50,*,err=20) hybrid
case('hybridc','hybmix')
  read(50,*,err=20) hybridc
  if ((hybridc.lt.0.d0).or.(hybridc.gt.1.d0)) then
    write(*,*)
    write(*,'("Error(readinput): invalid hybridc : ",G18.10)') hybridc
    write(*,*)
    stop
  end if
case('ecvcut')
  read(50,*,err=20) ecvcut
case('esccut')
  read(50,*,err=20) esccut
case('nvbse')
  read(50,*,err=20) nvbse0
  if (nvbse0.lt.0) then
    write(*,*)
    write(*,'("Error(readinput): nvbse < 0 : ",I8)') nvbse0
    write(*,*)
    stop
  end if
case('ncbse')
  read(50,*,err=20) ncbse0
  if (ncbse0.lt.0) then
    write(*,*)
    write(*,'("Error(readinput): ncbse < 0 : ",I8)') ncbse0
    write(*,*)
    stop
  end if
case('istxbse')
  do i=1,maxxbse
    read(50,'(A256)',err=20) str
    if (trim(str).eq.'') then
      if (i.eq.1) then
        write(*,*)
        write(*,'("Error(readinput): empty BSE extra valence state list")')
        write(*,*)
        stop
      end if
      nvxbse=i-1
      goto 10
    end if
    read(str,*,iostat=iostat) istxbse(i)
    if (iostat.ne.0) then
      write(*,*)
      write(*,'("Error(readinput): error reading BSE valence state list")')
      write(*,'("(blank line required after istxbse block)")')
      write(*,*)
      stop
    end if
  end do
  write(*,*)
  write(*,'("Error(readinput): BSE extra valence state list too long")')
  write(*,*)
  stop
case('jstxbse')
  do i=1,maxxbse
    read(50,'(A256)',err=20) str
    if (trim(str).eq.'') then
      if (i.eq.1) then
        write(*,*)
        write(*,'("Error(readinput): empty BSE extra conduction state list")')
        write(*,*)
        stop
      end if
      ncxbse=i-1
      goto 10
    end if
    read(str,*,iostat=iostat) jstxbse(i)
    if (iostat.ne.0) then
      write(*,*)
      write(*,'("Error(readinput): error reading BSE conduction state list")')
      write(*,'("(blank line required after jstxbse block)")')
      write(*,*)
      stop
    end if
  end do
  write(*,*)
  write(*,'("Error(readinput): BSE extra conduction state list too long")')
  write(*,*)
  stop
case('bsefull')
  read(50,*,err=20) bsefull
case('hxbse')
  read(50,*,err=20) hxbse
case('hdbse')
  read(50,*,err=20) hdbse
case('gmaxrf','gmaxrpa')
  read(50,*,err=20) gmaxrf
  if (gmaxrf.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): gmaxrf < 0 : ",G18.10)') gmaxrf
    write(*,*)
    stop
  end if
case('emaxrf')
  read(50,*,err=20) emaxrf
  if (emaxrf.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): emaxrf < 0 : ",G18.10)') emaxrf
    write(*,*)
    stop
  end if
case('fxctype')
  read(50,'(A256)',err=20) str
  str=trim(str)//' 0 0'
  read(str,*,err=20) fxctype
case('fxclrc')
  read(50,'(A256)',err=20) str
  str=trim(str)//' 0.0'
  read(str,*,err=20) fxclrc(:)
case('ntemp')
  read(50,*,err=20) ntemp
  if (ntemp.lt.1) then
    write(*,*)
    write(*,'("Error(readinput): ntemp < 1 : ",I8)') ntemp
    write(*,*)
    stop
  end if
case('trimvg')
  read(50,*,err=20)
  write(*,*)
  write(*,'("Info(readinput): variable ''trimvg'' is no longer used")')
  write(*,'(" setting msmooth=4 instead")')
  msmooth=4
case('rndstate','rndseed')
  read(50,*,err=20) rndstate(0)
  rndstate(0)=abs(rndstate(0))
case('taubdg')
  read(50,*,err=20) taubdg
case('rndatposc')
  read(50,*,err=20) rndatposc
case('rndbfcmt')
  read(50,*,err=20) rndbfcmt
case('rndavec')
  read(50,*,err=20) rndavec
case('emaxbdg')
  read(50,*,err=20) emaxbdg
  if (emaxbdg.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): emaxbdg <= 0 : ",G18.10)') emaxbdg
    write(*,*)
    stop
  end if
case('c_tb09')
  read(50,*,err=20) c_tb09
! set flag to indicate Tran-Blaha constant has been read in
  tc_tb09=.true.
case('rndachi')
  read(50,*,err=20) rndachi
case('lowq','highq','vhighq','uhighq')
  read(50,*,err=20) lv
  if (lv) then
    if (trim(block).eq.'lowq') then
      rgkmax=6.5d0
      gmaxvr=10.d0
      lmaxapw=6
      lmaxo=6
      fracinr=0.02d0
      nxlo=2
      radkpt=25.d0
      vkloff(:)=0.5d0
      nempty0=4.d0
      epspot=1.d-5
      epsengy=5.d-4
      epsforce=1.d-2
      autolinengy=.true.
      gmaxrf=2.5d0
      lradstp=6
    else if (trim(block).eq.'highq') then
! parameter set for high quality calculation
      rgkmax=max(rgkmax,8.d0)
      gmaxvr=max(gmaxvr,16.d0)
      msmooth=max(msmooth,4)
      lmaxapw=max(lmaxapw,8)
      lmaxo=max(lmaxo,8)
      fracinr=min(fracinr,0.005d0)
      nrmtscf=max(nrmtscf,1.5d0)
      nxlo=max(nxlo,2)
      radkpt=max(radkpt,50.d0)
      autokpt=.true.
      vkloff(:)=0.d0
      nempty0=max(nempty0,10.d0)
      epspot=min(epspot,1.d-7)
      epsengy=min(epsengy,1.d-5)
      epsforce=min(epsforce,1.d-4)
      autolinengy=.true.
      gmaxrf=max(gmaxrf,4.d0)
    else if (trim(block).eq.'vhighq') then
! parameter set for very high quality calculation
      rgkmax=max(rgkmax,9.d0)
      gmaxvr=max(gmaxvr,18.d0)
      msmooth=max(msmooth,4)
      lmaxapw=max(lmaxapw,12)
      lmaxo=max(lmaxo,10)
      fracinr=min(fracinr,0.002d0)
      nrmtscf=max(nrmtscf,2.d0)
      nxlo=max(nxlo,3)
      radkpt=max(radkpt,90.d0)
      autokpt=.true.
      vkloff(:)=0.d0
      nempty0=max(nempty0,20.d0)
      epspot=min(epspot,1.d-7)
      epsengy=min(epsengy,1.d-6)
      epsforce=min(epsforce,1.d-4)
      autolinengy=.true.
      gmaxrf=max(gmaxrf,5.d0)
    else
! parameter set for ultra high quality calculation
      rgkmax=max(rgkmax,10.d0)
      gmaxvr=max(gmaxvr,20.d0)
      msmooth=max(msmooth,4)
      lmaxapw=max(lmaxapw,12)
      lmaxo=max(lmaxo,12)
      fracinr=min(fracinr,0.001d0)
      nrmtscf=max(nrmtscf,4.d0)
      nxlo=max(nxlo,3)
      lorbcnd=.true.
      radkpt=max(radkpt,120.d0)
      autokpt=.true.
      vkloff(:)=0.d0
      nempty0=max(nempty0,40.d0)
      epspot=min(epspot,1.d-7)
      epsengy=min(epsengy,1.d-6)
      epsforce=min(epsforce,1.d-4)
      autolinengy=.true.
      gmaxrf=max(gmaxrf,6.d0)
    end if
    if (mp_mpi) then
      write(*,*)
      write(*,'("Info(readinput): parameters set by ",A," option")') trim(block)
      write(*,'(" rgkmax : ",G18.10)') rgkmax
      write(*,'(" gmaxvr : ",G18.10)') gmaxvr
      write(*,'(" msmooth : ",I4)') msmooth
      write(*,'(" lmaxapw : ",I4)') lmaxapw
      write(*,'(" lmaxo : ",I4)') lmaxo
      write(*,'(" fracinr : ",G18.10)') fracinr
      write(*,'(" nrmtscf : ",G18.10)') nrmtscf
      write(*,'(" nxlo : ",I4)') nxlo
      write(*,'(" radkpt : ",G18.10)') radkpt
      write(*,'(" autokpt : ",L1)') autokpt
      write(*,'(" vkloff : ",3G18.10)') vkloff
      write(*,'(" nempty0 : ",G18.10)') nempty0
      write(*,'(" epspot : ",G18.10)') epspot
      write(*,'(" epsengy : ",G18.10)') epsengy
      write(*,'(" epsforce : ",G18.10)') epsforce
      write(*,'(" autolinengy : ",L1)') autolinengy
      write(*,'(" gmaxrf : ",G18.10)') gmaxrf
      if (trim(block).eq.'uhighq') then
        write(*,'(" lorbcnd : ",L1)') lorbcnd
      end if
      if (trim(block).eq.'lowq') then
        write(*,'(" lradstp : ",I4)') lradstp
      end if
    end if
  end if
case('hmaxvr')
  read(50,*,err=20) hmaxvr
  if (hmaxvr.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): hmaxvr < 0 : ",G18.10)') hmaxvr
    write(*,*)
    stop
  end if
case('hkmax')
  read(50,*,err=20) hkmax
  if (hkmax.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): hkmax <= 0 : ",G18.10)') hkmax
    write(*,*)
    stop
  end if
case('lorbcnd')
  read(50,*,err=20) lorbcnd
case('lorbordc')
  read(50,*,err=20) lorbordc
  if (lorbordc.lt.2) then
    write(*,*)
    write(*,'("Error(readinput): lorbordc < 2 : ",I8)') lorbordc
    write(*,*)
    stop
  end if
  if (lorbordc.gt.maxlorbord) then
    write(*,*)
    write(*,'("Error(readinput): lorbordc too large : ",I8)') lorbordc
    write(*,'("Adjust maxlorbord in modmain and recompile code")')
    write(*,*)
    stop
  end if
case('nrmtscf')
  read(50,*,err=20) nrmtscf
  if (nrmtscf.lt.0.5d0) then
    write(*,*)
    write(*,'("Error(readinput): nrmtscf < 0.5 : ",G18.10)') nrmtscf
    write(*,*)
    stop
  end if
case('lmaxdos')
  read(50,*,err=20) lmaxdos
  if (lmaxdos.lt.0) then
    write(*,*)
    write(*,'("Error(readinput): lmaxdos < 0 : ",I8)') lmaxdos
    write(*,*)
    stop
  end if
case('epsdev')
  read(50,*,err=20) epsdev
  if (epsdev.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): epsdev <= 0 : ",G18.10)') epsdev
    write(*,*)
    stop
  end if
case('msmooth')
  read(50,*,err=20) msmooth
  if (msmooth.lt.0) then
    write(*,*)
    write(*,'("Error(readinput): msmooth < 0 : ",I8)') msmooth
    write(*,*)
    stop
  end if
case('npmae')
  read(50,*,err=20) npmae0
case('wrtvars')
  read(50,*,err=20) wrtvars
case('ftmtype')
  read(50,*,err=20) ftmtype
case('tmomfix')
  read(50,*,err=20) ntmfix
  if (ntmfix.le.0) then
    write(*,*)
    write(*,'("Error(readinput): ntmfix <= 0 : ",I8)') ntmfix
    write(*,*)
    stop
  end if
  if (allocated(itmfix)) deallocate(itmfix)
  allocate(itmfix(8,ntmfix))
  if (allocated(tmfix)) deallocate(tmfix)
  allocate(tmfix(ntmfix))
  if (allocated(rtmfix)) deallocate(rtmfix)
  allocate(rtmfix(3,3,2,ntmfix))
  do i=1,ntmfix
    read(50,*,err=20) is,ia,l,n
    if ((is.le.0).or.(ia.le.0).or.(l.lt.0).or.((n.ne.2).and.(n.ne.3))) then
      write(*,*)
      write(*,'("Error(readinput): invalid is, ia, l or n in tmomfix block : ",&
       &4I8)') is,ia,l,n
      write(*,*)
      stop
    end if
    itmfix(1,i)=is
    itmfix(2,i)=ia
    itmfix(3,i)=l
    itmfix(4,i)=n
! read k, p, x, y for the 2-index tensor or k, p, r, t for the 3-index tensor
    read(50,*,err=20) itmfix(5:8,i)
! read tensor component
    read(50,*,err=20) a,b
    tmfix(i)=cmplx(a,b,8)
! read parity and Euler angles of spatial and spin rotation matrices
    read(50,'(A256)',err=20) str
    str=trim(str)//' 0.0 0.0 0.0'
    read(str,*,err=20) p,v1(:),v2(:)
    if (abs(p).ne.1) then
      write(*,*)
      write(*,'("Error(readinput): parity should be -1 or 1 in tmomfix &
       &block : ",I8)') p
      write(*,*)
      stop
    end if
! convert Euler angles from degrees to radians
    v1(:)=v1(:)*pi/180.d0
    v2(:)=v2(:)*pi/180.d0
! compute the spatial and spin 3x3 rotation matrices from the Euler angles
    call eulerrot(v1,rtmfix(:,:,1,i))
    call eulerrot(v2,rtmfix(:,:,2,i))
! multiply the spatial rotation matrix by the parity
    rtmfix(:,:,1,i)=dble(p)*rtmfix(:,:,1,i)
  end do
case('tauftm')
  read(50,*,err=20) tauftm
  if (tauftm.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): tauftm < 0 : ",G18.10)') tauftm
    write(*,*)
    stop
  end if
case('ftmstep')
  read(50,*,err=20) ftmstep
  if (ftmstep.lt.1) then
    write(*,*)
    write(*,'("Error(readinput): ftmstep < 1 : ",I8)') ftmstep
    write(*,*)
    stop
  end if
case('cmagz','forcecmag')
  read(50,*,err=20) cmagz
case('rotavec')
  read(50,*,err=20) axang(:)
case('tstime')
  read(50,*,err=20) tstime
  if (tstime.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): tstime <= 0 : ",G18.10)') tstime
    write(*,*)
    stop
  end if
case('dtimes')
  read(50,*,err=20) dtimes
  if (dtimes.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): dtimes <= 0 : ",G18.10)') dtimes
    write(*,*)
    stop
  end if
case('pulse')
  read(50,*,err=20) npulse
  if (npulse.lt.1) then
    write(*,*)
    write(*,'("Error(readinput): npulse < 1 : ",I8)') npulse
    write(*,*)
    stop
  end if
  if (allocated(pulse)) deallocate(pulse)
  allocate(pulse(8,npulse))
  do i=1,npulse
    read(50,*,err=20) pulse(:,i)
  end do
case('ramp')
  read(50,*,err=20) nramp
  if (nramp.lt.1) then
    write(*,*)
    write(*,'("Error(readinput): nramp < 1 : ",I8)') nramp
    write(*,*)
    stop
  end if
  if (allocated(ramp)) deallocate(ramp)
  allocate(ramp(6,nramp))
  do i=1,nramp
    read(50,*,err=20) ramp(:,i)
  end do
case('ncgga')
  read(50,*,err=20)
  write(*,*)
  write(*,'("Info(readinput): variable ''ncgga'' is no longer used")')
case('dncgga')
  read(50,*,err=20) dncgga
  if (dncgga.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): dncgga < 0 : ",G18.10)') dncgga
    write(*,*)
    stop
  end if
case('ntswrite')
  read(50,*,err=20) ntswrite
case('nxoapwlo','nxapwlo')
  read(50,*,err=20) nxoapwlo
  if (nxoapwlo.lt.0) then
    write(*,*)
    write(*,'("Error(readinput): nxoapwlo < 0 : ",I8)') nxoapwlo
    write(*,*)
    stop
  end if
case('nxlo')
  read(50,*,err=20) nxlo
  if (nxlo.lt.0) then
    write(*,*)
    write(*,'("Error(readinput): nxlo < 0 : ",I8)') nxlo
    write(*,*)
    stop
  end if
case('tdrho1d')
  read(50,*,err=20) tdrho1d
case('tdrho2d')
  read(50,*,err=20) tdrho2d
case('tdrho3d')
  read(50,*,err=20) tdrho3d
case('tdmag2d')
  read(50,*,err=20) tdmag2d
case('tdmag3d')
  read(50,*,err=20) tdmag3d
case('tdcd2d')
  read(50,*,err=20) tdcd2d
case('tdcd3d')
  read(50,*,err=20) tdcd3d
case('tddos')
  read(50,*,err=20) tddos
case('epseph')
  read(50,*,err=20)
  write(*,*)
  write(*,'("Info(readinput): variable ''epseph'' is no longer used")')
case('rndevt0')
  read(50,*,err=20) rndevt0
case('ssxc','rstsf')
  read(50,*,err=20) ssxc
case('tempk')
  read(50,*,err=20) tempk
  if (tempk.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): tempk <= 0 : ",G18.10)') tempk
    write(*,*)
    stop
  end if
! set Fermi-Dirac smearing
  stype=3
! set the smearing width
  swidth=kboltz*tempk
case('avecu')
  read(50,*,err=20) avecu(:,1)
  read(50,*,err=20) avecu(:,2)
  read(50,*,err=20) avecu(:,3)
case('kpamax')
  read(50,*,err=20) kpamax
  if (kpamax.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): kpamax <= 0 : ",G18.10)') kpamax
    write(*,*)
    stop
  end if
case('evtype')
  read(50,*,err=20) evtype
case('wmaxgw')
  read(50,*,err=20) wmaxgw
case('gwdiag')
  read(50,*,err=20) gwdiag
case('actype')
  read(50,*,err=20) actype
case('npole')
  read(50,*,err=20) npole
  if (npole.le.0) then
    write(*,*)
    write(*,'("Error(readinput): npole <= 0 : ",I8)') npole
    write(*,*)
    stop
  end if
case('nspade')
  read(50,*,err=20) nspade
  if (nspade.le.0) then
    write(*,*)
    write(*,'("Error(readinput): nspade <= 0 : ",I8)') nspade
    write(*,*)
    stop
  end if
case('tfav0')
  read(50,*,err=20) tfav0
case('nqssjx')
  read(50,*,err=20) nqssjx
  if (nqssjx.le.0) then
    write(*,*)
    write(*,'("Error(readinput): nqssjx <= 0 : ",I8)') nqssjx
    write(*,*)
    stop
  end if
case('thssjx')
  read(50,*,err=20) thssjx(:)
! convert from radians to degrees
  thssjx(:)=thssjx(:)*pi/180.d0
case('rmtall')
  read(50,*,err=20) rmtall
case('maxthd','omp_num_threads','OMP_NUM_THREADS')
  read(50,*,err=20) maxthd
case('maxlvl','omp_max_active_levels','OMP_MAX_ACTIVE_LEVELS')
  read(50,*,err=20) maxlvl
case('emaxulr')
  read(50,*,err=20) emaxulr
  if (emaxulr.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): emaxulr <= 0 : ",G18.10)') emaxulr
    write(*,*)
    stop
  end if
case('stable')
  read(50,*,err=20) lv
  if (lv) then
    nxoapwlo=1
    mixtype=3
    broydpm(1)=min(broydpm(1),0.01d0)
    broydpm(2)=min(broydpm(2),0.04d0)
    msmooth=max(msmooth,4)
  end if
case('metagga')
  read(50,*,err=20) lv
  if (lv) then
    nempty0=max(nempty0,10.d0)
    lradstp=1
    nrmtscf=max(nrmtscf,2.d0)
    msmooth=max(msmooth,4)
    if (mp_mpi) then
      write(*,*)
      write(*,'("Info(readinput): parameters set by metagga option")')
      write(*,'(" nempty0 : ",G18.10)') nempty0
      write(*,'(" lradstp : ",I4)') lradstp
      write(*,'(" nrmtscf : ",G18.10)') nrmtscf
      write(*,'(" msmooth : ",I4)') msmooth
    end if
  end if
case('taudft')
  read(50,*,err=20) taudft
case('t0tdlr')
  read(50,*,err=20) t0tdlr
  if (t0tdlr.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): t0tdlr <= 0 : ",G18.10)') t0tdlr
    write(*,*)
    stop
  end if
case('wannier')
  read(50,'(A256)',err=20) str ! wann_seedname
  if (trim(str).ne.'') then
    wann_seedname = adjustl(trim(str))
  end if
  read(50,*,err=20) wann_ngridk(:) ! Fine k-point grid
  if ((wann_ngridk(1).le.0).or.(wann_ngridk(2).le.0).or.(wann_ngridk(3).le.0)) then
    write(*,*)
    write(*,'("Error(readinput): invalid ngridk for Wannier : ",3I8)') wann_ngridk
    write(*,*)
    stop
  end if
  read(50,*,err=20) wann_nwf   ! wann_nwf
  read(50,'(A256)',err=20) str ! wann_bands
  call getw90bands(str)
  read(50,'(A256)',err=20) str ! wann_projections
  if (trim(str).eq.'') then
    write(*,*)
    write(*,'("Error(readinput): no projections given in wann_projections string")')
    write(*,*)
    stop
  end if
  wann_projlines = 0
  do
    wann_projlines = wann_projlines + 1
    wann_projstr(wann_projlines) = trim(str)
    read(50,'(A256)',err=20) str
    if (trim(str).eq.'') then
      goto 10
    else
      if(index(trim(str),'Ang').eq.0.and.index(trim(str),'Bohr').eq.0.and.&
        scan(trim(str),':').eq.0.and.index(trim(str),'random').eq.0) then
        write(*,*)
        write(*,'("Error(readinput): invalid projection : ",A)') trim(str)
        write(*,*)
        stop
      end if
    end if
  end do
case('wannierExtra')
  do
    read(50,'(A256)',err=20) str
    if (trim(str).eq.'') then
      goto 10
    else
      wann_inputlines = wann_inputlines + 1
      wann_input(wann_inputlines) = trim(str)
    end if
  end do
case('')
  goto 10
case default
  write(*,*)
  write(*,'("Error(readinput): invalid block name : ",A)') trim(block)
  write(*,*)
  stop
end select
goto 10
20 continue
write(*,*)
write(*,'("Error(readinput): error reading from elk.in")')
write(*,'("Problem occurred in ''",A,"'' block")') trim(block)
write(*,'("Check input convention in manual")')
write(*,*)
stop
30 continue
close(50)
! scale the speed of light
solsc=sol*solscf
! scale and rotate the lattice vectors (not referenced again in code)
avec(:,1)=sc1*avec(:,1)
avec(:,2)=sc2*avec(:,2)
avec(:,3)=sc3*avec(:,3)
avec(:,:)=sc*avec(:,:)
t1=axang(4)
if (t1.ne.0.d0) then
  t1=t1*pi/180.d0
  call axangrot(axang(:),t1,rot)
  do i=1,3
    v1(:)=avec(:,i)
    call r3mv(rot,v1,avec(:,i))
  end do
end if
! randomise lattice vectors if required
if (rndavec.gt.0.d0) then
  do i=1,3
    do j=1,3
      t1=rndavec*(randomu()-0.5d0)
      avec(i,j)=avec(i,j)+t1
    end do
  end do
end if
! case of isolated molecule
if (molecule) then
! convert atomic positions from Cartesian to lattice coordinates
  call r3minv(avec,ainv)
  do is=1,nspecies
    do ia=1,natoms(is)
      call r3mv(ainv,atposl(:,ia,is),v1)
      atposl(:,ia,is)=v1(:)
    end do
  end do
end if
! randomise atomic positions if required
if (rndatposc.gt.0.d0) then
  call r3minv(avec,ainv)
  do is=1,nspecies
    do ia=1,natoms(is)
      call r3mv(avec,atposl(:,ia,is),v1)
      do i=1,3
        t1=rndatposc*(randomu()-0.5d0)
        v1(i)=v1(i)+t1
      end do
      call r3mv(ainv,v1,atposl(:,ia,is))
    end do
  end do
end if
! randomise the muffin-tin magnetic fields if required
if (rndbfcmt.gt.0.d0) then
  do is=1,nspecies
    do ia=1,natoms(is)
      do i=1,3
        t1=rndbfcmt*(randomu()-0.5d0)
        bfcmt0(i,ia,is)=bfcmt0(i,ia,is)+t1
      end do
    end do
  end do
end if
! set fxctype to fxctype if required
if (fxctype(1).eq.-1) then
  fxctype(:)=xctype(:)
end if
! find primitive cell if required
if (primcell) call findprimcell
! read in atomic species data
call readspecies
return
end subroutine
!EOC
