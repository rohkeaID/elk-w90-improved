
include ../make.inc

# modules
SRC_modules = modmain.f90 $(SRC_MPI) modmpi.f90 $(SRC_libxc) modxcifc.f90 \
 modfxcifc.f90 moddftu.f90 modrdm.f90 modphonon.f90 modscdft.f90 modtest.f90 \
 modrandom.f90 modstore.f90 modpw.f90 modvars.f90 modtddft.f90 modgw.f90 \
 modulr.f90 modjx.f90 modw90.f90 modw90overlap.f90 $(SRC_OMP) modomp.f90

# Elk program
SRC_elk = elk.f90

# main subroutines and functions
SRC_main0 = \
 zftwfir.f90 zfinp.f90 genwfsv.f90 genwfsvp.f90 match.f90 gengclg.f90 \
 gengclgq.f90 allatoms.f90 genzrho.f90 potcoul.f90 rfint.f90 gensfacgp.f90 \
 genjlgqr.f90 gengvec.f90 wavefmt.f90 rbsht.f90 rfsht.f90 rhomagsh.f90 \
 zfmtconj.f90 rhomag.f90 genbs.f90 genws.f90 zgerci.f90 spline.f90 sphcrd.f90 \
 gencfun.f90 potnucl.f90 genffacgp.f90 potxc.f90 zbsht.f90 zfsht.f90 \
 zfmtinp.f90 zftzf.f90 rhonorm.f90 rfinp.f90 writelat.f90 findkpt.f90 \
 potks.f90 fsmfield.f90 gradzf.f90 gradrf.f90 checkfsm.f90 eveqnsv.f90 \
 brzint.f90 rfmtsm.f90 rfirsm.f90

SRC_main1 = \
 gentpmae.f90 genstrain.f90 writeiad.f90 symveca.f90 symvec.f90 nuclei.f90 \
 genstress.f90 writestress.f90 latvstep.f90 putoccsv.f90 getoccsv.f90 \
 linengy.f90 elnes.f90 potefield.f90 writeinfo.f90 readinput.f90 hermite.f90 \
 writelinen.f90 writesf.f90 wigner3j.f90 wigner3jf.f90 wigner6j.f90 \
 clebgor.f90 gaunt.f90 gauntyry.f90 roteuler.f90 \
 eulerrot.f90 reciplat.f90 fermisurfbxsf.f90 writeeval.f90 readspecies.f90 \
 gencore.f90 genhvec.f90 writeforces.f90 pade.f90 pades.f90 \
 genjlgprmt.f90 symmetry.f90 findsym.f90 findsymcrys.f90 findsymsite.f90 \
 writekpts.f90 olpistl.f90 olpaa.f90 olpalo.f90 olplolo.f90 hmlaa.f90 \
 hmlalo.f90 writefsm.f90 plotpt1d.f90 testcheck.f90 \
 atom.f90 geomopt.f90 atpstep.f90 bandstr.f90 findsymlat.f90 \
 putevalfv.f90 getevalfv.f90 mixerifc.f90 axangsu2.f90 wfplot.f90 \
 gengclq.f90 writegclq.f90 writeefg.f90 gensocfr.f90 nesting.f90 \
 putevalsv.f90 getevalsv.f90 dos.f90 writedos.f90 putevecfv.f90 getevecfv.f90 \
 putevecsv.f90 getevecsv.f90 fermisurf.f90 writesym.f90 genidxlo.f90 \
 addlorbcnd.f90 sdelta.f90 stheta.f90 sdelta_mp.f90 stheta_mp.f90 \
 sdelta_fd.f90 stheta_fd.f90 sdelta_sq.f90 stheta_sq.f90 sdelta_lr.f90 \
 stheta_lr.f90 vecfbz.f90 mixlinear.f90 mixadapt.f90 mixbroyden.f90 rlmrot.f90 \
 rotrflm.f90 radnucl.f90 findngkmax.f90 writegeom.f90 effmass.f90 \
 mossbauer.f90 findband.f90 geomplot.f90 genlmirep.f90 gndstate.f90 \
 findqpt.f90 ssfext.f90 sstask.f90 spiralsc.f90 genscss.f90 \
 findscq.f90 writeqpts.f90 writelsj.f90 getcfgq.f90 proj2d.f90 dielectric.f90 \
 nonlinopt.f90 vecplot.f90 gengkvec.f90 wfcrplot.f90 energykncr.f90 \
 rhoplot.f90 eveqn.f90 eveqnit.f90 \
 eveqnss.f90 genevfsv.f90 rhoinit.f90 maginit.f90 grad2rfmt.f90 gendmat.f90 \
 symdmat.f90 gensdmat.f90 ggamt_1.f90 ggair_1.f90 ggamt_sp_1.f90 \
 ggair_sp_1.f90 ggamt_2a.f90 ggair_2a.f90 ggamt_2b.f90 ggair_2b.f90 \
 ggamt_sp_2a.f90 ggair_sp_2a.f90 ggamt_sp_2b.f90 ggair_sp_2b.f90 \
 genspecies.f90 writeemd.f90 rhomagk.f90 wavefcr.f90 exxengy.f90 \
 exxengyk.f90 eveqnhf.f90 gentau.f90 gentauk.f90 gentaucr.f90 \
 xc_c_tb09.f90 elfplot.f90 genexpmt.f90 gradzfmt.f90 \
 gradrfmt.f90 projsbf.f90 genpmatk.f90 putpmat.f90 \
 writepmat.f90 getpmat.f90 gengqrf.f90 moke.f90 sfacinit.f90 sfacrho.f90 \
 sfacmag.f90 gradwf2.f90 gradwfcr2.f90 force.f90 forcek.f90 symrf.f90 \
 symrfir.f90 symrvfir.f90 symrvf.f90 oepmain.f90 oepresk.f90 \
 oepvcl.f90 oepvclk.f90 genvchi0.f90 genspchi0.f90 epsinv.f90 dbxcplot.f90 \
 genvmat.f90 genvmatk.f90 genvbmatk.f90 vclcore.f90 genzrm.f90 genzfrm.f90 \
 zftrf.f90 writehmlbse.f90 hmlxbse.f90 hmlxbsek.f90 hmldbse.f90 \
 hmldbsek.f90 dielectric_bse.f90 genidxbse.f90 writeevbse.f90 hartfock.f90 \
 hflocal.f90 putkmat.f90 getkmat.f90 genwfpw.f90 writewfpw.f90 getwfpw.f90 \
 genfspecies.f90 writespecies.f90 massnucl.f90 genexpmat.f90 \
 writeexpmat.f90 rvfcross.f90 rfplot.f90 rotdmat.f90 hrmdmat.f90 rzfmtinp.f90 \
 rzfinp.f90 genpmat.f90 reademd.f90 emdplot.f90 rfhkintp.f90 \
 emdplot3d.f90 emdplot2d.f90 emdplot1d.f90 plotpt3d.f90 plotpt2d.f90 \
 genvmmtsv.f90 writeevsp.f90 sctovec.f90 rndevsv.f90 gengvc.f90 \
 straingkq.f90 eveqnz.f90 ztpmm.f90 writeepsinv.f90 \
 torque.f90 ssfsmjx.f90 gencfrc.f90 curden.f90 curdenk.f90 closefiles.f90 \
 curdenplot.f90 wxcplot.f90 nfftifc.f90 zfcmtinp.f90 sbessel.f90 \
 hmlrad.f90 olprad.f90 checkmt.f90 genylmg.f90 genzvclmt.f90 \
 ztorflm.f90 fderiv.f90 fintgt.f90 factnm.f90 factr.f90 mtdmin.f90 \
 fsmooth.f90 zpotclmt.f90 rfmtpack.f90 zfmtpack.f90 zfpack.f90 rfmtinp.f90 \
 charge.f90 chargemt.f90 rfmtlm.f90 zfmtlm.f90 genlofr.f90 \
 rdiracint.f90 rdirac.f90 genrmesh.f90 energynn.f90 rhocore.f90 \
 rtozfmt.f90 ztorfmt.f90 init0.f90 init1.f90 init2.f90 init3.f90 init4.f90 \
 rfinpc.f90 gcd.f90 sortidx.f90 erf.f90 zpotcoul.f90 \
 genylm.f90 genrlm.f90 rtozflm.f90 moment.f90 occupy.f90 \
 rfmtftoc.f90 rschrodint.f90 timesec.f90 lopzflm.f90 hmlistl.f90 sphcover.f90 \
 delevec.f90 hmllolo.f90 genapwfr.f90 mixpack.f90 rfpack.f90 polynom.f90 \
 genkmat.f90 genshtmat.f90 eveqnfv.f90 eveqnfvz.f90 \
 writeengy.f90 gridsize.f90 genvsig.f90 sbesseldm.f90 \
 writefermi.f90 readfermi.f90 zfmtctof.f90 rotrfmt.f90 zfmtftoc.f90 \
 rfinterp.f90 ylmrot.f90 ylmroty.f90 eveqnfvr.f90 findprimcell.f90 genppts.f90 \
 energy.f90 rfmtctof.f90 writechg.f90 potplot.f90 plot1d.f90 plot2d.f90 \
 plot3d.f90 findswidth.f90 mae.f90 writestate.f90 readstate.f90 rotaxang.f90 \
 axangrot.f90 i3minv.f90 i3mdet.f90 i3mtv.f90 r3mv.f90 r3mtv.f90 r3cross.f90 \
 r3minv.f90 r3mdet.f90 r3frac.f90 r3mm.f90 r3mtm.f90 r3mmt.f90 zflmconj.f90 \
 rotzflm.f90 z2mm.f90 z2mctm.f90 z2mmct.f90 symmat.f90 r3vo.f90

SRC_main = $(SRC_main0) $(SRC_main1)

SRC_ULR = \
 initulr.f90 genkpakq.f90 eveqnulr.f90 genevulr.f90 gndstulr.f90 rhomaguk.f90 \
 genidxulr.f90 genzvmatk.f90 genzvbmatk.f90 potclulr.f90 rhomagq.f90

SRC_phonon0 = \

SRC_phonon1 = \
 phdelete.f90 gendcfun.f90 dmatch.f90 gendvsig.f90 gradrhomt.f90 dpotxc.f90 \
 dolpalo.f90 dhmlalo.f90 dolpaa.f90 dhmllolo.f90 dwavefmt.f90 phonon.f90 \
 dpotks.f90 dhmlaa.f90 dhmlrad.f90 gendwfsv.f90 drhomagk.f90 drhomagsh.f90 \
 hmlaaq.f90 hmlaloq.f90 olpaaq.f90 olpaloq.f90 gengvsmt.f90 dpotcoul.f90 \
 dforcek.f90 gengkqvec.f90 gengqvec.f90 genmcph.f90 deveqnfv.f90 dhmlistl.f90 \
 dolpistl.f90 phmixpack.f90 putdevecfv.f90 getdevecfv.f90 doccupy.f90 \
 putdevecsv.f90 getdevecsv.f90 phononsc.f90 dyntask.f90 genscph.f90 phdisp.f90 \
 readdyn.f90 dynqtor.f90 dynrtoq.f90 dynsymapp.f90 dynsym.f90 dynev.f90 \
 dynevs.f90 phdos.f90 sumrule.f90 writephn.f90 phfext.f90 writegamma.f90 \
 readgamma.f90 phlwidth.f90 alpha2f.f90 writelambda.f90 phscdelete.f90 \
 ephcouple.f90 genephmat.f90 eliashberg.f90 phscdvs.f90 writedvs.f90 \
 readdvs.f90 readalpha2f.f90 mcmillan.f90 dforce.f90 gendsocfr.f90

SRC_phonon = $(SRC_phonon0) $(SRC_phonon1)

SRC_SCDFT = \
 hmlbdg.f90 eveqnbdg.f90 scdft.f90 genachi.f90 bdginv.f90 achiinit.f90 \
 genidxbdg.f90

SRC_DFTU = \
 gendmatmt.f90 genveedu.f90 genvmatmt.f90 writedftu.f90 readalphadu.f90 \
 fyukawa.f90 fyukawa0.f90 genfdufr.f90 energyfdu.f90 writeefdu.f90 \
 zbessela.f90 zbesselb.f90 readdmatmt.f90 genfdu.f90 findlambdadu.f90 \
 readfdu.f90 writetmdu.f90 writetm2du.f90 writetm3du.f90 dmtotm2.f90 \
 dmtotm3.f90 tm2todm.f90 tm3todm.f90 pottm2.f90 pottm3.f90 tm2pol.f90 \
 tm3pol.f90 gendmftm.f90 ftmfield.f90 writeftm.f90

SRC_RDMFT = \
 rdmft.f90 rdmdkdc.f90 rdmminc.f90 rdmvaryc.f90 rdmdedc.f90 rdmenergy.f90 \
 rdmengyxc.f90 rdmwritededn.f90 rdmwriteengy.f90 rdmminn.f90 rdmvaryn.f90 \
 rdmdedn.f90 rdmdexcdn.f90 genvclijji.f90 genvclijjk.f90 getvclijji.f90 \
 getvclijjk.f90 writevclijji.f90 writevclijjk.f90 rdmdexcdc.f90 rdmdtsdn.f90 \
 rdmentropy.f90 rdmeval.f90

SRC_GW0 = \
 gwlocal.f90 genwgw.f90 putgwsefm.f90 getgwsefm.f90 gwsefm.f90 gwsefmk.f90

SRC_GW1 = \
 dysonr.f90 gwspecf.f90 writegwsf.f90 acgwse.f90 minf_nm.f90 fmin_nm.f90 \
 polefit.f90 acpole.f90 zfpole.f90 gwbandstr.f90 gwefermi.f90

SRC_GW = $(SRC_GW0) $(SRC_GW1)

SRC_TDDFT0 = \
 tddft.f90 potkst.f90 timestep.f90 genhmlt.f90 tddftlr.f90 genvfxc.f90 \
 genfxcr.f90 tdtemp.f90 writetimes.f90 readtimes.f90

SRC_TDDFT1 = \
 genspfxcr.f90 genvfxcg.f90 gentimes.f90 \
 tddftsplr.f90 genspfxcg.f90 genafieldt.f90 readafieldt.f90 writetddft.f90 \
 writeafpdt.f90 writetddos.f90 dielectric_tdrt.f90 readcurrent.f90

SRC_TDDFT = $(SRC_TDDFT0) $(SRC_TDDFT1)

SRC_XC = \
  xc_pzca.f90 xc_pwca.f90 xc_pbe.f90 xc_am05.f90 xc_xalpha.f90 xc_wc06.f90 \
  x_wc06.f90 x_pbe.f90 c_pbe.f90 c_pbe_gcor.f90 xc_vbh.f90 fxc_pwca.f90

SRC_WANN = \
 getw90bands.f90 genw90twf.f90 \
 writew90win.f90 genw90amn.f90 \
 setupw90lib.f90 genw90input.f90 runw90.f90

SRC = $(SRC_modules) $(SRC_FFT) $(SRC_elk) $(SRC_main) $(SRC_phonon) \
 $(SRC_SCDFT) $(SRC_XC) $(SRC_DFTU) $(SRC_RDMFT) $(SRC_GW) $(SRC_TDDFT) \
 $(SRC_ULR) $(SRC_WANN)

OBJ = $(SRC:.f90=.o)
EXE = elk

#-------------------------------------------------------------------------------
# Suffix rules
#-------------------------------------------------------------------------------
.SUFFIXES: .o .f90
.f90.o:
	$(F90) $(F90_OPTS) -c $<

#-------------------------------------------------------------------------------
# Targets
#-------------------------------------------------------------------------------
elk:	$(OBJ)
	$(F90) $(F90_OPTS) -o $(EXE) $(OBJ) $(LIB_W90) $(LIB_LPK) $(LIB_FFT) $(LIB_libxc) \
         $(LIB_SYS)  

blas:
	cd BLAS; $(MAKE); cp blas.a ..

lapack:
	cd LAPACK; $(MAKE); cp lapack.a ..

fft:
	cd fftlib; $(MAKE); cp fftlib.a ..

all:	blas lapack fft elk

clean:
	rm -f *.o *.mod *~ ifc* *.gcno gmon.out *.aux *.dvi *.log \
         elk.pdf *.tex *.toc $(EXE)

cleanall:
	cd BLAS; $(MAKE) clean
	cd LAPACK; $(MAKE) clean
	cd fftlib; $(MAKE) clean
	$(MAKE) clean

doc:
	rm -f elk.aux elk.bbl elk.blg elk.log elk.pdf elk.tex elk.toc elk.lst
	ls $(SRC_modules) $(SRC_main) $(SRC_XC) $(SRC_phonon) $(SRC_DFTU) \
         $(SRC_RDMFT) $(SRC_TDDFT) > elk.lst
	./protex -F -s $(SRC_elk) $$(cat elk.lst) > elk.tex
	pdflatex elk;pdflatex elk;pdflatex elk

backup:
	tar -czf elk.tgz *.f90 fftlib eos spacegroup junk BLAS LAPACK \
         Makefile notes.txt ../docs ../examples ../tests protex ../COPYING \
         ../README elk_silhouette.pdf depend gentest genalltests \
         checklist.txt rmspaces grepelk vimelk

VERSION = $$(awk -F"/" '/data version/ {print $$2}' modmain.f90 | sed 's/ //g;s/,/./g')$
NAME = elk-$(VERSION)
RELEASE = ../release/$(NAME)

release:
	rm -rf $(RELEASE)
	mkdir $(RELEASE)
	cd BLAS; $(MAKE) clean
	cd LAPACK; $(MAKE) clean
	cd fftlib; $(MAKE) clean
	cd eos; $(MAKE) clean
	cd spacegroup; $(MAKE) clean
	cd ..; find . -type f -name *~ -exec rm -f {} \;
	cd ..; find . -type f -name gmon.out -exec rm -f {} \;
	cd ..; find . -type f -name RUNNING -exec rm -f {} \;
	mkdir $(RELEASE)/src
	cp -r BLAS $(RELEASE)/src
	cp -r LAPACK $(RELEASE)/src
	cp -r fftlib $(RELEASE)/src
	cp -r eos $(RELEASE)/src
	cp -r spacegroup $(RELEASE)/src
	cp $(SRC) Makefile protex elk_silhouette.pdf $(RELEASE)/src
	cp $(SRC) rmspaces grepelk vimelk $(RELEASE)/src
	cp mpi_stub.f90 omp_stub.f90 $(RELEASE)/src
	cp libxc_funcs.f90 libxc.f90 libxcifc.f90 libxcifc_stub.f90 $(RELEASE)/src
	cp zfftifc.f90 zfftifc_fftw.f90 zfftifc_mkl.f90 $(RELEASE)/src
	cp -r ../examples $(RELEASE)
	cp -r ../tests $(RELEASE)
	mkdir $(RELEASE)/species
	cp ../species/*.in $(RELEASE)/species
	cp -r ../utilities $(RELEASE)
	cp ../COPYING $(RELEASE)
	cp ../README $(RELEASE)
	cp ../release_notes.txt $(RELEASE)
	cp ../Makefile $(RELEASE)
	cp ../setup $(RELEASE)
	cp ../make.def $(RELEASE)
	$(MAKE) doc
	cp elk.pdf ../docs
	cd spacegroup;$(MAKE) doc;cp spacegroup.pdf ../../docs;$(MAKE) clean
	cp -r ../docs $(RELEASE)
	mkdir $(RELEASE)/release
	tar -C ../release -czf ../release/$(NAME).tgz $(NAME)
	cp ../docs/elk.pdf ../release
	cp ../docs/spacegroup.pdf ../release

lines:
	cat $(SRC) | wc -l

spaces:
	./rmspaces $(SRC)

vim:
	./vimelk

