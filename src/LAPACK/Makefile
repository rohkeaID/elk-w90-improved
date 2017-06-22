
AR = ar

include ../../make.inc

#-------------------------------------------------------------------------------
# Suffix Rules
#-------------------------------------------------------------------------------
.SUFFIXES: .o .f
.f.o:
	$(F77) $(F77_OPTS) -c $<

#-------------------------------------------------------------------------------
# File dependencies
#-------------------------------------------------------------------------------
SRC = \
 dgesv.f     dgetf2.f  dgetrf.f	 dgetri.f  dgetrs.f  disnan.f  dlabad.f \
 dlacpy.f    dladiv.f  dlae2.f	 dlaebz.f  dlaev2.f  dlagtf.f  dlagts.f \
 dlaisnan.f  dlamch.f  dlanst.f	 dlansy.f  dlapy2.f  dlapy3.f  dlarfb.f \
 dlarf.f     dlarfg.f  dlarft.f	 dlarnv.f  dlartg.f  dlaruv.f  dlascl.f \
 dlaset.f    dlasr.f   dlasrt.f	 dlassq.f  dlaswp.f  dlasyf.f  dlatrd.f \
 dorg2l.f    dorg2r.f  dorgql.f	 dorgqr.f  dorgtr.f  dorm2l.f  dorm2r.f \
 dormql.f    dormqr.f  dormtr.f	 dpotf2.f  dpotrf.f  dstebz.f  dstein.f \
 dsteqr.f    dsterf.f  dsyconv.f  dsyev.f   dsyevx.f  dsygs2.f  dsygst.f \
 dsygvx.f    dsysv.f   dsytd2.f	 dsytf2.f  dsytrd.f  dsytrf.f  dsytrs2.f \
 dsytrs.f    dtrti2.f  dtrtri.f	 ieeeck.f  iladlc.f  iladlr.f  ilaenv.f \
 ilazlc.f    ilazlr.f  iparmq.f	 xerbla.f  zgebak.f  zgebal.f  zgeev.f \
 zgehd2.f    zgehrd.f  zgesv.f	 zgetf2.f  zgetrf.f  zgetri.f  zgetrs.f \
 zheev.f     zheevx.f  zhegs2.f	 zhegst.f  zhegv.f   zhegvx.f  zhetd2.f \
 zhetrd.f    zhpevx.f  zhpgst.f	 zhpgvx.f  zhptrd.f  zhseqr.f  zlacgv.f \
 zlacpy.f    zladiv.f  zlahqr.f	 zlahr2.f  zlange.f  zlanhe.f  zlanhp.f \
 zlaqr0.f    zlaqr1.f  zlaqr2.f	 zlaqr3.f  zlaqr4.f  zlaqr5.f  zlarfb.f \
 zlarf.f     zlarfg.f  zlarft.f	 zlartg.f  zlascl.f  zlaset.f  zlasr.f \
 zlassq.f    zlaswp.f  zlatrd.f	 zlatrs.f  zpotf2.f  zpotrf.f  zpptrf.f \
 zrot.f	    zstein.f  zsteqr.f	 ztrevc.f  ztrexc.f  ztrti2.f  ztrtri.f \
 zung2l.f    zung2r.f  zunghr.f	 zungql.f  zungqr.f  zungtr.f  zunm2l.f \
 zunm2r.f    zunmhr.f  zunmql.f	 zunmqr.f  zunmtr.f  zupgtr.f  zupmtr.f \
 dbdsqr.f    dgebd2.f  dgebrd.f  dgelq2.f  dgelqf.f  dgeqr2.f  dgeqrf.f \
 dgesvd.f    dlabrd.f  dlange.f  dlas2.f   dlasq1.f  dlasq2.f  dlasq3.f \
 dlasq4.f    dlasq5.f  dlasq6.f  dlasv2.f  dorgbr.f  dorgl2.f  dorglq.f \
 dormbr.f    dorml2.f  dormlq.f  lsame.f

OBJ = $(SRC:.f=.o)

lapack:	$(OBJ)
	$(AR) -rc lapack.a $(OBJ)

clean:
	rm -f *.o *.mod *~ *.a ifc* *.gcno gmon.out

ls:
	ls -x --width=80 *.f
