
include ../../make.inc

#-------------------------------------------------------------------------------
# Suffix rules
#-------------------------------------------------------------------------------
.SUFFIXES: .o .f90
.f90.o:
	$(F90) $(F90_OPTS) -c $<

#-------------------------------------------------------------------------------
# Source files
#-------------------------------------------------------------------------------

SRC_modules = modmain.f90

SRC_main = main.f90

SRC_routines = readinput.f90 gencrystal.f90 sgsymb.f90 seitzgen.f90 \
 gengroup.f90 seitzeq.f90 seitzmul.f90 writegeom.f90 geomplot.f90 \
 findprimcell.f90 r3frac.f90 r3mv.f90 r3cross.f90 r3minv.f90 r3mm.f90

SRC = $(SRC_modules) $(SRC_main) $(SRC_routines)

OBJ = $(SRC:.f90=.o)
EXE = spacegroup

spacegroup:	$(OBJ)
	$(F90) $(F90_OPTS) -o $(EXE) $(OBJ)

clean:
	rm -f *.o *.mod *~ fort.* ifc* *.gcno *.OUT *.xsf *.ascii gmon.out \
         *.aux *.dvi *.log *.pdf *.tex *.toc $(EXE)

doc:
	../protex -s $(SRC_main) $(SRC_modules) $(SRC_routines) > spacegroup.tex
	pdflatex spacegroup;pdflatex spacegroup;pdflatex spacegroup

