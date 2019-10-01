#
#  Makefile by Scott Prahl in Nov 2010
#
#  requires libmie.a and  libmie.h
#
VERSION      = 8.0
SYS          = MacOSX

INSTALL_DIR  = /Users/prahl/Library/Mathematica/Applications/Optics
MLINKDIR     = /Applications/Mathematica.app/Contents/SystemFiles/Links/MathLink/DeveloperKit
CADDSDIR     = ${MLINKDIR}/MacOSX-x86-64/CompilerAdditions
INC_DIR      = ${CADDSDIR}
LIB_DIR      = ${CADDSDIR}
EXTRA_LIBS   = -stdlib=libc++ -lc++ -framework Foundation -lmie
MPREP        = ${CADDSDIR}/mprep

all : Mie

.c.o :
	${CC} -c -I${INC_DIR}  -I/usr/local/include $<

Mietm.c : Mie.tm
	${MPREP} $? -o $@

Mie : Mietm.o
	${CXX} -I${INC_DIR} -I/usr/local/include Mietm.o -L${LIB_DIR} -lMLi4 ${EXTRA_LIBS} -o $@

install:
	mkdir -p $(INSTALL_DIR)
	mkdir -p $(INSTALL_DIR)/External
	cp Mie.m $(INSTALL_DIR)
	cp Mie   $(INSTALL_DIR)/External

realclean:
clean:
	rm -f Mietm.o Mietm.c Mie

.PHONY : clean realclean
