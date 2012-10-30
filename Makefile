#
#  Makefile by Scott Prahl Jan 2012
#
VERSION = 2-3-4

CFLAGS = -Wall -ansi -dynamic -fno-common -g

LIB_EXT = .dylib		#MacOS X shared lib
#LIB_EXT = .so			#linux shared lib
LIB_EXT = .a			#static lib

#Base directory - adapt as needed
PREFIX=/usr/local
BIN_INSTALL=$(PREFIX)/bin
LIB_INSTALL=$(PREFIX)/lib
INC_INSTALL=$(PREFIX)/include
MMA_INSTALL=/Users/prahl/Library/Mathematica/Applications/Optics

#SHARED_LIB_OPT = -shared       #for redhat
SHARED_LIB_OPT = -bundle -flat_namespace -undefined suppress
DYNAMIC_LIB_OPT = -dynamiclib -install_name /usr/local/lib/libmie.dylib \
                  -compatibility_version 2.0 -current_version 2.0.0

MAIN = Makefile README

DOCS = 	doc/INSTALL doc/CHANGEFILE doc/mie_src.pdf

MMA = mma/Mietm.c mma/Mie.m mma/Mie.tm mma/Makefile

WSRC = 	src/mie_array.w   src/mie.w            src/test_latex.w  src/test_mie.w  \
		src/mie_root.w    src/mie_legendre.w   src/mie_doc.w     src/mie_main.w  \
		src/mie_lobatto.w src/mie_complex.w    src/mie_cylinder.w src/mie_cylinder_main.w

CSRC  = src/mie_array.c   src/mie_lobatto.c    src/mie_complex.c  \
        src/mie_root.c    src/mie_legendre.c   src/mie.c         src/mie_main.c  \
        src/test_mie.c    src/test_mie_array.c src/test_latex.c  src/test_mie_lobatto.c \
        src/mie_cylinder.c src/mie_cylinder_main.c

HSRC  = src/mie_array.h   src/libmie.h         src/mie.h         src/mie_complex.h \
        src/mie_root.h    src/mie_legendre.h   src/mie_lobatto.h src/mie_cylinder.h

OSRC  = src/system.bux    src/mie_doc.bux      src/cobweb.pl     src/version.pl \
        src/mygetopt.c    src/mygetopt.h       src/version.h     src/version.c \
        src/Makefile

INDENT_OPT = -bad -bap -c33 -ci1 -d1 -di20 -i4 -l78 -nbc -nfc1 -psl

all : mie

test : 
	cd src ; make test
	src/test_mie
	src/test_latex
	
doc docs doc/mie_src.pdf: $(WSRC)
	cd src ; make doc
	
lib : libmie.h libmie$(LIB_EXT)

mma: mma/Mie.exe

mie: $(WSRC) $(OSRC)
	cd src; make
	cp src/mie .
	cp src/miecyl .

mma/Mie.exe:
	cd mma ; make

install: mie
	mkdir -p $(BIN_INSTALL)
	cp mie  $(BIN_INSTALL)

libmie.h libmie$(LIB_EXT) : 
	cd src ; make libmie$(LIB_EXT)
	cp src/libmie.h .
	cp src/libmie$(LIB_EXT) .

install-lib: libmie$(LIB_EXT) libmie.h
	mkdir -p $(LIB_INSTALL)
	mkdir -p $(INC_INSTALL)
	cp src/libmie.h $(INC_INSTALL)
	cp src/libmie$(LIB_EXT) $(LIB_INSTALL)
	
install-mma: mma/Mie.m mma/Mie.exe
	mkdir -p $(MMA_INSTALL)
	mkdir -p $(MMA_INSTALL)/External
	cp mma/Mie.m $(MMA_INSTALL)
	cp mma/Mie.exe $(MMA_INSTALL)/External
	
install-all: mie libmie$(LIB_EXT) libmie.h mma/Mie.m mma/Mie.exe
	make install
	make install-lib
	make install-mma

old/mie_src.pdf mie_doc_old : $(WSRC)
	cd src ; ctwill -bhp mie_doc.w
	tex -output-directory src mie_doc.tex
	cd src ; refsort < mie_doc.ref > mie_doc.sref
	pdftex -output-directory src mie_doc.tex
	mv src/mie_doc.pdf src/mie_body.pdf
	cd src ; cweave -bhp mie_doc
	pdftex -output-directory src mie_doc.tex
	texexec --pdfselect --paper=letter --selection=51:57 --result src/mie_index.pdf src/mie_doc.pdf
	texexec --pdfselect --paper=letter --selection=58 --result src/mie_toc.pdf   src/mie_doc.pdf
	gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=doc/mie_src.pdf src/mie_toc.pdf src/mie_body.pdf src/mie_index.pdf &> /dev/null 
	rm -f src/mie_body.pdf src/mie_toc.pdf src/mie_index.pdf src/mie_doc.pdf
	rm -f src/mie_index.log src/mie_to.log
	rm -f texexec-mpgraph.mp texexec.tex texexec.tmp texexec.tui mpgraph.mp

tidy:
	cd src ; make tidy

dist: 
	make
	make test
	make doc
	make lib
	make mma
	mkdir -p       mie-$(VERSION)
	mkdir -p       mie-$(VERSION)/doc
	mkdir -p       mie-$(VERSION)/src
	mkdir -p       mie-$(VERSION)/mma
	ln $(MAIN)     mie-$(VERSION)
	ln $(DOCS)     mie-$(VERSION)/doc
	ln $(WSRC)     mie-$(VERSION)/src
	ln $(HSRC)     mie-$(VERSION)/src
	ln $(CSRC)     mie-$(VERSION)/src
	ln $(OSRC)     mie-$(VERSION)/src
	ln $(MMA)      mie-$(VERSION)/mma
#	tar cvf - mie-$(VERSION) | gzip  > mie-$(VERSION).tar.gz
	zip -r mie-$(VERSION) mie-$(VERSION)
	rm -rf mie-$(VERSION)
	
zip: 
	make test
	make tidy
	make doc
	make src
	mkdir -p	       xp-mie-$(VERSION)
	mkdir -p	       xp-mie-$(VERSION)/doc
	mkdir -p	       xp-mie-$(VERSION)/src
	ln $(MAIN)         xp-mie-$(VERSION)
	ln $(DOCS)         xp-mie-$(VERSION)/doc
	ln $(WSRC)         xp-mie-$(VERSION)/src
	ln $(HSRC)         xp-mie-$(VERSION)/src
	ln $(CSRC)         xp-mie-$(VERSION)/src
	ln $(OSRC)         xp-mie-$(VERSION)/src
	ln mie.exe	       xp-mie-$(VERSION)
	`perl -pi.bak -e 's/\n/\015\012/' xp-mie-$(VERSION)/src/*.c`
	`perl -pi.bak -e 's/\n/\015\012/' xp-mie-$(VERSION)/src/*.h`
	rm xp-mie-$(VERSION)/src/*.bak 
	zip -r xp-mie-$(VERSION) xp-mie-$(VERSION)
	rm -rf xp-mie-$(VERSION)

clean:
	cd mma  ; make clean
	cd src  ; make clean
	rm -f mie libmie$(LIB_EXT) libmie.h
	rm -f miecyl
	rm -f mpgraph.mp mprun.mp texexec.tmp texexec.tui
	rm -f test_mie test_mie_lobatto test_latex
	rm -f mie_doc.aux     mie_doc.idx     mie_doc.scn     mie_doc.toc
	rm -f mie_doc.log     mie_doc.sref 
	rm -f mie_doc.dvi     mie_doc.ref     mie_doc.tex     texexec.tex
	
realclean:
	make clean
	cd mma  ; make realclean
	cd src  ; make realclean
	rm -f *.tex *.dvi *.log *.aux *.sref *.o
	rm -f *.idx *.scn *.ref *.toc 
	rm -f *.pdf mie test_mie
	rm -f *.dylib
	rm -f $(CSRC) $(HSRC)
	rm -f mie_doc.pdf mie_body.pdf mie_toc.pdf mie_index.pdf
	rm -f doc/mie_src.pdf
	
help::
	@echo;\
	echo "Targets available for this Makefile:";\
	echo "clean       -- remove most generated objects";\
	echo "docs        -- generate TEX out of all files";\
	echo "dist        -- gzipped tarball of all files";\
	echo "install     -- install ad and iad programs";\
	echo "install-lib -- install interface and library programs";\
	echo "install-mma -- install Mathematica files";\
	echo "install-all -- install all the above";\
	echo "lib         -- library binary and interface files";\
	echo "mie         -- compile Mie scattering program";\
	echo "mma         -- adding-doubling files for Mathematica interface";\
	echo "realclean   -- remove all generated objects";\
	echo "test        -- run iad program on a bunch of test files";\
	echo "tidy        -- generate .c and .h files"

src/version.c :
	cd src ; ./version.pl
	
.PHONY: clean tidy dist doc docs lib mma test install install-all mie_doc

.SECONDARY : $(CSRC) $(HSRC)
