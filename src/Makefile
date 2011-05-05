#
#  Makefile by Scott Prahl May 2006
#

# flags for debugging
#CFLAGS = -Wall -ansi -fno-common -pedantic -W -Wconversion -Wshadow -g
# typical flags
CFLAGS = -fno-common -O3

LOCAL_LIBRARIES = -lm

SHARED_LIB_OPT = -shared       #for redhat
DYNAMIC_LIB_OPT = -dynamiclib -install_name /usr/local/lib/libmie.dylib \
                  -compatibility_version 2.0 -current_version 2.0.0

WSRC = 	mie_array.w   mie.w            test_latex.w  test_mie.w  \
		mie_root.w    mie_legendre.w   mie_doc.w     mie_main.w  \
		mie_lobatto.w mie_complex.w    mie_cylinder.w mie_cylinder_main.w

CSRC  = mie_array.c   mie_lobatto.c    mie_complex.c  \
        mie_root.c    mie_legendre.c   mie.c         mie_main.c  \
        test_mie.c    test_mie_array.c test_latex.c  test_mie_lobatto.c \
        mie_cylinder.c mie_cylinder_main.c

HSRC  = mie_array.h   libmie.h         mie.h         mie_complex.h \
        mie_root.h    mie_legendre.h   mie_lobatto.h mie_cylinder.h   

OSRC  = system.bux    mie_doc.bux      cobweb.pl     version.pl \
        mygetopt.c    mygetopt.h       version.h     version.c

OBJ   = mie_array.o   mie.o            mie_complex.o mie_legendre.o  \
        mie_root.o    mie_lobatto.o    mygetopt.o    version.o

INDENT_OPT = -bad -bap -c33 -ci1 -d1 -di20 -i4 -l78 -nbc -nfc1 -psl

all : mie miecyl

test: test_mie test_mie_lobatto test_latex

mie: version $(OBJ)  mie_main.o	
	$(CC) -o mie $(OBJ) mie_main.o $(LOCAL_LIBRARIES)

miecyl : mie_cylinder_main.o mie_cylinder.o mie_complex.o mie_array.o
	$(CC) -o miecyl mie_cylinder_main.o mie_cylinder.o mie_complex.o mie_array.o
	
libmie.h: mie.w
	ctangle -bph mie
	./cobweb.pl libmie.h
	indent libmie.h 
	rm -f libmie.h.bak libmie.h.BAK
	
libmie.a: $(OBJ) libmie.h
	ar rs libmie.a $(OBJ)

libmie.so: $(OBJ) libmie.h
	cc $(SHARED_LIB_OPT) -o libmie.so $(OBJ)

libmie.dylib: $(OBJ) libmie.h
	cc $(DYNAMIC_LIB_OPT) -o libmie.dylib $(OBJ)

test_mie: version $(OBJ) test_mie.o test_mie.c
	$(CC) -o test_mie $(OBJ) test_mie.o $(LOCAL_LIBRARIES)

test_mie_lobatto: version $(OBJ) test_mie_lobatto.o
	$(CC) -o test_mie_lobatto $(OBJ) test_mie_lobatto.o $(LOCAL_LIBRARIES)

test_latex: version $(OBJ) test_latex.o
	$(CC) -o test_latex $(OBJ) test_latex.o $(LOCAL_LIBRARIES)

tidy:
	apply ctangle *.w
	./version.pl
	apply ./cobweb.pl $(CSRC) $(HSRC)
	apply 'indent $(INDENT_OPT)' mie_array.c mie_lobatto.c mie_complex.c mie_root.c
	apply 'indent $(INDENT_OPT)' mie_legendre.c mie.c mie_main.c test_mie.c
	apply 'indent $(INDENT_OPT)' test_mie_array.c test_latex.c test_mie_lobatto.c
	apply 'indent $(INDENT_OPT)' mie_array.h mie_lobatto.h mie_complex.h mie_root.h
	apply 'indent $(INDENT_OPT)' mie_legendre.h mie.h 
	rm -f *.bak *.BAK

clean:
	rm -f libmie$(LIB_EXT) libmie.h
	rm -f mpgraph.mp mprun.mp texexec.tmp texexec.tui
	rm -f $(OBJ) test_latex.o test_mie_lobatto.o mie_main.o test_mie.o
	rm -f test_latex test_mie_lobatto test_mie mie
	rm -f mie_cylinder.o mie_cylinder_main.o miecyl
	
realclean:
	make clean
	rm -f *.tex *.dvi *.log *.aux *.sref *.o
	rm -f *.idx *.scn *.ref *.toc 
	rm -f *.dylib
	rm -f $(CSRC) $(HSRC)
	
doc docs: $(WSRC)
	cweave mie_doc
	pdftex mie_doc
	mv mie_doc.pdf ../doc/mie_src.pdf

mie_array.h :  mie_array.w
	ctangle -bph mie_array.w
	
mie_complex.h :  mie_complex.w
	ctangle -bph mie_complex.w

mie.h mie.c:  mie.w
	ctangle -bph mie.w
	
mie.o : mie_complex.h mie.c mie.h 
	$(CC) $(CFLAGS)   -c -o mie.o mie.c
	
mie_cylinder.h : mie_cylinder.w
	ctangle -bph mie_cylinder.w

mie_cylinder.o: mie_cylinder.h mie_complex.h mie_array.h
	$(CC) $(CFLAGS) -c -o mie_cylinder.o mie_cylinder.c

mie_cylinder_main.o: mie_cylinder_main.c mie_cylinder.h mie_complex.h mie_array.h
	$(CC) $(CFLAGS) -c -o mie_cylinder_main.o mie_cylinder_main.c

test_mie_lobatto.c : mie_lobatto.w
	ctangle -bph mie_lobatto

test_mie.c : test_mie.w
	ctangle -bph test_mie
	
version :
	./version.pl
	
.PHONY : clean dist doc docs lib mma test install install-all test mie_doc

.SECONDARY : $(CSRC) $(HSRC)
