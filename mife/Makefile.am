AUTOMAKE_OPTIONS = foreign -Wall
SUBDIRS = jsmn

MY_SOURCES = mife.c mife_io.c mife_internals.c flint_raw_io.c parse.c \
             mbp_glue.c util.c mbp_types.c jsmn/jsmn.c cmdline.c

AM_CFLAGS = ${DEBUG_CFLAGS} -I$(top_srcdir) -I$(top_srcdir)/dgs -I$(top_srcdir)/flint \
			-I$(top_builddir)/mife/jsmn -fopenmp

MY_LIB_LDFLAGS=-release 0.0.$(RELEASE) -no-undefined

LDADD = \
	$(top_builddir)/dgs/libdgs.la \
	$(top_builddir)/oz/liboz.la \
	$(top_builddir)/dgsl/libdgsl.la \
	-L$(top_builddir)/mife/jsmn \
	$(top_builddir)/gghlite/libgghlite.la \
	-L$(top_builddir)/mmap \
	-laesrand -lflint -ljsmn -lgmp -lmpfr -lm -lssl -lcrypto -lclt13 -lmmap

pkgincludesubdir = $(includedir)/mife
pkgincludesub_HEADERS = mife.h mife_io.h mife_internals.h flint_raw_io.h \
	mife_defs.h parse.h mbp_glue.h util.h cmdline.h mbp_types.h

AM_LDFLAGS = -rpath $(abs_top_builddir)/flint

bin_PROGRAMS = keygen encrypt eval
keygen_SOURCES  =  keygen.c $(MY_SOURCES)
encrypt_SOURCES = encrypt.c $(MY_SOURCES)
eval_SOURCES    =    eval.c $(MY_SOURCES)
