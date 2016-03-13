AUTOMAKE_OPTIONS = foreign -Wall
SUBDIRS = jsmn

AM_CFLAGS  = ${DEBUG_CFLAGS} -I$(top_srcdir) -I$(top_srcdir)/dgs -I$(top_srcdir)/flint -I$(top_srcdir)/mife/jsmn -fopenmp

lib_LTLIBRARIES=libmife.la

MY_LIB_LDFLAGS=-release 0.0.$(RELEASE) -no-undefined

libmife_la_SOURCES = mife.c mife_io.c mife_internals.c flint_raw_io.c parse.c mbp_glue.c util.c cmdline.c keygen.c
libmife_la_LDFLAGS = $(MY_LIB_LDFLAGS)
libmife_la_INCLUDEDOR = $(includedir)/mife
libmife_la_LIBADD = $(top_builddir)/gghlite/libgghlite.la \
	$(top_builddir)/dgs/libdgs.la \
	$(top_builddir)/oz/liboz.la \
        $(top_builddir)/dgsl/libdgsl.la \
	-L$(top_builddir)/flint -lflint \
        -lgmp -lmpfr -lm



include_HEADERS = mife.h mife_io.h mife_internals.h flint_raw_io.h mife_defs.h parse.h mbp_glue.h util.h cmdline.h mbp_types.h
