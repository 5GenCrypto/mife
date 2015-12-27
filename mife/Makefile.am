AUTOMAKE_OPTIONS = foreign -Wall

AM_CFLAGS  = ${DEBUG_CFLAGS} -I$(top_srcdir) -I$(top_srcdir)/dgs -fopenmp

lib_LTLIBRARIES=libmife.la

MY_LIB_LDFLAGS=-release 0.0.$(RELEASE) -no-undefined

libmife_la_SOURCES = mife.c
libmife_la_LDFLAGS = $(MY_LIB_LDFLAGS)
libmife_la_INCLUDEDOR = $(includedir)/mife
libmife_la_LIBADD = $(top_builddir)/gghlite/libgghlite.la \
	$(top_builddir)/dgs/libdgs.la \
	$(top_builddir)/oz/liboz.la \
        $(top_builddir)/dgsl/libdgsl.la \
	-L$(top_builddir)/flint -lflint \
        -lgmp -lmpfr -lm



include_HEADERS = mife.h
