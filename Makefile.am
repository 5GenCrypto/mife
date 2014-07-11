AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I m4
SUBDIRS = flint dgs

AM_CFLAGS  = ${DEBUG_CFLAGS} -I$(srcdir)/dgs/ -I$(srcdir)/flint #-I$(srcdir)/flint/flint

lib_LTLIBRARIES = libgpv.la

libgpv_la_SOURCES = gpv/gpv.c
libgpv_la_LDFLAGS = -release 0.0.$(RELEASE) -no-undefined

include_HEADERS = gpv/gpv.h

bin_PROGRAMS = gghlite
gghlite_SOURCES = src/main.c src/gghlite.c src/misc.c src/gghlite.h src/misc.h
gghlite_LDADD = $(srcdir)/dgs/libdgs.la $(srcdir)/libgpv.la $(srcdir)/flint/libflint.a -lgmp -lmpfr -lm
gghlite_LDFLAGS = -no-install
