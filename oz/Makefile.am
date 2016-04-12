AUTOMAKE_OPTIONS = foreign -Wall

AM_CFLAGS  = ${DEBUG_CFLAGS} -I$(top_srcdir) -fopenmp

lib_LTLIBRARIES=liboz.la

MY_LIB_LDFLAGS=-release 0.0.$(RELEASE) -no-undefined

liboz_la_SOURCES = oz.c flint-addons.c util.c sqrt.c invert.c mul.c ntt.c norm.c rem.c
liboz_la_LDFLAGS = $(MY_LIB_LDFLAGS)
liboz_la_INCLUDEDIR = $(includedir)/oz
liboz_la_LIBADD = -lflint -lgmp -lmpfr -lm -lgomp

pkgincludesubdir = $(includedir)/oz
pkgincludesub_HEADERS = oz.h flags.h flint-addons.h sqrt.h invert.h mul.h \
	norm.h rem.h ntt.h
noinst_HEADERS = util.h
