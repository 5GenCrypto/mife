AUTOMAKE_OPTIONS = foreign -Wall -Wno-extra-portability

AM_CFLAGS  = ${DEBUG_CFLAGS} -I$(top_srcdir) -fopenmp

lib_LTLIBRARIES = liboz.la

MY_LIB_LDFLAGS=-release 0.0.$(RELEASE) -no-undefined

liboz_la_SOURCES = oz.c flint-addons.c util.c sqrt.c invert.c mul.c ntt.c norm.c
liboz_la_LDFLAGS = $(MY_LIB_LDFLAGS)
liboz_la_INCLUDEDOR = $(includedir)/oz

include_HEADERS = oz.h flint-addons.h sqrt.h invert.h mul.h
noinst_HEADERS = util.h
