AUTOMAKE_OPTIONS = foreign
AM_CFLAGS = ${DEBUG_CFLAGS} -I$(top_srcdir) -I$(top_srcdir)/dgs -fopenmp

lib_LTLIBRARIES = libdgsl.la

libdgsl_la_SOURCES = dgsl.c gso.c

pkgincludesubdir = $(includedir)/dgsl
pkgincludesub_HEADERS = dgsl.h gso.h

libdgsl_la_LDFLAGS = -release 0.0.$(RELEASE) -no-undefined
libdgsl_la_INCLUDEDIR = $(includedir)/dgsl
libdgsl_la_LIBADD = -L$(top_builddir)/dgs -ldgs -lgmp  -lflint -L$(top_builddir)/oz -loz

# include_HEADERS = dgsl.h gso.h
