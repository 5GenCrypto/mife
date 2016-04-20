AUTOMAKE_OPTIONS = foreign -Wall

AM_CFLAGS  = ${DEBUG_CFLAGS} -I$(top_srcdir) -I$(top_srcdir)/dgs

lib_LTLIBRARIES=libmmap.la

MY_LIB_LDFLAGS=-release 0.0.$(RELEASE) -no-undefined

libmmap_la_SOURCES = mmap_clt.c mmap_gghlite.c
libmmap_la_LDFLAGS = $(MY_LIB_LDFLAGS)
libmmap_la_INCLUDEDIR = $(includedir)/oz
libmmap_la_LIBADD = -lgghlite -ldgsl

pkgincludesubdir = $(includedir)/mmap
pkgincludesub_HEADERS = mmap.h mmap_clt.h mmap_gghlite.h
