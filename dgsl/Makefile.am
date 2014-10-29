AUTOMAKE_OPTIONS = foreign
AM_CFLAGS = ${DEBUG_CFLAGS} -I$(top_srcdir) -I$(top_srcdir)/dgs -fopenmp

lib_LTLIBRARIES = libdgsl.la

libdgsl_la_SOURCES = dgsl.c gso.c
libdgsl_la_LDFLAGS = -release 0.0.$(RELEASE) -no-undefined
libdgsl_la_INCLUDEDIR = $(includedir)/dgsl

include_HEADERS = dgsl.h gso.h
