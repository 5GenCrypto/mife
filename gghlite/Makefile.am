AUTOMAKE_OPTIONS = foreign
AM_CFLAGS = ${DEBUG_CFLAGS} -I$(top_srcdir) -I$(top_srcdir)/dgs -fopenmp

lib_LTLIBRARIES = libgghlite.la

include_HEADERS = gghlite.h \
                  gghlite-defs.h \
                  gghlite-internals.h \
                  misc.h

libgghlite_la_SOURCES = gghlite.c \
                        gghlite_pk.c \
                        misc.c \
                        lattice_reduction.c \
                        gghlite.h \
                        misc.h \
                        ggh-defs.h \
                        ggh-internals.h \
                        api.c
libgghlite_la_LIBADD = -lgmp -lmpfr \
                       -L$(top_builddir)/flint -lflint \
                       -L$(top_builddir)/oz -loz \
                       -L$(top_builddir)/dgs -ldgs \
                       -L$(top_builddir)/dgsl -ldgsl

libgghlite_la_LDFLAGS = -release 0.0.$(RELEASE) -no-undefined
libgghlite_la_INCLUDEDIR = $(includedir)/gghlite
